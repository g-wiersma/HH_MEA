#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:36:23 2026

@author: Guido
"""
import brian2 as b2

import os
import sys
import pickle
from HH_utils.utils_simulator import simulate_HH_net
from brian2 import (second, ms, umeter, ufarad, pfarad, mV, cm, msiemens, nsiemens, kHz, 
                pA, defaultclock)
import numpy as np
import shutil
import h5py

simulation_num = int(float(sys.argv[1]))


base_cache_dir = os.getcwd() # Update this with a real directory
unique_cache_dir = os.path.join(base_cache_dir, str(simulation_num))    
os.makedirs(unique_cache_dir, exist_ok=True)
b2.prefs.codegen.runtime.cython.cache_dir = unique_cache_dir  

print("cache directory created succesfully")

##############################################################################

#%% select right simulation: 
    
missing_sims_indices = np.load("/home/cnph-slurm/Guido/missing_sims.npy")

simulation_num = int(missing_sims_indices[simulation_num])


#%% ### Define parameters and network configuration: ###

try: 
    
    # load parameters to simulate
    PARAM_FILE = "/home/cnph-slurm/Guido/HH_EI_pilot_params.pkl"
    with open(PARAM_FILE, 'rb') as f: 
        df_parameters_all = pickle.load(f)
    
    print("Prior parameters loaded succesfully")    
    df_parameters_run = df_parameters_all.iloc[simulation_num, :]
    print("Parameters of this run selected succesfully")    
    
    ### Standard parameters: ###
    
    Maxdelay = 25 * ms
    electrode_grid_dist = 135 * umeter
    electrode_field_r = 0.35*electrode_grid_dist
    radius_cirle_around_elec = electrode_field_r*np.sqrt(100/3) # makes sure that 36/100 of the electrode field captures the total neuron area # np.sqrt(electrode_field_r**2 + (3*electrode_field_r)**2) + electrode_field_r  
    Vmax = radius_cirle_around_elec*2 / Maxdelay          # maybe define in nicer way later; maximum velocity of propagation
    
    area_E = 1000 * umeter ** 2                                         # area of a neuron
    area_I = 1000 * umeter ** 2
    N_neurons_tot = 100
    
    ###############################################################################
    
    ### Model Configurations: ###
    include_I_neurons = df_parameters_run["include_I_neurons"]
    include_AHP = df_parameters_run["include_AHP"]
    
    
    ### Parameter Configurations: ###
    if include_AHP:
        g_AHP_E = (df_parameters_run.get("g_AHP", 0) * msiemens * cm ** -2) * area_E # [0.2, 10]
        g_AHP_I = g_AHP_E
    else:
        g_AHP_E = g_AHP_I = 0
    
    if include_I_neurons: 
        p_I_neuron =  p_I_neuron = df_parameters_run.get("p_I_neuron", 0)  
        S_connect_i2i = df_parameters_run.get("S_connect_i2i", 0) # [0.1, 10]
        noise_sigma_I = df_parameters_run.get("noise_sigma_I", 0) * mV # [1.5, 8]
        if p_I_neuron != 1:
            
            S_connect_e2e = df_parameters_run.get("S_connect_e2e", 0) # [0.05, 10]
            S_connect_e2i = df_parameters_run.get("S_connect_e2i", 0) # [0.1, 10]
            S_connect_i2e = df_parameters_run.get("S_connect_i2e", 0)# [0.1, 10]            
            
            noise_sigma_E = df_parameters_run.get("noise_sigma_E", 0) * mV # [1.5, 8]                                          # how much noise fluctuates over time
            
        else: # only I
            S_connect_e2e = 0
            S_connect_e2i = 0
            S_connect_i2e = 0
            
            noise_sigma_E = 0 * mV
    else: # Only E
        S_connect_e2e = df_parameters_run.get("S_connect_e2e", 0) # [0.05, 10]
        noise_sigma_E = df_parameters_run.get("noise_sigma_E", 0) * mV # [1.5, 8]          
       
        S_connect_i2i = 0
        S_connect_e2i = 0
        S_connect_i2e = 0
        
        noise_sigma_I = 0 * mV
    
    
    ###############################################################################
    
    ### Set p_connect based on %I neurons: ###
    
    if include_I_neurons: 
        N_i = int(np.round(N_neurons_tot * p_I_neuron))
        N_e = int(N_neurons_tot - N_i)
        if p_I_neuron == 0.2: 
            p_e2e = 0.2
            p_e2i = p_e2e
            p_i2e = 0.32
            p_i2i = p_i2e
        elif p_I_neuron == 0.4:
            p_e2e = 0.25
            p_e2i = p_e2e
            p_i2e = 0.2
            p_i2i = p_i2e 
        elif p_I_neuron == 0.6:
            p_e2e = 0.32
            p_e2i = p_e2e
            p_i2e = 0.18
            p_i2i = p_i2e 
        elif p_I_neuron == 0.8:
            p_e2e = 0.53
            p_e2i = p_e2e
            p_i2e = 0.15
            p_i2i = p_i2e 
        elif p_I_neuron == 1: 
            p_e2e = 0
            p_e2i = 0
            p_i2e = 0
            p_i2i = 0.11
    else:
        N_e = N_neurons_tot
        N_i = 0
        p_e2e = 0.18
        p_e2i = 0
        p_i2e = 0
        p_i2i = 0
    
    
    
    ###############################################################################
    
    simulation_dict = dict()
    
    runsettings_dict = {
        "noise_seed": simulation_num,
        "output_dir": None,                                             # to save/store results and/or figures
        "sim_name": f"HH_MEA_pilot_{simulation_num}",
        "plot_figs": False,                                              # plots the output of variables in recording_dict
        "save_figs": False,
        "calculate_features": True,      
        "save_features": False,                                         
        "save_simulated_data": False,
        "sim_time": 180 * second,
        "sim_transient": 5 * second,                                    # first seconds to discard
        "time_step": defaultclock.dt,
        "electrode_grid_distance": electrode_grid_dist,
        "radius_detect_neurons": electrode_field_r,                         # the radius in which an electrode measures neurons
        "radius_MEA_well": radius_cirle_around_elec                 # the radius of the well in which neurons are placed
        }   
    
    network_dict = {    
        "n_e_neurons": N_e,
        "n_i_neurons": N_i,
        
        "p_connect_e2e": p_e2e * 100/N_neurons_tot, #10/N_e,                                  # connection probability
        "p_connect_e2i": p_e2i * 100/N_neurons_tot, #20/N_i,
        "p_connect_i2e": p_i2e * 100/N_neurons_tot, #30/N_i,
        "p_connect_i2i": p_i2i * 100/N_neurons_tot, #10/N_i,
        
        "S_connect_e2e": S_connect_e2e,                                          # weights of connections
        "S_connect_e2i": S_connect_e2i,
        "S_connect_i2e": S_connect_i2e,
        "S_connect_i2i": S_connect_i2i,
        
        "p_connect_type": "random",                                       # random or small_world
        "p_connect_type_args": {"decay_constant": 150*umeter},
        "w_connect_type": "distributed",                                       # distributed or small world
        "w_connect_type_args": {"sigma": 0.7, "lower_boundary": 0, "upper_boundary": 2}, 
        
        "distance_delays": True,
        "distance_delays_args": {"V_max": Vmax},
    
        "neuron_position_type": "random",                                    # method to place the neurons over the electrodes: random, or grid
        "neuron_position_args": {"neuron_grid_dist": None},         # average distance between E neurons (placed on a grid)
    }
    
    excitatory_dict = {
        # Surface / Capacitance
        "area": area_E,
        "Cm": 74 * pfarad,   
    
        # Reversal Potentials
        "E_leak": -58 * mV,                                             # Leak reversal potential,
        "E_K": -90 * mV,                                                # Potassium reversal potential,
        "E_Na": 55 * mV,                                                # Sodium reversal potential
    
        # (Maximal) Conductances
        "g_Na": (50 * msiemens * cm ** -2) * area_E,
        "g_K": (5 * msiemens * cm ** -2) * area_E,
        "g_leak": (1/0.627) * nsiemens,
    
        # Excitability / firing threshold / noise standard deviation
        "I_sigma" : 19 * pA,                                            # difference in neuronal excitability, differs per neuron
        "V_threshold": -50 * mV,
        "noise_sigma": noise_sigma_E,                                          # how much noise fluctuates over time
        
        # refractory period post-spiking
        "refractory": 2 * ms,
        
        # Adaptation
        "include_adaptation": include_AHP,
        "E_AHP": -90 * mV,
        "g_AHP": g_AHP_E,
        "tau_Ca": 4 * second,
        "alpha_Ca": 0.00035,                                            # strength of the spike-frequency adaptation 
    
        # AMPA
        "include_ampa": True,
        "g_ampa": 1 * nsiemens,                                         # Maximal AMPA conductance
        "E_ampa": 0 * mV,                                               # AMPA reversal potential
        "tau_ampa": 2 * ms,                                             # AMPA decay time constant
        
        # NMDA
        "include_nmda": False,
        "g_nmda": None,                                         # maximal conductance of NMDA channels
        "E_nmda": None,
        "taus_nmda": None,                                          # decay time constant of nmda conductance
        "taux_nmda": None,                                            # rise time constant of nmda conductance
        "alpha_nmda": None,
        
        # GABA    
        "include_gaba": include_I_neurons, 
        "g_gaba": 1 * nsiemens,                                         # Maximal GABA conductance
        "E_gaba": -61.85 * mV,                                             # GABA reversal potential
        "tau_gaba": 50 * ms,                                            # GABA decay time constant             
    }
    
    inhibitory_dict = {
        # Surface / Capacitance 
        "area": area_I,
        "Cm": 63 * pfarad,
    
        # Reversal Potentials
        "E_leak": -50 * mV, # Leak reversal potential,
        "E_K": -90 * mV,    # Potassium reversal potential,
        "E_Na": 55 * mV,    # Sodium reversal potential,
    
        # (Maximal) Conductances
        "g_Na": (50 * msiemens * cm ** -2) * area_I,
        "g_K": (5 * msiemens * cm ** -2) * area_I,
        "g_leak": (1/0.739) * nsiemens,
    
        #  Excitability / firing threshold / noise standard deviation
        "I_sigma" : 19 * pA,
        "V_threshold": -42 * mV,
        "noise_sigma": noise_sigma_I,
        
        # refractory period post-spiking
        "refractory": 2 * ms,
        
        # Adaptation
        "include_adaptation": include_AHP, 
        "E_AHP": -90 * mV,
        "g_AHP": g_AHP_I,
        "tau_Ca": 4 * second,
        "alpha_Ca": 0.00035,                                            # strength of the spike-frequency adaptation 
    
        # AMPA
        "include_ampa": True,
        "g_ampa": 2 * nsiemens,                                         # Maximal AMPA conductance
        "E_ampa": 0 * mV,                                               # AMPA reversal potential
        "tau_ampa": 2 * ms,                                             # AMPA decay time constant
        
        # NMDA
        "include_nmda": False,
        "g_nmda": None,                                         # maximal conductance of NMDA channels
        "E_nmda": None,
        "taus_nmda": None,                                          # decay time constant of nmda conductance
        "taux_nmda": None,                                            # rise time constant of nmda conductance
        "alpha_nmda": None,
        
        # GABA    
        "include_gaba": include_I_neurons, 
        "g_gaba": 1 * nsiemens,                                         # Maximal GABA conductance
        "E_gaba": -61.85 * mV,                                             # GABA reversal potential
        "tau_gaba": 50 * ms,                                            # GABA decay time constant             
    }
    
    synapse_dict = {
    
        # Short-term synaptic depression
        "include_depression_e": True,
        "include_depression_i": True,
        
        "tau_depression_e2e": 500 * ms,
        "tau_depression_e2i": 500 * ms,
        "tau_depression_i2e": 500 * ms,
        "tau_depression_i2i": 500 * ms,
        
        "strength_depression_e2e": 0.025,
        "strength_depression_e2i": 0.025,
        "strength_depression_i2e": 0.025,
        "strength_depression_i2i": 0.025,
    
        # Short-term synaptic facilitation (time constants)
        "include_facilitation_e": False,
        "include_facilitation_i": False,
        
        "tau_facilitation_e2e": 10 * ms,
        "tau_facilitation_e2i": 10 * ms,
        "tau_facilitation_i2e": 10 * ms,
        "tau_facilitation_i2i": 10 * ms,
        
        # Asynchronous release E
        "include_asynchr_e": False,
        "tau_asynchr_e": 700 * ms,                                      # time constant of asynchronous release
        "strength_asynchr_e": 0.1,                                      # strength of asynchronous release
        "sat_asynchr_e": 0.1 / ms,                                      # saturation level of asynchronous release
        "NT_per_vesicle_e": 5,                                          # number of neurotransmitters in a vesicle
        
        # Asynchronous release I
        "include_asynchr_i": False,
        "tau_asynchr_i": 700 * ms,                                      # time constant of asynchronous release
        "strength_asynchr_i": 0.1,                                      # strength of asynchronous release
        "sat_asynchr_i": 0.5 / ms,                                      # saturation level of asynchronous release
        "NT_per_vesicle_i": 5,                                          # number of neurotransmitters in a vesicle    
    }
    
    recording_dict = {
        # Measures directly from the neurons:    
        "Spikes_e": True, 
        "Spikes_i": True,
        "Voltage_e": True, 
        "Voltage_i": True, 
        "I_syn": False,
        "I_ampa": False, 
        "I_nmda": False, 
        "I_gaba": False, 
        "I_AHP": False,
        "depression": False, 
        "facilitation": False,
        
        # Voltage and spiketimes measured by the electrodes: 
        "electrode_info": True,
    }
    
    plot_dict = {
        "neuron_voltplot": False,                         # plots the membrane potential of a neuron in the network over time
        "neuron_rasterplot": True,                       # makes a raster plot of all neurons
        "electrode_voltplot": False,                      # plots the voltage measured by the 12 virtual electrodes
        "electrode_rasterplot": False,                    # raster plot of APs and network bursts detected at electrodes (requires electrode_voltplot=True)
        "STDplot": False,                                # plots the amount of short-term depression over time
        "synapseplot": False,                            # plots AMPA, NMDA and GABA currents
        "adaptationplot": False,                         # plots the afterhyperpolarization current
        "summary_mechanismplot": False,                   # post-synaptic currents, firing rate, STD and AHP in one plot
        "topologyplot": False,                            # shows neuron and electrode placement
        "one_electrode_voltplot": False,                 # plots the voltage signal of one MEA electrode
        "neuron_spikerateplot": False,                    # plots spike rate over time based on Brian spike detection (not detected, but directly from the neurons)
        "calciumplot": False,                            # plots calcium imaging intensity over time
        "patchplot": False,                              # shows simulated patch-clamp recording of a neuron
    }
    
    analysis_dict = {
        "calc_spike_features": True, 
        "calc_burst_features": True, 
        "analysis_args": {      
            "time_bin": 25 * ms,                        # size of time bins to compute firing rate in seconds
            "smoothing": {"width": 7, "sigma": 3.0},    # parameters to smooth the firing rate with gaussian kernel
            "burst_detection": {                        # burst detection parameters
                "Startth_factor": 0.25,                 # factor * max(firingrate) = threshold to start the NB detection
                "Stopth_factor": 1/30,                  # factor * max(firingrate) = threshold to stop NB detection
                "time_th": 50 * ms,                     # minimal NB duration in seconds
                "Eth_factor": 0.5,                      # What fraction of electrodes needs to be active within a burst
                "Act_thres": 0.02,                      # minimal mean firingrate and electrode needs to have to be considered active
                "FBth_factor": 1 / 20,                  # factor * max(firingrate) = threshold to detect fragment peaks in FR
                "FBprom_factor": 1 / 30,                # factor * max(firingrate) = minimal prominence of fragments in FR
                "detect_mega_network_bursts": True,     # in this case, network bursts are grouped in megaburts where possible
                "megaburst": 
                    {"lowerbound": 25e-3, "upperbound": 1.5, "resetvalue": 50e-3},   #bounds for allowed INBI threshold  
            }
        }
    }
    
    
    # Collect all input settings into one dictionary:
    simulation_dict = {
        "runsettings_dict": runsettings_dict,
        "network_dict": network_dict,
        "excitatory_dict": excitatory_dict,
        "inhibitory_dict": inhibitory_dict,
        "synapse_dict": synapse_dict,
        "recording_dict": recording_dict,
        "plot_dict": plot_dict,
        "analysis_dict": analysis_dict, 
    }
        
        
        
    output_dict, electrode_info_dict, feature_df = simulate_HH_net(simulation_dict)
    
    APs = electrode_info_dict["APs"]
    
    fs = 1 / (defaultclock.dt / second)
    
    
    ### Save as HDF5 ###
    
    result_dir = sys.argv[2]
    feat_dir = sys.argv[3]
    
    hdf5_AP_filename = os.path.join(result_dir, f"APs_{simulation_num}.h5")
    feat_filename = os.path.join(feat_dir, f"features_{simulation_num}.pkl")
    
    
    with h5py.File(hdf5_AP_filename, "w") as h5f:
        h5f.create_dataset("APs", data=APs, compression="gzip")    
        h5f.create_dataset("noise_seed", data=simulation_num)
        h5f.create_dataset("fs", data=fs)
    
    with open(feat_filename, 'wb') as f: 
        pickle.dump(feature_df, f)

        
                
        
    shutil.rmtree(unique_cache_dir)
        
    
    
except Exception as e:
    
    print("-----")
    print(f"Simulation {simulation_num}: failed")
    print("-----")
    print(e)
    print("-----")
    
    if os.path.exists(unique_cache_dir):
        shutil.rmtree(unique_cache_dir)
    