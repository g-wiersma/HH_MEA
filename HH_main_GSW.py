#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:04:18 2026

Main file to construct and simulate the HH model.

@author: Guido
"""

from HH_utils.utils_simulator import simulate_HH_net
from brian2 import (second, ms, umeter, ufarad, mV, cm, msiemens, nsiemens, kHz, 
                pA, defaultclock)

#%% ### Define parameters and network configuration: ###


N_e = 80
N_i = 20
Maxdelay = 25 * ms
neuron_grid_dist = 405/7 * umeter 
Vmax = ((N_e+N_i + N_e+N_i) * neuron_grid_dist) / Maxdelay          # define in nicer way later; maximum velocity of propagation
print("\nCheck if Vmax is correctly implemented!!!")
area_E = 1000 * umeter ** 2                                         # area of a neuron
area_I = 1000 * umeter ** 2


simulation_dict = dict()

runsettings_dict = {
    "noise_seed": 111,
    "output_dir": None,
    "sim_name": "Test_config_x",
    "plot_figs": True,                                              # plots the output of variables in recording_dict
    "save_figs": False,
    "calculate_features": False,      
    "save_features": False,                                         # do / do not calculate MEA features.
    "save_simulated_data": False,
    "sim_time": 45 * second,
    "sim_transient": 5 * second,
    "time_step": defaultclock.dt,
    "electrode_grid_distance": 135 * umeter,
    "radius_detect_neurons": 135/2 * umeter # the radius in which an electrode measures neurons
    }   

network_dict = {    
    "n_e_neurons": N_e,
    "n_i_neurons": N_i,
    
    "p_connect_e2e": 0.3, #10/N_e,                                          # connection probability
    "p_connect_e2i": 0.2, #20/N_i,
    "p_connect_i2e": 0.2, #30/N_i,
    "p_connect_i2i": 0.1, #10/N_i,
    
    "S_connect_e2e": 0.2,                                          # weights of connections
    "S_connect_e2i": 0.2,
    "S_connect_i2e": 0.2,
    "S_connect_i2i": 0.2,
    
    "p_connect_type": "random",                                       # or random, small_world
    "p_connect_type_args": {"decay_constant": 150*umeter},
    "w_connect_type": "distributed",                                       # or small world
    "w_connect_type_args": {"sigma": 0.7, "lower_boundary": 0, "upper_boundary": 2}, 
    
    "distance_delays": True,
    "distance_delays_args": {"V_max": Vmax},

    "neuron_position_type": "random",                                    # method to place the neurons: random, or grid
    "neuron_position_args": {"neuron_grid_dist": neuron_grid_dist},         # average distance between E neurons (placed on a grid)
}

excitatory_dict = {
    # Surface / Capacitance
    "area": area_E,
    "Cm": (2.36 * ufarad * cm ** -2) * area_E,   

    # Reversal Potentials
    "E_leak": -58 * mV,                                             # Leak reversal potential,
    "E_K": -90 * mV,                                                # Potassium reversal potential,
    "E_Na": 55 * mV,                                                # Sodium reversal potential

    # (Maximal) Conductances
    "g_Na": (50 * msiemens * cm ** -2) * area_E,
    "g_K": (5 * msiemens * cm ** -2) * area_E,
    "g_leak": (0.1 * msiemens * cm ** -2) * area_E,

    # Excitability / firing threshold / noise standard deviation
    "I_sigma" : 19 * pA,                                            # difference in neuronal excitability, differs per neuron
    "V_threshold": -50 * mV,
    "noise_sigma": 5 * mV,                                          # how much noise fluctuates over time
    
    # refractory period post-spiking
    "refractory": 2 * ms,
    
    # Adaptation
    "include_adaptation": True,
    "E_AHP": -90 * mV,
    "g_AHP": 2.5 * nsiemens,
    "tau_Ca": 7 * second,
    "alpha_Ca": 0.00035*3,                                            # strength of the spike-frequency adaptation 

    # AMPA
    "include_ampa": True,
    "g_ampa": 3 * nsiemens,                                         # Maximal AMPA conductance
    "E_ampa": 0 * mV,                                               # AMPA reversal potential
    "tau_ampa": 2 * ms,                                             # AMPA decay time constant
    
    # NMDA
    "include_nmda": True,
    "g_nmda": 2 * nsiemens,                                         # maximal conductance of NMDA channels
    "E_nmda": 0 * mV,
    "taus_nmda": 100 * ms,                                          # decay time constant of nmda conductance
    "taux_nmda": 2 * ms,                                            # rise time constant of nmda conductance
    "alpha_nmda": 0.5 * kHz,
    
    # GABA    
    "include_gaba": True, 
    "g_gaba": 1 * nsiemens,                                         # Maximal GABA conductance
    "E_gaba": -62 * mV,                                             # GABA reversal potential
    "tau_gaba": 50 * ms,                                            # GABA decay time constant             
}

inhibitory_dict = {
    # Surface / Capacitance 
    "area": area_I,
    "Cm": (2.36 * ufarad * cm ** -2) * area_I,

    # Reversal Potentials
    "E_leak": -50 * mV, # Leak reversal potential,
    "E_K": -90 * mV,    # Potassium reversal potential,
    "E_Na": 55 * mV,    # Sodium reversal potential,

    # (Maximal) Conductances
    "g_Na": (50 * msiemens * cm ** -2) * area_I,
    "g_K": (5 * msiemens * cm ** -2) * area_I,
    "g_leak": (0.05 * msiemens * cm ** -2) * area_I,

    #  Excitability / firing threshold / noise standard deviation
    "I_sigma" : 19 * pA,
    "V_threshold": -42 * mV,
    "noise_sigma": 5 * mV,
    
    # refractory period post-spiking
    "refractory": 2 * ms,
    
    # Adaptation
    "include_adaptation": True, 
    "E_AHP": -90 * mV,
    "g_AHP": 2.5 * nsiemens,
    "tau_Ca": 7 * second,
    "alpha_Ca": 0.00035,                                            # strength of the spike-frequency adaptation 

    # AMPA
    "include_ampa": True,
    "g_ampa": 2 * nsiemens,                                         # Maximal AMPA conductance
    "E_ampa": 0 * mV,                                               # AMPA reversal potential
    "tau_ampa": 2 * ms,                                             # AMPA decay time constant
    
    # NMDA
    "include_nmda": True,
    "g_nmda": 2 * nsiemens,                                         # maximal conductance of NMDA channels
    "E_nmda": 0 * mV,
    "taus_nmda": 100 * ms,                                          # decay time constant of nmda conductance
    "taux_nmda": 2 * ms,                                            # rise time constant of nmda conductance
    "alpha_nmda": 0.5 * kHz,
    
    # GABA    
    "include_gaba": True, 
    "g_gaba": 1 * nsiemens,                                         # Maximal GABA conductance
    "E_gaba": -62 * mV,                                             # GABA reversal potential
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
    
    "strength_depression_e2e": 0.1,
    "strength_depression_e2i": 0.1,
    "strength_depression_i2e": 0.1,
    "strength_depression_i2i": 0.1,

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
    "Spikes_e": True, 
    "Spikes_i": True,
    "Voltage_e": True, 
    "Voltage_i": True, 
    "I_syn": True,
    "I_ampa": True, 
    "I_nmda": True, 
    "I_gaba": True, 
    "I_AHP": True,
    "depression": True, 
    "facilitation": False,
}

plot_dict = {
    "voltplot": True,                               # plots the membrane potential of a neuron in the network over time
    "rasterplot": True,                             # makes a raster plot of all neurons
    "electrodeplot": True,                          # plots the voltage measured by the 12 virtual electrodes
    "rasterlecplot": True,                          # raster plot of APs and network bursts detected at electrodes (requires electrodeplot=True)
    "STDplot": False,                                # plots the amount of short-term depression over time
    "synapseplot": False,                            # plots AMPA, NMDA and GABA currents
    "adaptationplot": False,                         # plots the afterhyperpolarization current
    "mechanismplot": True,                          # excitatory raster plot with EPSC, STD and adaptation
    "topologyplot": True,                           # shows neuron and electrode placement
    "onechannelplot": False,                         # plots the voltage signal of one MEA electrode
    "spikerateplot": True,                          # plots spike rate over time based on Brian spike detection
    "calciumplot": False,                            # plots calcium imaging intensity over time
    "patchplot": False,                              # shows simulated patch-clamp recording of a neuron
    "burstplot": True,                              # plot the detected bursts on the rasterplot
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
    


    
    
output_dict, feature_df = simulate_HH_net(simulation_dict)

# print("\nFeatures: \n")
# for col, val in feature_df.iloc[0].items():
#     print(col, val)


#%% Save simulation dict

import pickle

run_number = 6

runsettings_dict_save = {
    "noise_seed": 111,
    "output_dir": None,
    "sim_name": f"Test_config_dict_{str(run_number)}",
    "plot_figs": False,                                              # plots the output of variables in recording_dict
    "save_figs": False,
    "calculate_features": True,      
    "save_features": False,                                         # do / do not calculate MEA features.
    "save_simulated_data": False,
    "sim_time": 185 * second,
    "sim_transient": 5 * second,
    "time_step": defaultclock.dt,
    "electrode_grid_distance": 135 * umeter,
    "radius_detect_neurons": 135/2 * umeter # the radius in which an electrode measures neurons
    }   

simulation_dict["runsettings_dict"] = runsettings_dict_save
save_name = f"/home/Guido/UT_VU/py_code/Data_Maurice/{runsettings_dict_save['sim_name']}.pkl"
with open(save_name, 'wb') as f: 
    pickle.dump(simulation_dict, f)








#%% plot mean membrane voltage of population

plt.figure()
trace_E = output_dict["trace_E"]
plt.plot(trace_E.t, 1000*np.mean(trace_E.V, axis=0))
plt.xlabel("t (s)")
plt.ylabel("Mean membrane voltage E-neurons (mV)")
plt.show()

V_BP = bandpass_filter(1000*np.mean(trace_E.V, axis=0), 10000, lowcut=0.1, highcut=30, order=2)

plt.figure()
plt.plot(trace_E.t, V_BP)
plt.xlabel("t (s)")
plt.ylabel("Mean membrane voltage E-neurons BP (mV)")
plt.show()



