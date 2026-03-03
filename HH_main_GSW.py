#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 15:04:18 2026

Main file to construct and simulate the HH model.

@author: Guido
"""

import sys
sys.path.append("/home/Guido/UT_VU/py_code/HH_model/HH_utils")
from utils_simulator import simulate_HH_net

import numpy as np
from brian2 import (second, ms, umeter, ufarad, mV, cm, msiemens, nsiemens, kHz, 
                pA, defaultclock)

#%% ### Define parameters and network configuration: ###


N_e = 80
N_i = 20
Maxdelay = 25 * ms
grid_dist = 45 * umeter 
Vmax = (np.sqrt(N_e+N_i + N_e+N_i) * grid_dist) / Maxdelay          # define in nicer way later; maximum velocity of propagation
print("Check if Vmax is correctly implemented!!!\n")
area_E = 1000 * umeter ** 2
area_I = 1000 * umeter ** 2


simulation_dict = dict()

runsettings_dict = {
    "noise_seed": 101,
    "output_dir": None,
    "sim_name": "TestTestTest",
    "plot_figs": True,                                              # plots the output of variables in recording_dict
    "save_figs": False,
    "save_simulated_data": False,
    "sim_time": 50 * second,
    "sim_transient": 5 * second,
    "time_step": defaultclock.dt
    }   

network_dict = {    
    "n_e_neurons": N_e,
    "n_i_neurons": N_i,
    
    "p_connect_e2e": 0.01,                                          # connection probability
    "p_connect_e2i": 0.01,
    "p_connect_i2e": 0.01,
    "p_connect_i2i": 0.1,
    
    "S_connect_e2e": 0.01,                                          # weights of connections
    "S_connect_e2i": 0.01,
    "S_connect_i2e": 0.01,
    "S_connect_i2i": 0.01,
    
    "connect_type": "random", 
    "connect_type_args": None,
    "S_connect_distributed": True, 
    "S_connect_distributed_args": {"sigma": 0.7, "lower_boundary": 0, "upper_boundary": 2}, 
    
    "distance_delays": True,
    "distance_delays_args": {"V_max": Vmax},

    "neuron_position_type": "grid",                                
    "neuron_position_args": {"grid_distance": 45 * umeter},         # average distance between E neurons (placed on a grid)
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
    "I_sigma" : 19 * pA,
    "V_threshold": -50 * mV,
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
    "include_nmda": False,
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
    "include_nmda": False,
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
    "sat_asynchr_e": 0.5 / ms,                                      # saturation level of asynchronous release
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

simulation_dict = {
    "runsettings_dict": runsettings_dict,
    "network_dict": network_dict,
    "excitatory_dict": excitatory_dict,
    "inhibitory_dict": inhibitory_dict,
    "synapse_dict": synapse_dict,
    "recording_dict": recording_dict,
}


output_dict = simulate_HH_net(simulation_dict)





















