#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 14:54:18 2026

File to run the simulator of the HH model proposed by Doorn; in different
configurations. Inputs are passed via a dictionary, which determines the model
components, parameter settings, simulation settings, and which output is 
stored. 

@author: Guido
"""


# IMPORT LIBRARIES
import numpy as np
from brian2 import devices, BrianLogger
from brian2 import second
from brian2 import StateMonitor, SpikeMonitor, run, Network, collect
import brian2.codegen.cpp_prefs
import pickle
import random
import sys
sys.path.append("/home/Guido/UT_VU/py_code/HH_model/HH_utils")
from utils_model_eqs import get_equations_HH
from utils_helper_func import build_population, build_connection


# Skip C99 support check
brian2.codegen.cpp_prefs._compiler_supports_c99 = True
BrianLogger.suppress_hierarchy('brian2.devices')
BrianLogger.suppress_hierarchy('brian2.parsing')



def simulate_HH_net(simulation_dict):
    """
    Construct and simulate a Hodgkin–Huxley network in Brian2.
    
    This function:
    - Extracts model and simulation parameters from `simulation_dict`
    - Builds excitatory and/or inhibitory neuron populations
    - Configures synaptic mechanisms (AMPA, NMDA, GABA, STP, etc.)
    - Constructs network connectivity
    - Sets up requested StateMonitor and SpikeMonitor objects
    - Runs the simulation
    - Optionally saves recorded monitor data
    - Returns a dictionary containing all created monitors
    
    Input:
        simulation_dict : dict
            Dictionary containing:
                - runsettings_dict
                - network_dict
                - synapse_dict
                - recording_dict
                - excitatory_dict (optional)
                - inhibitory_dict (optional)
    
    Output:
        output_monitors : dict
            Dictionary mapping monitor names to Brian2 monitor objects.
    """
    runsettings_dict = simulation_dict["runsettings_dict"]
    network_dict = simulation_dict["network_dict"]    
    if "excitatory_dict" in simulation_dict:
        excitatory_dict = simulation_dict["excitatory_dict"]
    if "inhibitory_dict" in simulation_dict:
        inhibitory_dict = simulation_dict["inhibitory_dict"]
    synapse_dict = simulation_dict["synapse_dict"]
    recording_dict = simulation_dict["recording_dict"]
    
    #%% ### Set noise seed: ###
    
    noise_seed = runsettings_dict["noise_seed"]
    devices.device.seed(noise_seed)                           # set the seed for all the random number realisations        
   
    sim_name = runsettings_dict.get("sim_name", f"Simulation_{noise_seed}")
    print(f"\nConstructing Brian2 model for: {sim_name}")
    #%% ### Collect Biological Model Components: ###
    
    E_mechanisms = set()
    I_mechanisms = set()
    
    if "excitatory_dict" in simulation_dict and network_dict["n_e_neurons"] > 0: 
        include_E = True
        if excitatory_dict.get("include_adaptation", False): 
            E_mechanisms.add("AHP")
        if excitatory_dict.get("noise_sigma", 0) > 0: 
            E_mechanisms.add("neuronal_noise")       
        if excitatory_dict.get("include_ampa", False):
            E_mechanisms.add("ampa")        
        if excitatory_dict.get("include_nmda", False):
            E_mechanisms.add("nmda")
        if excitatory_dict.get("include_gaba", False):
            E_mechanisms.add("gaba")            
        if synapse_dict.get("include_depression_e", False):
            E_mechanisms.add("std")   
        if synapse_dict.get("include_facilitation_e", False):
            E_mechanisms.add("stf")      
        if synapse_dict.get("include_asynchr_e", False):
            E_mechanisms.add("asynchr")

    if "inhibitory_dict" in simulation_dict and network_dict["n_i_neurons"] > 0: 
        include_I = True        
        if inhibitory_dict.get("include_adaptation", False): 
            I_mechanisms.add("AHP")
        if inhibitory_dict.get("noise_sigma", 0) > 0: 
            I_mechanisms.add("neuronal_noise")
        if inhibitory_dict.get("include_ampa", False):
            I_mechanisms.add("ampa")        
        if inhibitory_dict.get("include_nmda", False):
            I_mechanisms.add("nmda")
        if inhibitory_dict.get("include_gaba", False):
            I_mechanisms.add("gaba")            
        if synapse_dict.get("include_depression_i", False):
            I_mechanisms.add("std")     
        if synapse_dict.get("include_facilitation_i", False):
            I_mechanisms.add("stf")     
        if synapse_dict.get("include_asynchr_i", False):
            I_mechanisms.add("asynchr")


    model_components_dict =  {
     "excitatory_neurons": include_E,
     "inhibitory_neurons": include_I,        
     "E_mechanisms": E_mechanisms,
     "I_mechanisms": I_mechanisms,
     } 
    
    equations_dict = get_equations_HH(model_components_dict)
    
    

    #%% ### Extract Equations and Construct Neuron Populations: ###    

    if include_E:
        P_E = build_population(
            pop_key="E",
            n_neurons=network_dict["n_e_neurons"],
            neuron_dict=excitatory_dict,
            equations_dict=equations_dict,
            synapse_dict=synapse_dict,
            model_mechanisms=E_mechanisms,
            include_inhibition=include_I
        )

    if include_I:
        P_I = build_population(
            pop_key="I",
            n_neurons=network_dict["n_i_neurons"],
            neuron_dict=inhibitory_dict,
            equations_dict=equations_dict,
            synapse_dict=synapse_dict,
            model_mechanisms=E_mechanisms,   # or I_mechanisms if you separate them
            include_inhibition=True    # By definiton
        )
            
    
    #%% ### Position of Neurons: ###
    
    placement_type = network_dict["neuron_position_type"]
    if placement_type == "grid":
        
        # position neurons on a grid
        grid_dist = network_dict["neuron_position_args"]["grid_distance"]      
        
        if include_E:         
            NlE = np.ceil(np.sqrt(network_dict["n_e_neurons"]))
            P_E.x = '(i % NlE) * grid_dist + 4 * umeter'
            P_E.y = '(i // NlE) * grid_dist + 2 * umeter'
        
        if include_I:         
            NlI = np.ceil(np.sqrt(network_dict["n_i_neurons"]))
            grid_dist_I = ((NlE-1)*grid_dist)/(NlI-1)       # also distribute inhibitory neurons homogeneously
            P_I.x = '(i % NlI) * grid_dist_I + 2 * umeter'
            P_I.y = '(i // NlI) * grid_dist_I+ 4 * umeter'


    #%% ### Connectivity and Synapses: ###
    
    if include_E:
        E_synapse_equation = equations_dict["E_synapse"]
        E_onpre_equation = equations_dict["E_onpre"]
        
        Conn_EE = build_connection(
            pre_group = P_E,
            post_group= P_E,
            conn_key= "e2e",                 # e.g. "e2e", "e2i", "i2e", "i2i"
            synapse_equation=E_synapse_equation,
            onpre_equation=E_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict= excitatory_dict,
            pre_model_mechanisms = E_mechanisms,
            post_model_mechanisms = E_mechanisms,
            network_dict=network_dict
        )
        
    if include_I: 
        I_synapse_equation = equations_dict["I_synapse"]
        I_onpre_equation = equations_dict["I_onpre"]
        Conn_II = build_connection(
            pre_group=P_I,
            post_group=P_I,
            conn_key="i2i",
            synapse_equation=I_synapse_equation,
            onpre_equation=I_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict=inhibitory_dict,
            pre_model_mechanisms=I_mechanisms,
            post_model_mechanisms=I_mechanisms,
            network_dict=network_dict
        )
    if include_E and include_I:
        Conn_EI = build_connection(
            pre_group=P_E,
            post_group=P_I,
            conn_key="e2i",
            synapse_equation=E_synapse_equation,
            onpre_equation=E_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict=inhibitory_dict,
            pre_model_mechanisms=E_mechanisms,
            post_model_mechanisms=I_mechanisms,
            network_dict=network_dict
        )
        Conn_IE = build_connection(
            pre_group=P_I,
            post_group=P_E,
            conn_key="i2e",
            synapse_equation=I_synapse_equation,
            onpre_equation=I_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict=excitatory_dict,
            pre_model_mechanisms=I_mechanisms,
            post_model_mechanisms=E_mechanisms,
            network_dict=network_dict
        )

    #%% ### Recording output_monitors: ###
    
    
    recordstring_E = []
    recordstring_I = []
    output_monitors = {}
    
    # Spikes and Voltage: 
    if include_E:
        if recording_dict.get("Spikes_e", False):
            output_monitors["spikes_E"] = SpikeMonitor(P_E, record=True)
        if recording_dict.get("Voltage_e", False):
            recordstring_E.append("V")
    if include_I:        
        if recording_dict.get("Spikes_i", False):
            output_monitors["spikes_I"] = SpikeMonitor(P_I, record=True)
        if recording_dict.get("Voltage_i", False):
            recordstring_I.append("V")
    
    # Synaptic currents
    if recording_dict.get("I_syn", False): 
        recordstring_E.append("I_syn") 
        recordstring_I.append("I_syn")
    
    currents_to_loop = ["I_ampa","I_gaba","I_nmda","I_AHP"]
    for current in currents_to_loop:
        if recording_dict.get(current, False):
            if current.split("_")[1] in E_mechanisms:     
                recordstring_E.append(current)
            if current.split("_")[1] in I_mechanisms:
                recordstring_I.append(current)
    
    # Set up neuron output_monitors
    dt2 = runsettings_dict["time_step"]

    
    
    if include_E:
        output_monitors["trace_E"] = StateMonitor(P_E, recordstring_E, record=True, dt=dt2)
        output_monitors["spikes_E"] = SpikeMonitor(P_E)
    if include_I:
        output_monitors["trace_I"] = StateMonitor(P_I, recordstring_I, record=True, dt=dt2)
        output_monitors["spikes_I"] = SpikeMonitor(P_I)            
    
    # Synapse output_monitors: 
    if recording_dict.get("depression", False): 
        if include_E:
            if len(Conn_EE) > 25:
                recordlist_E = random.sample(range(1, len(Conn_EE) - 1), 25)            
            else:
                recordlist_E = range(1, len(Conn_EE) - 1)
           
            if "std" in E_mechanisms:
                if "stf" in E_mechanisms and recording_dict.get("facilitation", False):
                    output_monitors["traceconn_EE"] = StateMonitor(Conn_EE, ['x_d', 'u_d'], record=recordlist_E)
                else:
                    output_monitors["traceconn_EE"] = StateMonitor(Conn_EE, ['x_d'], record=recordlist_E)
        if include_I:
            if len(Conn_II) > 25:
                recordlist_I = random.sample(range(1, len(Conn_II) - 1), 25)            
            else:
                recordlist_I = range(1, len(Conn_II) - 1)
            if "std" in I_mechanisms:
                if "stf" in I_mechanisms and recording_dict.get("facilitation", False):
                    output_monitors["traceconn_II"] = StateMonitor(Conn_II, ['x_d', 'u_d'], record=recordlist_I)
                else:
                    output_monitors["traceconn_II"] = StateMonitor(Conn_II, ['x_d'], record=recordlist_I)
            
        if include_E and include_I:
            if len(Conn_EI) > 25:
                recordlist_EI = random.sample(range(1, len(Conn_EI) - 1), 25)            
            else:
                recordlist_EI = range(1, len(Conn_EI) - 1)
                      
            if "std" in E_mechanisms:
                if "stf" in E_mechanisms and recording_dict.get("facilitation", False):
                    output_monitors["traceconn_EI"] = StateMonitor(Conn_EI, ['x_d', 'u_d'], record=recordlist_EI)
                else:
                    output_monitors["traceconn_EI"] = StateMonitor(Conn_EI, ['x_d'], record=recordlist_EI)
                    
            if len(Conn_IE) > 25:
                recordlist_IE = random.sample(range(1, len(Conn_IE) - 1), 25)            
            else:
                recordlist_IE = range(1, len(Conn_IE) - 1)
                         
            if "std" in I_mechanisms:
                if "stf" in I_mechanisms and recording_dict.get("facilitation", False):
                    output_monitors["traceconn_IE"] = StateMonitor(Conn_IE, ['x_d', 'u_d'], record=recordlist_IE)
                else:
                    output_monitors["traceconn_IE"] = StateMonitor(Conn_IE, ['x_d'], record=recordlist_IE)
    
    #%% ### Simulate! ###
    
    network_to_run = Network(collect())
    for key, monitor in output_monitors.items(): 
        network_to_run.add(monitor)
        
    network_to_run.run( runsettings_dict["sim_time"], report='text' )
    
                               
    if runsettings_dict.get("save_simulated_data", False):
        output_dir = runsettings_dict.get("output_dir", "/home/")                       
        sim_name = runsettings_dict.get("sim_name", f"Simulation_{noise_seed}")
        with open(output_dir +sim_name, 'wb') as f: 
            pickle.dump(output_monitors, f)
    
   
    #%% ### Construct Plots: NOT IMPLEMENTED YET ###
   
    save_figs = runsettings_dict.get("save_figs", None)                           # whether you want to save your figures. NOT IMplemented yET!
    sim_transient = runsettings_dict.get("sim_transient", 0 * second)             # transient time to discard when plotting/computing summstats
    
    
    return output_monitors
