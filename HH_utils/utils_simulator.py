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
import timeit
import pickle
import random

from brian2 import devices, BrianLogger
from brian2 import StateMonitor, SpikeMonitor, Network, collect
from brian2 import second
import brian2.codegen.cpp_prefs

from HH_utils.utils_model_eqs import get_equations_HH
from HH_utils.utils_helper_func import build_population, build_connection
from HH_utils.utils_plot import get_plots_HH
from HH_utils.utils_analysis import get_electrode_info, get_features_HH

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
    analysis_dict = simulation_dict.get("analysis_dict", False)
    plot_dict = simulation_dict.get("plot_dict", False)
    
    
    
    #%% ### Set noise seed: ###
    
    noise_seed = runsettings_dict["noise_seed"]
    devices.device.seed(noise_seed)                           # set the seed for all the random number realisations        
   
    sim_name = runsettings_dict.get("sim_name", f"Simulation_{noise_seed}")
    print(f"\nConstructing Brian2 model for: {sim_name}")
    time0 = timeit.default_timer()
    
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
    else:
        include_E = False

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
    else: 
        include_I = False


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
            model_mechanisms=I_mechanisms,   
            include_inhibition=True    # By definiton
        )
            
    
    #%% ### Position of Neurons: ###
    
    placement_type = network_dict["neuron_position_type"]
    if placement_type == "grid":
        
        # position neurons on a grid
        grid_dist = network_dict["neuron_position_args"]["neuron_grid_dist"]      
        
        if include_E:                             
            NlE = np.ceil(np.sqrt(network_dict["n_e_neurons"]))          
            offsetE = (NlE - 1) * grid_dist / 2            
            P_E.x = '(i % NlE) * grid_dist - offsetE'
            P_E.y = '(i // NlE) * grid_dist - offsetE'
        
        if include_I:        
            NlI = np.ceil(np.sqrt(network_dict["n_i_neurons"]))
            grid_dist_I = ((NlE-1)*grid_dist)/(NlI-1)       # also distribute inhibitory neurons homogeneously
            offsetI = (NlI - 1) * grid_dist_I / 2    
            P_I.x = '(i % NlI) * grid_dist_I - offsetI'
            P_I.y = '(i // NlI) * grid_dist_I - offsetI'
    
    
    
    if placement_type == "random": 
        # The idea is to create a cirle around the electrodes, and randomly place
        # the neurons in this circle.
        np.random.seed(noise_seed)
        
        electrode_dist = runsettings_dict["electrode_grid_distance"]
        span_width_electrodes = 3 * electrode_dist
        radius_cirle_around_elec = np.sqrt((0.5*span_width_electrodes)**2 + (0.5*span_width_electrodes)**2)
        N_E = network_dict.get("n_e_neurons", 0)
        N_I = network_dict.get("n_i_neurons", 0)
        N_tot = N_E + N_I
        
        theta = 2 * np.pi * np.random.rand(N_tot)
        rad = radius_cirle_around_elec * np.sqrt(np.random.rand(N_tot))  # radius *sqrt(randn) to prevent 'overpopulation' in center
        theta_E = theta[0:N_E]
        theta_I = theta[N_E:]
        rad_E = rad[0:N_E]
        rad_I = rad[N_E:]
        
        if include_E:            
            P_E.x = rad_E*np.cos(theta_E)
            P_E.y = rad_E*np.sin(theta_E)       
        
        if include_I: 
            P_I.x = rad_I*np.cos(theta_I)
            P_I.y = rad_I*np.sin(theta_I)       
        


    #%% ### Connectivity and Synapses: ###
    
    if include_E:
        EE_synapse_equation = equations_dict["EE_synapse"]
        EE_onpre_equation = equations_dict["EE_onpre"]
        
        Conn_EE = build_connection(
            pre_group = P_E,
            post_group= P_E,
            conn_key= "e2e",                 # e.g. "e2e", "e2i", "i2e", "i2i"
            synapse_equation=EE_synapse_equation,
            onpre_equation=EE_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict= excitatory_dict,
            pre_model_mechanisms = E_mechanisms,
            post_model_mechanisms = E_mechanisms,
            network_dict=network_dict
        )
    
    if include_I: 
        II_synapse_equation = equations_dict["II_synapse"]
        II_onpre_equation = equations_dict["II_onpre"]
        Conn_II = build_connection(
            pre_group=P_I,
            post_group=P_I,
            conn_key="i2i",
            synapse_equation=II_synapse_equation,
            onpre_equation=II_onpre_equation,
            synapse_dict=synapse_dict,
            post_neuron_dict=inhibitory_dict,
            pre_model_mechanisms=I_mechanisms,
            post_model_mechanisms=I_mechanisms,
            network_dict=network_dict
        )
    if include_E and include_I:
        EI_synapse_equation = equations_dict["EI_synapse"]
        IE_synapse_equation = equations_dict["IE_synapse"]
        EI_onpre_equation = equations_dict["EI_onpre"]
        IE_onpre_equation = equations_dict["IE_onpre"]
        
        Conn_EI = build_connection(
            pre_group=P_E,
            post_group=P_I,
            conn_key="e2i",
            synapse_equation=EI_synapse_equation,
            onpre_equation=EI_onpre_equation,
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
            synapse_equation=IE_synapse_equation,
            onpre_equation=IE_onpre_equation,
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
        if include_E: 
            recordstring_E.append("I_syn") 
        if include_I:
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
    if include_I:
        output_monitors["trace_I"] = StateMonitor(P_I, recordstring_I, record=True, dt=dt2)
              
    
    # Synapse output_monitors: 
    if recording_dict.get("depression", False): 
        if include_E:
            if len(Conn_EE) > 25:       # select 25 random E-neurons to study.
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
    
    time1 = timeit.default_timer()
    print(f"\nNetwork constructed in {(time1-time0):.0f} s. \nInitializing simulation...\n" )
    network_to_run.run( runsettings_dict["sim_time"], report='text' )
    
                               
    if runsettings_dict.get("save_simulated_data", False):
        output_dir = runsettings_dict.get("output_dir", "/home/")                       
        sim_name = runsettings_dict.get("sim_name", f"Simulation_{noise_seed}")
        with open(output_dir +sim_name + "_Monitors.pkl", 'wb') as f: 
            pickle.dump(output_monitors, f)
    
   
    #%% ### Calculate features: ###
    
     
    if include_E: 
        N_E = network_dict["n_e_neurons"]
    else: 
        N_E = None
        P_E = None
    if include_I:
        N_I = network_dict["n_i_neurons"]
    else: 
        N_I = None
        P_I = None
        
    if plot_dict: 
        electrodeplot= plot_dict.get("electrodeplot", False)
        rasterlecplot = plot_dict.get("rasterlecplot", False)
        topologyplot = plot_dict.get("topologyplot", False)
        onechannelplot = plot_dict.get("onechannelplot", False)
   
    if (electrodeplot or rasterlecplot or topologyplot or onechannelplot) or (
        analysis_dict is not False and runsettings_dict.get("calculate_features", True) ):
        
        electrode_info_dict = get_electrode_info(
                runsettings_dict, output_monitors, excitatory_pop=P_E, 
                inhibitory_pop=P_I, N_E=N_E, N_I=N_I)
    else:
        electrode_info_dict = None # variable is needed, but can be empty
    
    if analysis_dict is not False:
        if runsettings_dict.get("calculate_features", True):
            time2 = timeit.default_timer()
            print("\nCalculating summary features...")
            APs = electrode_info_dict["APs"]
            rec_time = runsettings_dict["sim_time"] / second
            transient = runsettings_dict.get("transient", 0 * second) / second
            dt2 = runsettings_dict["time_step"]                
            fs = 1 / (dt2 / second)
            
            NBs, feature_df = get_features_HH(analysis_dict, APs, rec_time, transient, fs, return_NBs=True)
            electrode_info_dict["NBs"] = NBs
            if runsettings_dict.get("save_features", False): 
                output_dir = runsettings_dict.get("output_dir", "/home/")                       
                sim_name = runsettings_dict.get("sim_name", f"Simulation_{noise_seed}")
                with open(output_dir +sim_name + "_features.pkl", 'wb') as f: 
                    pickle.dump(feature_df, f)
            time3 = timeit.default_timer()
            print(f"\nSummary features calculated in {(time3-time2):.0f}s" )
        else: 
            feature_df = None
            print("\nNo data-describing features calculated ('calculate_features' in the runsettings dictionary is set to 'False')")
    else: 
        feature_df = None
        print("\nNo data-describing features calculated (no analysis dictionary given)")
        
    #%% ### Construct Plots: ###
    
    
    if runsettings_dict.get("plot_figs", True):
        if plot_dict is not False:
            get_plots_HH(simulation_dict, output_monitors, model_components_dict, 
                     electrode_info_dict, P_E, P_I)
        else: 
            print("\nNo figures generated ('plot_dict' not found in the simulation dictionary)")
    else: 
        print("\nNo figures generated ('plot_figs' in the runsettings is set to 'False')")
                
    time4= timeit.default_timer()
    print(f"\nTotal runtime: {(time4-time0):.0f} s" )    
    
    return output_monitors, feature_df
