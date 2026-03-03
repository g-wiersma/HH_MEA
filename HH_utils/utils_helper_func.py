#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 15:23:03 2026



@author: Guido
"""

from brian2 import NeuronGroup, Synapses, second, pA

###############################################################################

#%% ### Helper Functions: ###


def build_population(
    pop_key,                  # "E" or "I"
    n_neurons,
    neuron_dict,
    equations_dict,
    synapse_dict,
    model_mechanisms,
    include_inhibition=False
):
    """
    Construct a Brian2 NeuronGroup population (excitatory or inhibitory)
    and assign all membrane and synaptic parameters.

    Input: 
        pop_key : str
            Population identifier ("E" or "I"), used to retrieve equations.
        n_neurons : int
            Number of neurons in the population.
        neuron_dict : dict
            Dictionary containing intrinsic neuron parameters
            (e.g., Cm, E_leak, g_Na, etc.).
        equations_dict : dict
            Dictionary containing neuron, synapse, and on_pre equations
            with keys formatted as "{pop_key}_neuron", "{pop_key}_synapse",
            "{pop_key}_onpre".
        synapse_dict : dict
            Dictionary containing synaptic conductances, reversal potentials,
            and time constants (AMPA, NMDA, GABA).
        model_mechanisms : set or list
            Enabled synaptic mechanisms (e.g., {"nmda"}).
        include_inhibition : bool, optional
            Whether inhibitory conductances (GABA) should be added
            to the population.

    Output:    
        P : NeuronGroup
            Constructed Brian2 neuron population.

    """
    # Extract equations
    neuron_eq = equations_dict.get(f"{pop_key}_neuron")
    
  
    refractory = neuron_dict.get("refractory", 0 * second)
    std_I      = neuron_dict.get("I_sigma", 0 * pA)
    
    if "AHP" in model_mechanisms: 
        alpha_Ca = neuron_dict["alpha_Ca"]        
        P = NeuronGroup(
            n_neurons,
            model=neuron_eq,
            threshold='V>0*mV',
            reset='Ca += alpha_Ca',
            refractory=refractory,
            method='exponential_euler'
        )
    else:
        P = NeuronGroup(
            n_neurons,
            model=neuron_eq,
            threshold='V>0*mV',            
            refractory=refractory,
            method='exponential_euler'
        )
        
    # Initial conditions
    P.V = neuron_dict["E_leak"]
    P.I = '(rand() - 0.5) * std_I'
    
    # Membrane parameters
    P.Cm  = neuron_dict["Cm"]
    P.VT  = neuron_dict["V_threshold"]
    P.El  = neuron_dict["E_leak"]
    P.EK  = neuron_dict["E_K"]
    P.ENa = neuron_dict["E_Na"]
    
    P.g_na = neuron_dict["g_Na"]
    P.g_kd = neuron_dict["g_K"]
    P.gl   = neuron_dict["g_leak"]
    
    # AMPA
    P.g_ampa  = neuron_dict["g_ampa"]
    P.E_ampa  = neuron_dict["E_ampa"]
    P.tau_ampa = neuron_dict["tau_ampa"]
    
    # Neuronal Noise (optional)
    if "AHP" in model_mechanisms: 
        P.tau_Ca = neuron_dict["tau_Ca"]
        P.g_AHP = neuron_dict["g_AHP"]
       
    if "neuronal_noise" in model_mechanisms: 
        P.noise_sigma = neuron_dict["noise_sigma"]
    else: 
        P.noise = 0
        
    # NMDA (optional)
    if "nmda" in model_mechanisms:
        P.E_nmda = neuron_dict["E_nmda"]
        P.g_nmda = neuron_dict["g_nmda"]
    
    # GABA (optional)
    if include_inhibition:
        P.E_gaba  = neuron_dict["E_gaba"]
        P.g_gaba  = neuron_dict["g_gaba"]
        P.tau_gaba = neuron_dict["tau_gaba"]
    
    return P


def build_connection(
    pre_group,
    post_group,
    conn_key,                 # e.g. "e2e", "e2i", "i2e", "i2i"
    synapse_equation,
    onpre_equation,
    synapse_dict,
    post_neuron_dict,
    pre_model_mechanisms,
    post_model_mechanisms,
    network_dict
):
    """
    Construct a Brian2 Synapses object between two neuron populations
    and assign all mechanism-dependent parameters.

    Input:    
        pre_group : NeuronGroup
            Presynaptic population.
        post_group : NeuronGroup
            Postsynaptic population.
        conn_key : str
            Connection identifier (e.g., "e2e", "e2i", "i2e", "i2i").
            Used to retrieve connection-specific parameters.
        synapse_equation : str
            Synapse model equations.
        onpre_equation : str
            on_pre event equations.
        synapse_dict : dict
            Dictionary containing synaptic parameters.
        model_mechanisms : set or list
            Enabled synaptic mechanisms (e.g., {"std", "stf", "nmda"}).
        network_dict : dict
            Network-level configuration (e.g., distance delays).

    Output:
        Conn : Synapses
            Constructed and configured Synapses object.
    """
    connect_type = network_dict["connect_type"] # or clustered, or small world, or hub-based, diectional, 
    connect_args = network_dict.get("connect_type_args", False)

    S_value= network_dict.get(f"S_connect_{conn_key}", False)   
    
    Conn = Synapses(
        pre_group,
        post_group,
        model=synapse_equation,
        on_pre=onpre_equation,
        method="euler"
    )

    # Connectivity pattern
    if connect_type == "random":
        Conn.connect(p=network_dict[f"p_connect_{conn_key}"]   )
        
        if S_value is not None:
            Conn.S = S_value

    # Short-term depression
    if "std" in pre_model_mechanisms:
        Conn.tau_d = synapse_dict[f"tau_depression_{conn_key}"]
        Conn.U     = synapse_dict[f"strength_depression_{conn_key}"]

        # Short-term facilitation (requires std in your logic)
        if "stf" in pre_model_mechanisms:
            Conn.tau_f = synapse_dict[f"tau_facilitation_{conn_key}"]

    # NMDA
    if "nmda" in post_model_mechanisms and conn_key in ['e2e', 'e2i']:
        Conn.taus_nmda  = post_neuron_dict["taus_nmda"]
        Conn.alpha_nmda = post_neuron_dict["alpha_nmda"]
        Conn.taux_nmda  = post_neuron_dict["taux_nmda"]

    # Asynchronous release
    if "asynchr" in pre_model_mechanisms:
        Conn.x0    = synapse_dict[f"NT_per_vesicle_{conn_key[0]}"]
        Conn.tau_ar = synapse_dict[f"tau_asynchr_{conn_key[0]}"]
        Conn.U_max  = synapse_dict[f"sat_asynchr_{conn_key[0]}"]
        Conn.U_ar   = synapse_dict[f"strength_asynchr_{conn_key[0]}"]

    # Distance-dependent delays
    if network_dict.get("distance_delays", False):
        Vmax = network_dict["distance_delays_args"]["V_max"]
        Conn.delay = '(sqrt((x_pre - x_post)**2 + (y_pre - y_post)**2))/Vmax'
    print("add w_distribution ")
    return Conn