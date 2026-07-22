#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 16:59:26 2026

@author: Guido
"""



import random
import numpy as np
import pandas as pd
import pickle

### Prior and sampling functions: ###
    
PRIORS = {
    # Boolean parameters (equal probability)
    "include_I_neurons": [True, False],
    "include_AHP": [True, False],

    # Discrete parameter (equal probability)
    "p_I_neuron": [0.2, 0.4, 0.6, 0.8, 1],

    # Continuous parameters (uniform)
    "S_connect_e2e": (0.05, 10.0),
    "noise_sigma_E": (1.5, 8.0),    # mV

    "S_connect_e2i": (0.1, 10.0),
    "S_connect_i2e": (0.1, 10.0),
    "S_connect_i2i": (0.1, 10.0),
    "noise_sigma_I": (1.5, 8.0),    # mV

    "g_AHP": (0.2, 10.0), #  msiemens * cm ** -2) * area_E (e.g. 1000 umeter ** 2  )
}

def sample_parameter(priors=PRIORS):
    params = {}

    # Boolean parameters
    params["include_I_neurons"] = random.choice(priors["include_I_neurons"])
    params["include_AHP"] = random.choice(priors["include_AHP"])

    # Conditional parameters
    if params["include_I_neurons"]:
        params["p_I_neuron"] = random.choice(priors["p_I_neuron"])

        params["S_connect_e2i"] = np.random.uniform(*priors["S_connect_e2i"])
        params["S_connect_i2e"] = np.random.uniform(*priors["S_connect_i2e"])
        params["S_connect_i2i"] = np.random.uniform(*priors["S_connect_i2i"])
        params["noise_sigma_I"] = (
            np.random.uniform(*priors["noise_sigma_I"]) 
        )

        if params["p_I_neuron"] != 1.0:
            params["S_connect_e2e"] = np.random.uniform(*priors["S_connect_e2e"])
            params["noise_sigma_E"] = (
                np.random.uniform(*priors["noise_sigma_E"]) 
            )
    else:
        # No inhibitory neurons
        params["p_I_neuron"] = None
        params["S_connect_e2i"] = None
        params["S_connect_i2e"] = None
        params["S_connect_i2i"] = None
        params["noise_sigma_I"] = None

        params["S_connect_e2e"] = np.random.uniform(*priors["S_connect_e2e"])
        params["noise_sigma_E"] = (
            np.random.uniform(*priors["noise_sigma_E"]) 
        )

    if params["include_AHP"]:
        params["g_AHP"] = np.random.uniform(*priors["g_AHP"])
    else:
        params["g_AHP"] = None

    return params



def generate_prior_samples(n_samples):
    samples = []

    for _ in range(n_samples):
        params = sample_parameter()
        samples.append(params)

    return pd.DataFrame(samples)



prior_samples = generate_prior_samples(800000)


### Store prior samples: ###
prior_samples_name = "EI_pilot_full_samples.pkl"
storage_dir = "/home/Guido/UT_VU/HH_MEA_EI_pilot/" + prior_samples_name

with open(storage_dir, 'wb') as f: 
    pickle.dump(prior_samples, f)



