#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 15:18:08 2025

Utility file containing different neurological ODE set-ups to provide to the 
MEADataAnalysis.py file. 

@author: Guido
"""

import sys
sys.path.append("/home/Guido/UT_VU/py_code/HH_model/HH_utils")



#%% Single Neuron Equations: 

### HH equations ### 

eqs_HH_basis = '''
dV/dt = noise + (gl*(El-V) + g_na*m**3*h*(ENa-V) + g_kd*n**4*(EK-V) + I_tot)/Cm : volt
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dhp/dt = 0.128*exp((17.*mV-V+VT)/(18.*mV))/ms*(1.-hp) - 4./(1+exp((30.*mV-V+VT)/(5.*mV)))/ms*h : 1
alpha_m = 0.32*(mV**-1)*4*mV/exprel((13*mV-V+VT)/(4*mV))/ms : Hz
beta_m = 0.28*(mV**-1)*5*mV/exprel((V-VT-40*mV)/(5*mV))/ms : Hz
alpha_h = 0.128*exp((17*mV-(V)+VT)/(18*mV))/ms : Hz
beta_h = 4./(1+exp((40*mV-(V)+VT)/(5*mV)))/ms : Hz
alpha_n = 0.032*(mV**-1)*5*mV/exprel((15*mV-V+VT)/(5*mV))/ms : Hz
beta_n = 0.5*exp((10*mV-V+VT)/(40*mV))/ms : Hz
I : amp
x : meter
y : meter
Cm : farad
gl : siemens
El : volt
EK : volt
ENa : volt
g_na : siemens
g_kd : siemens
VT : volt
'''


### AHP equations ###

eqs_I_AHP = '''
I_tot = I + I_AHP + I_syn : amp
I_AHP = g_AHP*Ca*(EK-V) : amp
dCa/dt = - Ca / tau_Ca : 1
tau_Ca : second
g_AHP: siemens
alpha_Ca : 1
'''

# If AHP is not included
eqs_I_no_AHP = '''
I_tot = I + I_syn : amp
'''




### Noise Equation: ###
eqs_neuronal_noise = '''
noise = noise_sigma*(2*gl/Cm)**0.5*randn()/sqrt(dt) : volt/second (constant over dt)
noise_sigma : volt
'''
eqs_neuronal_no_noise = '''
noise : volt/second 
'''



###############################################################################

### I_syn = ... equations ###
eqs_I_syn_ampa_tot = '''
I_syn = I_ampa : amp
'''
eqs_I_syn_ampa_nmda_tot = '''
I_syn = I_ampa + I_nmda : amp
'''
eqs_I_syn_ampa_gaba_tot = '''
I_syn = I_ampa + I_gaba : amp
'''
eqs_I_syn_ampa_nmda_gaba_tot = '''
I_syn = I_ampa + I_nmda + I_gaba : amp
'''


###############################################################################

### I_syn differential equations: ###
        
### AMPA basis ###
eqs_I_ampa = '''
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa : 1 
g_ampa : siemens
E_ampa : volt
tau_ampa : second
'''

### AMPA + asynchr_E basis ###
eqs_I_ampa_asynchr = '''
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_E_tot : 1 
qar_E_tot : Hz                        # (the asynchronous release rate (all (pre)synapses toghether))
g_ampa : siemens
E_ampa : volt
tau_ampa : second
'''

### NMDA basis ###
eqs_I_nmda = '''
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
s_nmda_tot :1
E_nmda : volt
g_nmda : siemens
'''

### GABA basis ###
eqs_I_gaba = '''
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_gaba/dt = -s_gaba/tau_gaba :1
g_gaba : siemens
E_gaba : volt
tau_gaba : second
'''


### GABA + asynchr_I basis ###
eqs_I_gaba_asynchr = '''
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_gaba/dt = -s_gaba/tau_gaba + qar_I_tot:1
qar_I_tot : Hz
g_gaba : siemens
E_gaba : volt
tau_gaba : second
'''





#%% Synapse Equations:

### Basic Equations ###    

eqs_synapse_basis = '''
w : 1
S : 1
'''
eqs_synapse_std = '''
dx_d/dt = (1-x_d)/tau_d :1 (clock-driven)
tau_d : second
U : 1
'''
eqs_synapse_stf = '''
du_d/dt = (U-u_d)/tau_f : 1 (clock-driven)
tau_f : second
'''
eqs_synapse_nmda = '''
s_nmda_tot_post = w * S * s_nmda : 1 (summed)           # connection weight * synaptic strength * part of op receptor at some time t, summed over all synapses
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
taus_nmda : second
alpha_nmda : Hz
taux_nmda : second
'''


###############################################################################

### Combine equations when combined mechanisms alter the ODEs ###

# Without NMDA

eqs_synapse_std_asynchr_E = '''
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt : Hz  (constant over dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
tau_d : second
tau_ar : second
U_max : Hz
U_ar : 1
U : 1
'''
eqs_synapse_std_asynchr_I = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt : Hz (constant over dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
tau_d : second
tau_ar : second
U_max : Hz
U_ar : 1
U : 1
'''

eqs_synapse_std_stf_asynchr_E = '''
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt : Hz  (constant over dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
x0 : 1
tau_d : second
tau_f : second
tau_ar : second
U_max : Hz
U_ar : 1
U : 1
'''
eqs_synapse_std_stf_asynchr_I = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt : Hz  (constant over dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
x0 : 1
tau_d : second
tau_f : second
tau_ar : second
U_max : Hz
U_ar : 1
U : 1
'''


###############################################################################

# With NMDA

eqs_synapse_nmda_std = '''
s_nmda_tot_post = w * S * s_nmda * x_d : 1 (summed)  # compared to nmda without std, here x_d is added 
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d : 1 (clock-driven)
tau_d : second
taus_nmda : second
alpha_nmda : Hz
taux_nmda : second
U : 1
'''

eqs_synapse_nmda_std_stf = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda  :1 (summed)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*x_nmda*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d :1 (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
tau_d : second
tau_f : second
taus_nmda : second
alpha_nmda : Hz
taux_nmda : second
U : 1
'''

eqs_synapse_nmda_std_asynchr = '''
s_nmda_tot_post = w * S * x_d * s_nmda : 1 (summed) 
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt : Hz  (constant over dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
taus_nmda : second
alpha_nmda : Hz
taux_nmda : second
tau_d : second
tau_ar : second
U_max : Hz
U_ar : 1
U : 1
'''

eqs_synapse_nmda_std_stf_asynchr = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda : 1 (summed)
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt : Hz  (constant over dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
U : 1
tau_ar : second
x0 : 1
tau_d : second
tau_f : second
taus_nmda : second
alpha_nmda : Hz
taux_nmda : second
U_max : Hz
U_ar : 1
'''



###############################################################################

#%% ### OnPre Equations: ###

### Basic Equations: ###

eqs_onpre_std = '''
x_d *= (1-U)
'''

eqs_onpre_std_stf = '''
x_d *= (1-u_d)
u_d += U*(1-u_d)
'''
eqs_onpre_asynchr = '''
uar += U_ar*(U_max-uar)
'''
eqs_onpre_nmda = '''
x_nmda += 1
'''

### AMPA mechanism-dependent equations: ###

eqs_onpre_ampa_basic = '''
s_ampa += w * S 
'''
eqs_onpre_ampa_std = '''
s_ampa += w * x_d * S
'''
eqs_onpre_ampa_std_stf = '''
s_ampa += w * x_d * u_d * S
'''

### GABA mechanism-dependent equations: ###

eqs_onpre_gaba_basic = '''
s_gaba += w * S 
'''
eqs_onpre_gaba_std = '''
s_gaba += w * x_d * S
'''
eqs_onpre_gaba_std_stf = '''
s_gaba += w * x_d * u_d * S
'''



#%% Import function 


def get_equations_HH(dict_model_config):
    """
    Function to return the model equations corresponding to the requested 
    model construction. All equations follow the Brian2 format.
    
    Input: 
        dict_model_config (dict)
            dictionary with the following structure: {
            "excitatory_neurons": True/False,
            "inhibitory_neurons": True/False,        
            "E_mechanisms": ["mechanism_1", "mechanism_2"],
            "I_mechanisms": ["mechanism_1", "mechanism_2"]
            }
            
        example of E only neurons with ampa and std: 
        dict_model_config = {"excitatory_neurons": True, 
                             "E_mechanisms": ["ampa", "std]}                                              
                            }  
  Output: 
      mod_equations (dict)
          dictionary with the model equations. It contains the neuron equations, 
          the synaptic equations, and the corresponding pre-synaptic equations. 
          It has the following structure: 
              {
              mod_equations["E_neuron"] = E_neuron equation
              mod_equations["E_synapse"] = E_synapse equation
              mod_equations["E_onpre"] = E_onpre equation
              (and/or)
              mod_equations["I_neuron"] = I_neuron equation
              mod_equations["I_synapse"] = I_neuron equation
              mod_equations["I_onpre"] = I_onpre equation 
              }                                             
    """
    
    ###########################################################################
    def get_neuron_eq(E_mechanisms, I_mechanisms):      
        """
        Construct brian2 neuron model equations for excitatory (E) and inhibitory (I)
        populations based on selected mechanisms.
    
        Input: 
            E_mechanisms : set of strings
                Set of enabled mechanisms for excitatory neurons.
                Possible entries include:
                    "AHP"       : afterhyperpolarization current
                    "ampa"      : AMPA synaptic current
                    "nmda"      : NMDA synaptic current
                    "asynchr"   : asynchronous release component
    
            I_mechanisms : set of strings
                Set of enabled mechanisms for inhibitory neurons.
                Possible entries include:
                    "AHP"       : afterhyperpolarization current
                    "gaba"      : GABA synaptic current
                    "asynchr"   : asynchronous release component
    
        Output:
            eqs_neuron_E : str
                Complete equation string for excitatory neurons, assembled
                from Hodgkin–Huxley base equations and selected mechanisms.
        
            eqs_neuron_I : str
                Complete equation string for inhibitory neurons, assembled
                from Hodgkin–Huxley base equations and selected mechanisms.
        """

        eqs_neuron_E = eqs_HH_basis
        eqs_neuron_I = eqs_HH_basis
        
        ### Select AHP ###
        if {"AHP"} <= E_mechanisms:
            eqs_neuron_E += eqs_I_AHP
        else:
            eqs_neuron_E += eqs_I_no_AHP
        
        if {"AHP"} <= I_mechanisms:
            eqs_neuron_I += eqs_I_AHP
        else:
            eqs_neuron_I += eqs_I_no_AHP
        
        if {"neuronal_noise"} <= E_mechanisms:
            eqs_neuron_E += eqs_neuronal_noise
        else: 
            eqs_neuron_E += eqs_neuronal_no_noise
            
        if {"neuronal_noise"} <= I_mechanisms:
            eqs_neuron_I += eqs_neuronal_noise
        else: 
            eqs_neuron_I += eqs_neuronal_no_noise    
            
            
            
        ### Select Synapse Types for I_syn (post-synaptic specific) ###

        # Excitatory neurons 
        if ({"ampa"} <= E_mechanisms
            and not ({"nmda"} & E_mechanisms)
            and not ({"gaba"} & E_mechanisms)):
            eqs_neuron_E += eqs_I_syn_ampa_tot
                
        elif ({"ampa", "nmda"} <= E_mechanisms
            and not ({"gaba"} & E_mechanisms)):            
            eqs_neuron_E += eqs_I_syn_ampa_nmda_tot       
        
        elif ({"ampa", "gaba"} <= E_mechanisms
            and not ({"nmda"} & E_mechanisms)):            
            eqs_neuron_E += eqs_I_syn_ampa_gaba_tot
                
        elif ({"ampa", "nmda"} <= E_mechanisms
            and {"gaba"} <= E_mechanisms):            
            eqs_neuron_E += eqs_I_syn_ampa_nmda_gaba_tot
        
        
        
        # Inhibitory neurons
        if ({"ampa"} <= I_mechanisms
            and not ({"nmda"} & I_mechanisms)
            and not ({"gaba"} & I_mechanisms)):            
            eqs_neuron_I += eqs_I_syn_ampa_tot
                
        elif ({"ampa", "nmda"} <= I_mechanisms
            and not ({"gaba"} & I_mechanisms)):            
            eqs_neuron_I += eqs_I_syn_ampa_nmda_tot
               
        elif ({"ampa", "gaba"} <= I_mechanisms
            and not ({"nmda"} & I_mechanisms)):        
            eqs_neuron_I += eqs_I_syn_ampa_gaba_tot
                
        elif ({"ampa", "nmda", "gaba"} <= I_mechanisms):            
            eqs_neuron_I += eqs_I_syn_ampa_nmda_gaba_tot
        
        else:
            raise ValueError("\nNo valid synaptic mechanisms/neurotransmitters selected as input!")
        
        
        ### Add ODEs for the synapses ###
        if {"nmda"} <= E_mechanisms:
            eqs_neuron_E += eqs_I_nmda
        if {"nmda"} <= I_mechanisms:
            eqs_neuron_I += eqs_I_nmda        
            
        if "ampa" in E_mechanisms and "asynchr" in E_mechanisms:
            eqs_neuron_E += eqs_I_ampa_asynchr
        if "ampa" in I_mechanisms and "asynchr" in E_mechanisms:
            eqs_neuron_I += eqs_I_ampa_asynchr            
        if "ampa" in E_mechanisms and "asynchr" not in E_mechanisms: 
            eqs_neuron_E += eqs_I_ampa
        if "ampa" in I_mechanisms and "asynchr" not in E_mechanisms: 
            eqs_neuron_I += eqs_I_ampa
            
        if "gaba" in E_mechanisms and "asynchr" in I_mechanisms:
            eqs_neuron_E += eqs_I_gaba_asynchr
        if "gaba" in I_mechanisms and "asynchr" in I_mechanisms:
            eqs_neuron_I += eqs_I_gaba_asynchr

        if "gaba" in E_mechanisms and "asynchr" not in I_mechanisms:
            eqs_neuron_E += eqs_I_gaba
        if "gaba" in I_mechanisms and "asynchr" not in I_mechanisms:
            eqs_neuron_I += eqs_I_gaba
            
        return eqs_neuron_E, eqs_neuron_I
        
    
    def get_synapse_eq(E_mechanisms, I_mechanisms):      
        """
        Construct brian2 synapse model equations for excitatory (E) and inhibitory (I)
        populations based on selected mechanisms.

        Input: 
            E_mechanisms : set of strings
                Set of enabled mechanisms for excitatory neurons.
                Possible entries include:
                    "std"       : short term depression component
                    "nmda"      : NMDA synaptic current
                    "asynchr"   : asynchronous release component

            I_mechanisms : set of strings
                Set of enabled mechanisms for inhibitory neurons.
                Possible entries include:
                    "std"       : short term depression component
                    "stf"       : short term facilitation component
                    "asynchr"   : asynchronous release component

        Output:
            eqs_synapse_E : str
                Complete equation string for excitatory synapses, assembled
                from selected mechanisms.
        
            eqs_synapse_I : str
                Complete equation string for inhibitory synapses, assembled
                from selected mechanisms.
        """

        eqs_synapse_E = eqs_synapse_basis
        eqs_synapse_I = eqs_synapse_basis
        
        ###########################################################################
        
        ### Check compatibility: ###
        
        if ( (E_mechanisms & {"stf", "asynchr"} and "std" not in E_mechanisms)
            or ("nmda" in E_mechanisms and "ampa" not in E_mechanisms)
            or (I_mechanisms & {"stf", "asynchr", "asynchr"} and "gaba" not in I_mechanisms) ):      
            raise ValueError(
                "\nNo valid/implemented synaptic mechanisms/neurotransmitters selected as input!"
            )
        
        
        ### Excitatory Synapse: ###
        
        # No NMDA/Asynchr. release
        if ( ({"std"} <= E_mechanisms)
            and not ({"nmda"} & E_mechanisms)
            and not ({"asynchr"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_std
            
            if {"stf"} <= E_mechanisms: 
                eqs_synapse_E += eqs_synapse_stf
        
        # No NMDA, incl. asynchr. release
        if ( ({"asynchr", "std"} <= E_mechanisms)
            and not ({"nmda"} & E_mechanisms)
            and not ({"stf"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_std_asynchr_E 
        
        if ( ({"asynchr", "std", "stf"} <= E_mechanisms)
            and not ({"nmda"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_std_stf_asynchr_E 
        
        # Incl. NMDA
        if ( ({"nmda"} <= E_mechanisms)
            and not ({"std"} & E_mechanisms)
            and not ({"asynchr"} & E_mechanisms) 
            and not ({"stf"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_nmda
            
        if ( ({"nmda", "std"} <= E_mechanisms)
            and not ({"asynchr"} & E_mechanisms) 
            and not ({"stf"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_nmda_std
        
        if ( ({"nmda", "std", "stf"} <= E_mechanisms)
            and not ({"asynchr"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_nmda_std_stf
        
        if ( ({"nmda", "std", "asynchr"} <= E_mechanisms)
            and not ({"stf"} & E_mechanisms) ):
            eqs_synapse_E += eqs_synapse_nmda_std_asynchr
        
        if ( {"nmda", "std", "stf", "asynchr"} <= E_mechanisms):
            eqs_synapse_E += eqs_synapse_nmda_std_stf_asynchr
        
       
        
        ###########################################################################
        
        ### Inhibitory Synapse: ###
        
        # No Asynchr. release
        if ( ({"std"} <= I_mechanisms)
            and not ({"asynchr"} & E_mechanisms) ):
            eqs_synapse_I += eqs_synapse_std
            
            if {"stf"} <= I_mechanisms: 
                eqs_synapse_I += eqs_synapse_stf
        
        if ( ({"std", "asynchr"} <= I_mechanisms)
            and not ({"stf"} & E_mechanisms) ):
            eqs_synapse_I += eqs_synapse_std_asynchr_I
        
        if ({"std", "asynchr", "stf"} <= I_mechanisms): 
            eqs_synapse_I += eqs_synapse_std_stf_asynchr_I
            
            
        return eqs_synapse_E, eqs_synapse_I
    

    def get_onpre_eq(E_mechanisms, I_mechanisms):      
        """
        Construct the stochastic brian2 pre-synaptic model equations for 
        excitatory (E) and inhibitory (I) populations based on selected mechanisms.
    
        Input: 
            E_mechanisms : set of strings
                Set of enabled mechanisms for excitatory neurons.
                Possible entries include:
                    "std"       : short term depression component
                    "nmda"      : NMDA synaptic current
                    "asynchr"   : asynchronous release component
    
            I_mechanisms : set of strings
                Set of enabled mechanisms for inhibitory neurons.
                Possible entries include:
                    "std"       : short term depression component
                    "stf"       : short term facilitation component
                    "asynchr"   : asynchronous release component
    
        Output:
            eqs_onpre_E : str
                Complete equation string for excitatory pre-synaptic equations,
                assembled from selected mechanisms.
        
            eqs_synapse_I : str
                Complete equation string for inhibitory pre-synaptic equations,
                assembled from selected mechanisms.
        """
    
        eqs_onpre_E = ''''''
        eqs_onpre_I = ''''''
        
        ###########################################################################
        
        ### Excitatory Onpre: ###
        if "nmda" in E_mechanisms: 
            eqs_onpre_E += eqs_onpre_nmda
 
        if "asynchr" in E_mechanisms:
            eqs_onpre_E += eqs_onpre_asynchr
        
        if not (E_mechanisms & {"std", "stf"}):
            eqs_onpre_E += eqs_onpre_ampa_basic
        
        if ( (E_mechanisms & {"std"})
            and not (E_mechanisms & {"stf"}) ):
            eqs_onpre_E += eqs_onpre_std
            eqs_onpre_E += eqs_onpre_ampa_std            
            
        if E_mechanisms <= {"std", "stf"}:
            eqs_onpre_E += eqs_onpre_std_stf
            eqs_onpre_E += eqs_onpre_ampa_std_stf
       
        
       ### Inhibitory Onpre: ###
        
        if "asynchr" in I_mechanisms:
            eqs_onpre_I += eqs_onpre_asynchr
            
        if not (I_mechanisms & {"std", "stf"}):
            eqs_onpre_I += eqs_onpre_gaba_basic
       
        if ( (I_mechanisms & {"std"})
            and not (I_mechanisms & {"stf"}) ):
            eqs_onpre_I += eqs_onpre_std
            eqs_onpre_I += eqs_onpre_gaba_std
      
        if I_mechanisms <= {"std", "stf"}:
            eqs_onpre_I += eqs_onpre_std_stf
            eqs_onpre_I += eqs_onpre_gaba_std_stf


        return eqs_onpre_E, eqs_onpre_I


    ###########################################################################
    
    ### Model components: ###
    incl_E_neurons = dict_model_config["excitatory_neurons"]   # this is a boolean    
    incl_I_neurons = dict_model_config["inhibitory_neurons"]   
    if incl_E_neurons: 
        E_mechanisms = dict_model_config["E_mechanisms"] # contains a set of mechanisms
    else: 
        E_mechanisms = set()
    if incl_I_neurons: 
        I_mechanisms = dict_model_config["I_mechanisms"] # contains a set of mechanisms
    else: 
        I_mechanisms = set()
    
    ### Construct equation-output-dict: ###
    mod_equations = {}

    E_neuron, I_neuron = get_neuron_eq(E_mechanisms, I_mechanisms)
    E_synapse, I_synapse = get_synapse_eq(E_mechanisms, I_mechanisms)
    E_onpre, I_onpre = get_onpre_eq(E_mechanisms, I_mechanisms)
    
    
    if incl_E_neurons:         
        mod_equations["E_neuron"] = E_neuron        
        mod_equations["E_synapse"] = E_synapse
        mod_equations["E_onpre"] = E_onpre
        
    if incl_I_neurons: 
        mod_equations["I_neuron"] = I_neuron
        mod_equations["I_synapse"] = I_synapse
        mod_equations["I_onpre"] = I_onpre
        

    return mod_equations