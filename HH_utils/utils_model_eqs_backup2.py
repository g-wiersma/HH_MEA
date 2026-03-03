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
from utils_helper_func import sort_tuple


#%% Single Neuron Equations: 

### HH equations ### 

eqs_HH_basis_E = '''
dV/dt = noise + (gl_E*(El_E-V)+g_na_E*(m*m*m)*h*(ENa_E-V)+g_kd_E*(n*n*n*n)*(EK_E-V)+I_tot)/Cm_E : volt
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dhp/dt = 0.128*exp((17.*mV-V+VT_E)/(18.*mV))/ms*(1.-hp)-4./(1+exp((30.*mV-V+VT_E)/(5.*mV)))/ms*h : 1
alpha_m = 0.32*(mV**-1)*4*mV/exprel((13*mV-V+VT_E)/(4*mV))/ms : Hz
beta_m = 0.28*(mV**-1)*5*mV/exprel(((V)-VT_E-40*mV)/(5*mV))/ms : Hz
alpha_h = 0.128*exp((17*mV-(V)+VT_E)/(18*mV))/ms : Hz
beta_h = 4./(1+exp((40*mV-(V)+VT_E)/(5*mV)))/ms : Hz
alpha_n = 0.032*(mV**-1)*5*mV/exprel((15*mV-V+VT_E)/(5*mV))/ms : Hz
beta_n = .5*exp((10*mV-V+VT_E)/(40*mV))/ms : Hz
I_AHP = g_AHP_E*Ca*(EK_E-V) : amp
dCa/dt = - Ca / tau_Ca_E : 1 
x : meter
y : meter
I_tot : amp
'''
eqs_HH_basis_I = '''
dV/dt = noise + (gl_I*(El_I-V)+g_na_I*(m*m*m)*h*(ENa_I-V)+g_kd_I*(n*n*n*n)*(EK_I-V)+I_tot_I)/Cm_I : volt
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dhp/dt = 0.128*exp((17.*mV-V+VT_I)/(18.*mV))/ms*(1.-hp)-4./(1+exp((30.*mV-V+VT_I)/(5.*mV)))/ms*h : 1
alpha_m = 0.32*(mV**-1)*4*mV/exprel((13*mV-(V)+VT_I)/(4*mV))/ms : Hz
beta_m = 0.28*(mV**-1)*5*mV/exprel(((V)-VT_I-40*mV)/(5*mV))/ms : Hz
alpha_h = 0.128*exp((17*mV-(V)+VT_I)/(18*mV))/ms : Hz
beta_h = 4./(1+exp((40*mV-(V)+VT_I)/(5*mV)))/ms : Hz
alpha_n = 0.032*(mV**-1)*5*mV/exprel((15*mV-V+VT_I)/(5*mV))/ms : Hz
beta_n = .5*exp((10*mV-V+VT_I)/(40*mV))/ms : Hz
I_AHP = g_AHP_I*Ca*(EK_I-V) : amp
dCa/dt = - Ca / tau_Ca_I : 1 
x : meter
y : meter
I : amp
'''

### AHP equations ###
eqs_I_AHP_E = '''
I_tot = I + I_AHP + I_syn : amp
I_AHP = g_AHP_E*Ca*(EK_E-V) : amp
dCa/dt = - Ca / tau_Ca_E : 1 
'''
eqs_I_AHP_I = '''
I_tot = I + I_AHP + I_syn : amp
I_AHP = g_AHP_I*Ca*(EK_I-V) : amp
dCa/dt = - Ca / tau_Ca_I : 1 
'''
eqs_I_no_AHP = '''
I_tot = I + I_syn : amp
'''



### Noise Equation: ###
eqs_neuronal_noise_E = '''
noise = sigma_E*(2*gl_E/Cm_E)**.5*randn()/sqrt(dt) : volt/second  # (more or less constant over dt)
'''

eqs_neuronal_noise_I = '''
noise = sigma_I*(2*gl_I/Cm_I)**.5*randn()/sqrt(dt) : volt/second  # (more or less constant over dt)
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
'''

### AMPA + asynchr_E basis ###
eqs_I_ampa_asynchr = '''
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_E_tot : 1 
qar_E_tot : Hz                        # (the asynchronous release rate (all (pre)synapses toghether))
'''

### NMDA basis ###
eqs_I_nmda = '''
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
s_nmda_tot :1
'''

### GABA basis ###
eqs_I_gaba = '''
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_gaba/dt = -s_gaba/tau_gaba :1
'''


### GABA + asynchr_I basis ###
eqs_I_gaba_asynchr = '''
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_gaba/dt = -s_gaba/tau_gaba + qar_I_tot:1
qar_I_tot : Hz
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
'''


###############################################################################

### Combine equations when combined mechanisms alter the ODEs ###

# Without NMDA

eqs_synapse_std_asynchr_E = '''
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
tau_d : second
tau_ar : second
U : 1
'''
eqs_synapse_std_asynchr_I = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
tau_d : second
tau_ar : second
U : 1
'''

eqs_synapse_std_stf_asynchr_E = '''
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
x0 : 1
tau_d : second
tau_f : second
tau_ar : second
U : 1
'''
eqs_synapse_std_stf_asynchr_I = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
x0 : 1
tau_d : second
tau_f : second
tau_ar : second
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
U : 1
'''

eqs_synapse_nmda_std_asynchr = '''
s_nmda_tot_post = w * S * x_d * s_nmda : 1 (summed) 
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
x0 : 1
tau_d : second
tau_ar : second
U : 1
'''

eqs_synapse_nmda_std_stf_asynchr = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda : 1 (summed)
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  # (scales correctly for different dt)
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
'''



###############################################################################

#%% ### OnPre Equations: ###

### Basic Equations: ###

eqs_onpre_std = '''
x_d *= (1-u_d)
'''
eqs_onpre_stf = '''
u_d += U*(1-u_d)
'''
eqs_onpre_asynchr = '''
uar += U_ar*(U_max-uar)
U_max : Hz
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

        eqs_neuron_E = eqs_HH_basis_E
        eqs_neuron_I = eqs_HH_basis_I
        
        ### Select AHP ###
        if {"AHP"} <= E_mechanisms:
            eqs_neuron_E += eqs_I_AHP_E
        else:
            eqs_neuron_E += eqs_I_no_AHP
    
        if {"AHP"} <= I_mechanisms:
            eqs_neuron_I += eqs_I_AHP_I
        else:
            eqs_neuron_I += eqs_I_no_AHP
        
        ### Select Synapse Types for I_syn = ... ###
        if ( ({"ampa"} <= E_mechanisms)
            and not ({"nmda"} & E_mechanisms)
            and not ({"gaba"} & I_mechanisms) ):
            eqs_neuron_E += eqs_I_syn_ampa_tot
            eqs_neuron_I += eqs_I_syn_ampa_tot
                 
        elif ( ({"ampa", "nmda"} <= E_mechanisms)            
            and not ({"gaba"} & I_mechanisms) ):
            eqs_neuron_E += eqs_I_syn_ampa_nmda_tot
            eqs_neuron_I += eqs_I_syn_ampa_nmda_tot
            
        elif ( ({"ampa"} <= E_mechanisms)
            and not ({"nmda"} & E_mechanisms)
            and ({"gaba"} <= I_mechanisms) ):
            eqs_neuron_E += eqs_I_syn_ampa_gaba_tot
            eqs_neuron_I += eqs_I_syn_ampa_gaba_tot
            
        elif ( {"ampa", "nmda"} <= E_mechanisms            
            and  ({"gaba"} <= I_mechanisms) ):
            eqs_neuron_E += eqs_I_syn_ampa_nmda_gaba_tot
            eqs_neuron_I += eqs_I_syn_ampa_nmda_gaba_tot
        
        else:
            raise ValueError("\nNo valid synaptic mechanisms/neurotransmitters selected as input!")
        
        ### Add ODEs for the synapses ###
        if {"nmda"} <= E_mechanisms:
            eqs_neuron_E += eqs_I_nmda
            eqs_neuron_I += eqs_I_nmda
        if "ampa" in E_mechanisms and "asynchr" in E_mechanisms:
            eqs_neuron_E += eqs_I_ampa_asynchr
            eqs_neuron_I += eqs_I_ampa_asynchr            
        if "ampa" in E_mechanisms and "asynchr" not in E_mechanisms: 
            eqs_neuron_E += eqs_I_ampa
            eqs_neuron_I += eqs_I_ampa
            
        if "gaba" in I_mechanisms and "asynchr" in I_mechanisms:
            eqs_neuron_E += eqs_I_gaba_asynchr
            eqs_neuron_I += eqs_I_gaba_asynchr
        if "gaba" in I_mechanisms and "asynchr" not in I_mechanisms:
            eqs_neuron_E += eqs_I_gaba
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
        if "std" in E_mechanisms: 
            eqs_onpre_E += eqs_onpre_std
        if "stf" in E_mechanisms: 
            eqs_onpre_E += eqs_onpre_stf
        
        # AMPA equation: 
        
        if not (E_mechanisms & {"std", "stf"}):
            eqs_onpre_E += eqs_onpre_ampa_basic
        
        if ( (E_mechanisms & {"std"})
            and not (E_mechanisms & {"stf"}) ):
            eqs_onpre_E += eqs_onpre_ampa_std
       
        if E_mechanisms & {"std", "stf"}:
            eqs_onpre_E += eqs_onpre_ampa_std_stf
       
        
        ### Inhibitory Onpre: ###
        
        if not (I_mechanisms & {"std", "stf"}):
            eqs_onpre_I += eqs_onpre_gaba_basic
       
        if ( (I_mechanisms & {"std"})
            and not (I_mechanisms & {"stf"}) ):
            eqs_onpre_I += eqs_onpre_gaba_std
      
        if I_mechanisms & {"std", "stf"}:
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
        mod_equations["E_synapse"] = E_synapse[E_mechanisms]
        mod_equations["E_onpre"] = E_onpre[E_mechanisms]
        
    if incl_I_neurons: 
        mod_equations["I_neuron"] = I_neuron
        mod_equations["I_synapse"] = I_synapse[I_mechanisms]
        mod_equations["I_onpre"] = I_onpre[I_mechanisms]
        

    return mod_equations