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


#%% Single Neuron equations

### HH equations: ###
eqs_HH_E = '''
dV/dt = noise + (gl_E*(El_E-V)+g_na_E*(m*m*m)*h*(ENa_E-V)+g_kd_E*(n*n*n*n)*(EK_E-V)+I+I_syn+I_AHP)/Cm_E : volt
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
I : amp
'''

eqs_HH_I = '''
dV/dt = noise + (gl_I*(El_I-V)+g_na_I*(m*m*m)*h*(ENa_I-V)+g_kd_I*(n*n*n*n)*(EK_I-V)+I+I_syn+I_AHP)/Cm_I : volt
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

### Noise Equation: ###
eqs_neuronal_noise_E = '''
noise = sigma_E*(2*gl_E/Cm_E)**.5*randn()/sqrt(dt) : volt/second (constant over dt)
'''

eqs_neuronal_noise_I = '''
noise = sigma_I*(2*gl_I/Cm_I)**.5*randn()/sqrt(dt) : volt/second (constant over dt)
'''

#%% Excitatory only cultures
###############################################################################

### only AMPA ### 

eqs_I_syn_ampa = '''
I_syn =  I_ampa : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa : 1 
'''
eqs_neuron_ampa = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa
eqs_synapse_ampa = '''
w : 1
S : 1
'''
eqs_onpre_ampa = '''
s_ampa += w * S
'''



###############################################################################

### AMPA + std ### 
eqs_neuron_ampa_stdE = eqs_neuron_ampa
eqs_synapse_ampa_stdE = '''
dx_d/dt = (1-x_d)/tau_d :1 (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_stdE = '''
x_d *= (1-U)
s_ampa += w * S * x_d 
'''



###############################################################################

### AMPA + std + stf ###
eqs_neuron_ampa_stdE_stfE = eqs_neuron_ampa
eqs_synapse_ampa_stdE_stfE = '''
dx_d/dt = (1-x_d)/tau_d : 1 (clock-driven)
du_d/dt = (U-u_d)/tau_f : 1 (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_stdE_stfE = '''
x_d *= (1-u_d)
u_d += U*(1-u_d)
s_ampa += w * x_d * u_d * S
'''


###############################################################################

### AMPA + std + asynchr. release ###
eqs_I_syn_ampa_stdE_asynchrE = '''
I_syn  = I_ampa : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_tot : 1 
qar_tot : Hz                        # the asynchronous release rate (all (pre)synapses toghether)
'''
eqs_neuron_ampa_stdE_asynchrE = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_stdE_asynchrE
eqs_synapse_ampa_stdE_asynchrE = '''
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_stdE_asynchrE = '''
x_d *= (1-U) 
uar += U_ar*(U_max-uar)
s_ampa += w * x_d * S
'''


###############################################################################

### AMPA + std + stf + asynchr. release ###
eqs_neuron_ampa_stdE_stfE_asynchrE = eqs_neuron_ampa_stdE_asynchrE
eqs_synapse_ampa_stdE_stfE_asynchrE = '''
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
w : 1
S : 1
U : 1
'''
eqs_onpre_ampa_stdE_stfE_asynchrE = '''
x_d *= (1-U) 
u_d += U*(1-u_d)
uar += U_ar*(U_max-uar)
s_ampa += w * x_d * u_d * S
'''



###############################################################################

### AMPA + NMDA ###
eqs_I_syn_ampa_nmda = '''
I_syn =  I_ampa + I_nmda: amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
ds_ampa/dt = -s_ampa/tau_ampa : 1 
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
s_nmda_tot :1
'''
eqs_neuron_ampa_nmda = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_nmda
eqs_synapse_ampa_nmda = '''
s_nmda_tot_post = w * S * s_nmda : 1 (summed)           # connection weight * synaptic strength * part of op receptor at some time t, summed over all synapses
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
w : 1
S : 1
'''
eqs_onpre_ampa_nmda = '''
s_ampa += w * S           # postsyn s_ampa (in neuronal equation)
x_nmda += 1               # increment per-synapse NMDA gating (in synaptic equation)
'''



###############################################################################

### AMPA + NMDA + short term depression ###
eqs_neuron_ampa_nmda_stdE = eqs_neuron_ampa_nmda
eqs_synapse_ampa_nmda_stdE = '''
s_nmda_tot_post = w * S * s_nmda * x_d :1 (summed) 
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d : 1 (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_nmda_stdE= '''
x_nmda += 1 
x_d *= (1-U) 
s_ampa += w * x_d * S
'''



###############################################################################

### AMPA + NMDA + short term depression + short term facilitation ###
eqs_neuron_ampa_nmda_stdE_stfE = eqs_neuron_ampa_nmda
eqs_synapse_ampa_nmda_stdE_stfE = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda  :1 (summed)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*x_nmda*(1-s_nmda) : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) :1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d :1 (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_nmda_stdE_stfE = '''
x_nmda += 1
x_d *= (1-u_d)
s_ampa += w * S * x_d * u_d
u_d += U*(1-u_d)
'''



###############################################################################

### AMPA + NMDA + short term depression + asynchr. release ###
eqs_I_syn_ampa_nmda_stdE_asynchrE = '''
I_syn  = I_ampa + I_nmda : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_tot : 1 
s_nmda_tot : 1
qar_tot : Hz                        # the asynchronous release rate (all (pre)synapses toghether)
'''
eqs_neuron_ampa_nmda_stdE_asynchrE = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_nmda_stdE_asynchrE
eqs_synapse_ampa_nmda_stdE_asynchrE = '''
s_nmda_tot_post = w * S * x_d * s_nmda : 1 (summed) 
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_ampa_nmda_stdE_asynchrE = '''
x_nmda += 1 
x_d *= (1-U) 
uar += U_ar*(U_max-uar)
s_ampa += w * x_d * S
'''


###############################################################################

### AMPA + NMDA + short term depression + short term facilitation + asynchr. release ###
eqs_neuron_ampa_nmda_stdE_stfE_asynchrE = eqs_neuron_ampa_nmda_stdE_asynchrE
eqs_synapse_ampa_nmda_stdE_stfE_asynchrE = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda : 1 (summed)
qar_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
w : 1
S : 1
U : 1
'''
eqs_onpre_ampa_nmda_stdE_stfE_asynchrE = '''
x_nmda += 1 
x_d *= (1-U) 
u_d += U*(1-u_d)
uar += U_ar*(U_max-uar)
s_ampa += w * x_d * u_d * S
'''








#%% Excitatory + inhibitory neuronal cultures

###############################################################################

### AMPA + GABA ### 
eqs_I_syn_ampa_gaba = '''
I_syn =  I_ampa + I_gaba : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_ampa/dt = -s_ampa/tau_ampa :1 
ds_gaba/dt = -s_gaba/tau_gaba :1
'''
eqs_neuron_E_ampa_gaba = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_gaba
eqs_neuron_I_ampa_gaba = eqs_HH_I + eqs_neuronal_noise_I + eqs_I_syn_ampa_gaba
eqs_synapse_E_ampa_gaba = '''
w : 1
S : 1
'''
eqs_synapse_I_ampa_gaba = '''
w : 1
S : 1
'''
eqs_onpre_E_ampa_gaba = eqs_onpre_ampa
eqs_onpre_I_ampa_gaba = '''
s_gaba += w * S 
'''


###############################################################################

### AMPA + GABA + std_E  ### 
eqs_neuron_E_ampa_gaba_stdE = eqs_neuron_E_ampa_gaba
eqs_neuron_I_ampa_gaba_stdE = eqs_neuron_I_ampa_gaba
eqs_synapse_E_ampa_gaba_stdE = eqs_synapse_ampa_stdE
eqs_synapse_I_ampa_gaba_stdE = eqs_synapse_I_ampa_gaba         
eqs_onpre_E_ampa_gaba_stdE = eqs_onpre_ampa_stdE
eqs_onpre_I_ampa_gaba_stdE = eqs_onpre_I_ampa_gaba

###############################################################################

### AMPA + GABA + std_E + stf_E ###

eqs_neuron_E_ampa_gaba_stdE = eqs_neuron_E_ampa_gaba
eqs_neuron_I_ampa_gaba_stdE = eqs_neuron_I_ampa_gaba
eqs_synapse_E_ampa_gaba_stdE = eqs_synapse_ampa_stdE_stfE
eqs_synapse_I_ampa_gaba_stdE = eqs_synapse_I_ampa_gaba         
eqs_onpre_E_ampa_gaba_stdE = eqs_onpre_ampa_stdE_stfE
eqs_onpre_I_ampa_gaba_stdE = eqs_onpre_I_ampa_gaba


###############################################################################

### AMPA + GABA + std_E + asynchr_E ###
### CHECK
eqs_I_syn_ampa_gaba_nmda_stdE_asynchrE = '''
I_syn  = I_ampa + I_gaba + I_nmda : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_gaba = g_gaba*(E_gaba-V)*s_gaba : amp
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_E_tot : 1 
ds_gaba/dt = -s_gaba/tau_gaba : 1
s_nmda_tot : 1
qar_E_tot : Hz                        # the asynchronous release rate (all exc. (pre)synapses toghether)
'''
eqs_neuron_E_ampa_gaba_nmda_stdE_asynchrE = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_gaba_nmda_stdE_asynchrE
eqs_neuron_E_ampa_gaba_stdE = eqs_neuron_E_ampa_gaba
eqs_neuron_I_ampa_gaba_stdE = eqs_neuron_I_ampa_gaba
eqs_synapse_E_ampa_gaba_stdE = eqs_synapse_ampa_stdE_asynchrE
eqs_synapse_I_ampa_gaba_stdE = eqs_synapse_I_ampa_gaba         
eqs_onpre_E_ampa_gaba_stdE = eqs_onpre_ampa_stdE_asynchrE
eqs_onpre_I_ampa_gaba_stdE = eqs_onpre_I_ampa_gaba



###############################################################################

### AMPA + GABA + std_E + std_I ### 
## #CHECK
eqs_neuron_E_ampa_gaba_stdE_stdI = eqs_neuron_E_ampa_gaba
eqs_neuron_I_ampa_gaba_stdE_stdI = eqs_neuron_I_ampa_gaba
eqs_synapse_E_ampa_gaba_stdE_stdI = eqs_synapse_ampa_stdE
eqs_synapse_I_ampa_gaba_stdE_stdI = eqs_synapse_ampa_stdE          # same ODE for another NT
eqs_onpre_E_ampa_gaba_stdE_stdI = eqs_onpre_ampa_stdE
eqs_onpre_I_ampa_gaba_stdE_stdI = '''
x_d *= (1-U)
s_gaba += w * S * x_d
'''


###############################################################################

### AMPA + GABA + stdE + stdI + asynchrE + asynchrI

eqs_I_syn_ampa_gaba_stdE_stdI_asynchrE_asynchrI = '''
I_syn  = I_ampa + I_gaba + I_nmda : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_E_tot : 1 
ds_gaba/dt = -s_gaba/tau_gaba + qar_I_tot:1
qar_E_tot : Hz                        # the asynchronous release rate (all exc. (pre)synapses toghether)
qar_I_tot : Hz
'''
eqs_neuron_E_ampa_gaba_stdE_stdI_asynchrE_asynchrI = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_gaba_stdE_stdI_asynchrE_asynchrI
eqs_neuron_I_ampa_gaba_stdE_stdI_asynchrE_asynchrI = eqs_HH_I + eqs_neuronal_noise_I + eqs_I_syn_ampa_gaba_stdE_stdI_asynchrE_asynchrI
eqs_synapse_E_ampa_gaba_stdE_stdI_asynchrE_asynchrI = eqs_synapse_ampa_stdE_asynchrE
eqs_synapse_I_ampa_gaba_stdE_stdI_asynchrE_asynchrI = eqs_synapse_ampa_stdE_asynchrE # Same ODE as for E-synapses


###############################################################################

### AMPA + GABA + stdE + stdI + stfE + stfI ###
## #CHECK
eqs_neuron_E_ampa_gaba_stdE_stdI_stfE_stfI = eqs_neuron_E_ampa_gaba
eqs_neuron_I_ampa_gaba_stdE_stdI_stfE_stfI = eqs_neuron_I_ampa_gaba
eqs_synapse_E_ampa_gaba_stdE_stdI_stfE_stfI = eqs_synapse_ampa_stdE_stfE
eqs_synapse_I_ampa_gaba_stdE_stdI_stfE_stfI = eqs_synapse_ampa_stdE_stfE # same ODE as for E-synapses
eqs_onpre_E_ampa_gaba_stdE_stdI_stfE_stfI = eqs_onpre_ampa_stdE_stfE
eqs_onpre_I_ampa_gaba_stdE_stdI_stfE_stfI = '''
x_d *= (1-u_d)
u_d += U*(1-u_d)
s_gaba += w * S * x_d * u_d 
'''




###############################################################################

# AMPA + GABA + NMDA
eqs_I_syn_ampa_gaba_nmda = '''
I_syn =  I_ampa + I_gaba + I_nmda : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
ds_ampa/dt = -s_ampa/tau_ampa :1 
ds_gaba/dt = -s_gaba/tau_gaba :1
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
s_nmda_tot :1
'''
eqs_neuron_E_ampa_gaba_nmda = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_gaba_nmda
eqs_neuron_I_ampa_gaba_nmda = eqs_HH_I + eqs_neuronal_noise_I + eqs_I_syn_ampa_gaba_nmda
eqs_synapse_E_ampa_gaba_nmda = eqs_synapse_ampa_nmda
eqs_synapse_I_ampa_gaba_nmda = eqs_synapse_ampa_nmda
eqs_onpre_E_ampa_gaba_nmda = eqs_onpre_ampa_nmda
eqs_onpre_I_ampa_gaba_nmda = eqs_onpre_I_ampa_gaba




###############################################################################

# AMPA + GABA + NMDA + stdE + stdI
#CHECK
eqs_neuron_E_ampa_gaba_nmda_stdE_stdI= eqs_neuron_E_ampa_gaba_nmda
eqs_neuron_I_ampa_gaba_nmda_stdE_stdI = eqs_neuron_I_ampa_gaba_nmda
eqs_synapse_E_ampa_gaba_nmda_stdE_stdI = eqs_synapse_ampa_nmda_stdE
eqs_synapse_I_ampa_gaba_nmda_stdE_stdI = eqs_synapse_I_ampa_gaba_stdE #HERE
eqs_onpre_E_ampa_gaba_nmda_stdE_stdI = eqs_onpre_ampa_nmda_stdE
eqs_onpre_I_ampa_gaba_nmda_stdE_stdI = eqs_onpre_I_ampa_gaba_stdE




###############################################################################

# AMPA + GABA + NMDA + stdE + stfE
#CHECK
eqs_neuron_E_ampa_gaba_nmda_stdE_stfE = eqs_neuron_E_ampa_gaba_nmda
eqs_neuron_I_ampa_gaba_nmda_stdE_stfE = eqs_neuron_I_ampa_gaba_nmda
eqs_synapse_E_ampa_gaba_nmda_stdE_stfE = eqs_synapse_ampa_nmda_stdE_stfE
eqs_synapse_I_ampa_gaba_nmda_stdE_stfE = eqs_synapse_I_ampa_gaba_stdE_stfE
eqs_onpre_E_ampa_gaba_nmda_stdE_stfE = eqs_onpre_ampa_nmda_stdE_stfE
eqs_onpre_I_ampa_gaba_nmda_stdE_stfE = eqs_onpre_I_ampa_gaba_stdE_stfE

    


###############################################################################

# AMPA + GABA + NMDA + stdE + stdI + asynchr. release_E + asynchr. release_I
## #CHECK
eqs_I_syn_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = '''
I_syn  = I_ampa + I_gaba + I_nmda : amp
I_ampa = g_ampa*(E_ampa-V)*s_ampa : amp
I_gaba = g_gaba*(E_gaba-V)*s_gaba :amp
I_nmda = g_nmda*(E_nmda-V)*s_nmda_tot/(1+exp(-0.062*V/mV)/3.57) : amp
ds_ampa/dt = -s_ampa/tau_ampa + qar_E_tot : 1 
ds_gaba/dt = -s_gaba/tau_gaba + qar_I_tot:1
s_nmda_tot : 1
qar_E_tot : Hz                        # the asynchronous release rate (all exc. (pre)synapses toghether)
qar_I_tot : Hz
'''
eqs_neuron_E_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = eqs_HH_E + eqs_neuronal_noise_E + eqs_I_syn_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI
eqs_neuron_I_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = eqs_HH_I + eqs_neuronal_noise_I + eqs_I_syn_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI
eqs_synapse_E_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = '''
s_nmda_tot_post = w * S * x_d * s_nmda : 1 (summed) 
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_synapse_I_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d/x0*uar*dt*(1-uar*dt))+uar*dt*x_d/x0, 0, 2*x_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_E_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = eqs_onpre_ampa_nmda_stdE_asynchrE
eqs_onpre_I_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI = '''
x_d *= (1-U)
s_gaba += w * S * x_d
uar += U_ar*(U_max-uar)
'''





###############################################################################

# AMPA + GABA + NMDA + STD + STF + asynchr. release
eqs_neuron_E_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = eqs_neuron_E_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI
eqs_neuron_I_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = eqs_neuron_I_ampa_gaba_nmda_stdE_stdI_asynchrE_asynchrI
eqs_synapse_E_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = '''
s_nmda_tot_post = w * S * x_d * u_d * s_nmda : 1 (summed)
qar_E_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
ds_nmda/dt = -s_nmda/(taus_nmda)+alpha_nmda*(x_nmda)*(1-s_nmda) + x0 * qar : 1 (clock-driven)
dx_nmda/dt = -x_nmda/(taux_nmda) : 1 (clock-driven)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
w : 1
S : 1
U : 1
'''
eqs_synapse_I_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = '''
qar_I_tot_post = w * x0 * qar :Hz (summed)
qar = clip(randn()*sqrt(x_d*u_d/x0*uar*dt*(1-uar*dt))+uar*dt*u_d*x_d/x0, 0, 2*x_d*u_d/x0*uar*dt)/dt :Hz  (scales correctly for different dt)
dx_d/dt = (1-x_d)/tau_d -qar : 1 (clock-driven)
du_d/dt = (U-u_d)/tau_f :1 (clock-driven)
duar/dt = -uar/tau_ar : Hz (clock-driven)
w : 1
U : 1
S : 1
'''
eqs_onpre_E_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = eqs_onpre_ampa_nmda_stdE_stfE_asynchrE
eqs_onpre_I_ampa_gaba_nmda_stdE_stdI_stfE_stfI_asynchrE_asynchrI = '''
x_d *= (1-u_d)
u_d += U*(1-u_d)
uar += U_ar*(U_max-uar)
s_gaba += w * S * x_d * u_d 
'''







#%%
###############################################################################

### Excitatory + inhibitory neuronal cultures, but GABAr blocked

# only AMPA + GABA=0



# AMPA + GABA=0 + NMDA



# AMPA + GABA=0 + NMDA + asynchr. release




#%%
###############################################################################

### Excitatory + inhibitory neuronal cultures, but AMPAr and NMDAr blocked

# only AMPA=0 + GABA



# AMPA =0 + GABA + NMDA=0



# AMPA=0 + GABA + NMDA=0 + asynchr. release










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
   
    ### Lookup tables: ### 
    E_neuron = {
                sort_tuple(("ampa",)): eqs_neuron_ampa, 
                sort_tuple(("ampa", "std")): eqs_neuron_ampa_std, 
                sort_tuple(("ampa", "std", "stf")): eqs_neuron_ampa_std_stf,
                sort_tuple(("ampa", "nmda")): eqs_neuron_ampa_nmda, 
                sort_tuple(("ampa", "nmda", "std")): eqs_neuron_ampa_nmda_std, 
                sort_tuple(("ampa", "nmda", "std", "stf")): eqs_neuron_ampa_nmda_std_stf, 
                sort_tuple(("ampa", "nmda", "std", "asynchr")): eqs_neuron_ampa_nmda_std_asynchr, 
                sort_tuple(("ampa", "nmda", "std", "stf", "asynchr")): eqs_neuron_ampa_nmda_std_stf_asynchr, 
                sort_tuple(("ampa", "gaba")): eqs_neuron_E_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_neuron_E_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_neuron_E_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_neuron_E_ampa_gaba_std_asynchr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_neuron_E_ampa_gaba_std_stf_asynchr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_neuron_E_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_neuron_E_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_neuron_E_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_neuron_E_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_neuron_E_ampa_gaba_nmda_std_stf_asynchr,                 
                }
    
    E_synapse = {
                sort_tuple(("ampa",)): eqs_synapse_ampa, 
                sort_tuple(("ampa", "std")): eqs_synapse_ampa_std, 
                sort_tuple(("ampa", "std", "stf")): eqs_synapse_ampa_std_stf,
                sort_tuple(("ampa", "nmda")): eqs_synapse_ampa_nmda,
                sort_tuple(("ampa", "nmda", "std")): eqs_synapse_ampa_nmda_std, 
                sort_tuple(("ampa", "nmda", "std", "stf")): eqs_synapse_ampa_nmda_std_stf, 
                sort_tuple(("ampa", "nmda", "std", "asynchr")): eqs_synapse_ampa_nmda_std_asynchr, 
                sort_tuple(("ampa", "nmda", "std", "stf", "asynchr")): eqs_synapse_ampa_nmda_std_stf_asynchr,
                sort_tuple(("ampa", "gaba")): eqs_synapse_E_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_synapse_E_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_synapse_E_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_synapse_E_ampa_gaba_std_asyncr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_synapse_E_ampa_gaba_std_stf_asyncr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_synapse_E_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_synapse_E_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_synapse_E_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_synapse_E_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_synapse_E_ampa_gaba_nmda_std_stf_asynchr, 
                
                }
    E_onpre = {
                sort_tuple(("ampa",)): eqs_onpre_ampa, 
                sort_tuple(("ampa", "std")): eqs_onpre_ampa_std, 
                sort_tuple(("ampa", "std", "stf")): eqs_onpre_ampa_std_stf,
                sort_tuple(("ampa", "nmda")): eqs_onpre_ampa_nmda,
                sort_tuple(("ampa", "nmda", "std")): eqs_onpre_ampa_nmda_std, 
                sort_tuple(("ampa", "nmda", "std", "stf")): eqs_onpre_ampa_nmda_std_stf, 
                sort_tuple(("ampa", "nmda", "std", "asynchr")): eqs_onpre_ampa_nmda_std_asynchr, 
                sort_tuple(("ampa", "nmda", "std", "stf", "asynchr")): eqs_onpre_ampa_nmda_std_stf_asynchr,
                sort_tuple(("ampa", "gaba")): eqs_onpre_E_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_onpre_E_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_onpre_E_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_onpre_E_ampa_gaba_std_asynchr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_onpre_E_ampa_gaba_std_stf_asynchr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_onpre_E_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_onpre_E_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_onpre_E_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_onpre_E_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_onpre_E_ampa_gaba_nmda_std_stf_asynchr, 
                
                }
    
    ###########################################################################
    
    I_neuron = {
                sort_tuple(("ampa", "gaba")): eqs_neuron_I_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_neuron_I_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_neuron_I_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_neuron_E_ampa_gaba_std_asynchr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_neuron_E_ampa_gaba_std_stf_asynchr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_neuron_I_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_neuron_I_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_neuron_I_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_neuron_I_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_neuron_I_ampa_gaba_nmda_std_stf_asynchr,
                
                }    
    I_synapse = {
                sort_tuple(("ampa", "gaba")): eqs_synapse_I_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_synapse_I_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_synapse_I_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_synapse_E_ampa_gaba_std_asyncr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_synapse_E_ampa_gaba_std_stf_asyncr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_synapse_I_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_synapse_I_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_synapse_I_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_synapse_I_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_synapse_I_ampa_gaba_nmda_std_stf_asynchr, 
                }
    
    I_onpre = {
                sort_tuple(("ampa", "gaba")): eqs_onpre_I_ampa_gaba, 
                sort_tuple(("ampa", "gaba", "std")): eqs_onpre_I_ampa_gaba_std, 
                sort_tuple(("ampa", "gaba", "std", "stf")): eqs_onpre_I_ampa_gaba_std_stf, 
                # sort_tuple(("ampa", "gaba", "std", "asynchr")): eqs_onpre_E_ampa_gaba_std_asynchr, # Not existing; CAVE
                # sort_tuple(("ampa", "gaba", "std", "stf", "asynchr")): eqs_onpre_E_ampa_gaba_std_stf_asynchr, # Not existing; CAVE
                sort_tuple(("ampa", "gaba", "nmda")): eqs_onpre_I_ampa_gaba_nmda, 
                sort_tuple(("ampa", "gaba", "nmda", "std")): eqs_onpre_I_ampa_gaba_nmda_std, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf")): eqs_onpre_I_ampa_gaba_nmda_std_stf, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "asynchr")): eqs_onpre_I_ampa_gaba_nmda_std_asynchr, 
                sort_tuple(("ampa", "gaba", "nmda", "std", "stf", "asynchr")): eqs_onpre_I_ampa_gaba_nmda_std_stf_asynchr, 
                }
    
    ###########################################################################
    
    ### Model components: ###
    incl_E_neurons = dict_model_config["excitatory_neurons"]   # this is a boolean
    incl_I_neurons = dict_model_config["inhibitory_neurons"]
   
    
    ### Construct equation-output-dict: ###
    mod_equations = {}
        
    if incl_E_neurons: 
        E_mechanisms = sort_tuple(dict_model_config["E_mechanisms"])
        try:
            mod_equations["E_neuron"] = E_neuron[E_mechanisms]
            mod_equations["E_synapse"] = E_synapse[E_mechanisms]
            mod_equations["E_onpre"] = E_onpre[E_mechanisms]
        except KeyError:
            raise ValueError(
                f"\nUnsupported/not readily implemented excitatory mechanisms combination: {E_mechanisms}.\n"                
                )

    
    if incl_I_neurons:
        I_mechanisms = sort_tuple(dict_model_config["I_mechanisms"])
        try:
            mod_equations["I_neuron"] = I_neuron[I_mechanisms]
            mod_equations["I_synapse"] = I_synapse[I_mechanisms]
            mod_equations["I_onpre"] = I_onpre[I_mechanisms]        
        except KeyError:
            raise ValueError(
                f"\nUnsupported/not readily implemented inhibitory mechanisms combination: {I_mechanisms}.\n"                
                )

    return mod_equations