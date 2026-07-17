#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:17:26 2026

Utilities plot file. 

@author: Guido
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import networkx as nx
from brian2 import (second, ms, umeter, mV, pA)
from HH_utils.utils_analysis import Calciumtrace

#%% ### Plot functions: ###

def volt_plot(runsettings_dict, output_monitors): 
    trace_E = output_monitors.get("trace_E", None)
    trace_I = output_monitors.get("trace_I", None)
    plt.figure(dpi=200)
    if trace_E is not None:
        plt.plot(trace_E.t / second, trace_E[2].V / mV, 'k', linewidth=0.7)
    if trace_I is not None:
        plt.plot(trace_I.t/second, trace_I[2].V/mV, 'r', linewidth=0.7)
    plt.xlabel('time (s)')
    plt.ylabel('Membrane Potential of an individual neuron (mV)')
    
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    if trace_I is not None:
        plt.legend(['Excitatory', 'Inhibitory'])
    if runsettings_dict.get("save_figs", False):
        outputdir = runsettings_dict.get("output_dir", "/home/")
        noise_idx = runsettings_dict["noise_seed"]
        simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
        plt.savefig(outputdir + simname + '_VoltSingleNeuron.png')
    plt.show()
    

def raster_plot(runsettings_dict, output_monitors, N_E=0): 
    spikes_E = output_monitors.get("spikes_E", None)
    spikes_I = output_monitors.get("spikes_I", None)
    transient = runsettings_dict.get("sim_transient", 0 * second)
    simtime = runsettings_dict["sim_time"]
    plt.figure(dpi=200)
    if spikes_E is not None:
        plt.plot(spikes_E.t / second, spikes_E.i, '.k', ms=0.7)
    if spikes_I is not None:
        plt.plot(spikes_I.t / second, spikes_I.i + N_E, '.r', ms = 0.7)
    plt.title('Rasterplot per neuron (E-neurons: black, I-neurons: red)')
    plt.xlabel('time (s)')
    plt.ylabel('neuron index')
    plt.xlim([transient/second, simtime/second])
    
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    if runsettings_dict.get("save_figs", False):
        outputdir = runsettings_dict.get("output_dir", "/home/")
        noise_idx = runsettings_dict["noise_seed"]
        simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
        plt.savefig(outputdir + simname + '_RasterPlot.png')

    plt.show()


def electrode_plot(runsettings_dict, t, Voltagefilt):

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    
    offset = 100
   
    for k in range(12):
        ax.plot(
            t,
            Voltagefilt[k, :] + k * offset,
            linewidth=0.75,
            color="black"
        )
   
    
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("Electrode number", fontsize=12)
    ax.set_yticks(np.arange(12) * offset)
    ax.set_yticklabels(np.arange(12), fontsize=12)
   
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

   
    plt.tight_layout()
    if runsettings_dict.get("save_figs", False):
        outputdir = runsettings_dict.get("output_dir", "/home/")
        noise_idx = runsettings_dict["noise_seed"]
        simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
        plt.savefig(outputdir + simname + '_ElectrodePlot.png')
    plt.show()


def rasterlec_plot(runsettings_dict, t, APs, NBs=None, time_bin=25 * ms):   
    simtime = runsettings_dict["sim_time"]
    transient = runsettings_dict.get("sim_transient", 0 * second)
    dt2 = runsettings_dict["time_step"]
    
    
    fig, ax = plt.subplots(figsize=(10, 2.2), dpi=300)
    ax.plot(APs[:, 1] * dt2, APs[:, 0], 'k|', ms=9, alpha=.5)
    ax.hlines([-0.3, 0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7, 10.7],
        transient/second - 1, simtime/second + 1, colors='black', linewidths=0.5)
    
    if NBs is not None: 
        # Shade network bursts
        for nb in NBs:           
            start = transient/second + nb[0] * (time_bin/second)
            end   = transient/second + nb[1] * (time_bin/second)
            ax.axvspan(start, end, color='red', alpha=0.2)
    # Labels
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Electrode")
    
    # Electrode ticks
    ax.set_yticks(np.arange(12))
    ax.set_yticklabels(np.arange(1, 13))   # or np.arange(12) if you prefer 0–11
    
    # Limits
    ax.set_xlim([transient/second - 1, simtime/second + 1])
    ax.set_ylim([-0.5, 11.5])
    
    # Despine
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    if runsettings_dict.get("save_figs", False):
        outputdir = runsettings_dict.get("output_dir", "/home/")
        noise_idx = runsettings_dict["noise_seed"]
        simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
        fig.savefig(outputdir + simname + '_RasterLecPlot.png', dpi=300)
    plt.show()
        
# runsettings_dict, output_monitors, network_dict, excitatory_pop=None, inhibitory_pop=None, N_E=0, N_I=0): 
def topology_plot(runsettings_dict, network_dict, x_electrodes, y_electrodes, excitatory_pop=None, inhibitory_pop=None, N_E=0, N_I=0):
    N_tot = N_E  
    Adj = np.full((N_E, N_E), 0)
    Network = nx.from_numpy_array(Adj)  # make graph of adjacency matrix
    nodess = list(range(0, N_E))
    position = [(excitatory_pop.x[i] / umeter, excitatory_pop.y[i] / umeter) for i in range(N_E)]
    pos = {nodess[i]: position[i] for i in range(len(nodess))}  # make position argument
    
    # define inhibitory neuron network
    if inhibitory_pop is not None and N_I>0:
        N_tot += N_I
        Adj2 = np.full((N_I, N_I), 0)
        Network2 = nx.from_numpy_array(Adj2)  # make graph of adjacency matrix
        nodess2 = list(range(0, N_I))
        position2 = [(inhibitory_pop.x[i] / umeter, inhibitory_pop.y[i] / umeter) for i in range(N_I)]
        pos2 = {nodess2[i]: position2[i] for i in range(len(nodess2))}  # make position argument
    
    # define electrode network
    elecNetwork = nx.from_numpy_array(np.zeros((12, 12)))
    elecnodess = list(range(0, 12))
    elecposition = [(x_electrodes[i] / umeter, y_electrodes[i] / umeter) for i in range(12)]
    elecpos = {elecnodess[j]: elecposition[j] for j in range(12)}
    radius_elec = runsettings_dict["radius_detect_neurons"] / umeter
    
    plt.figure(figsize=(10, 10), dpi=200)
    ax = plt.gca()
    # draw electrodes    
    nx.draw(elecNetwork, elecpos, node_size=9000/N_tot, node_color='black', with_labels=False, font_color='whitesmoke', ax=ax)    
    for node, (x, y) in elecpos.items():
        circ = Circle((x, y), radius_elec, fill=True, color='grey', alpha=0.4, linewidth=2)
        ax.add_patch(circ)
        
        ax.text(x, y, str(int(node)+1),
            ha='center', va='center',
            fontsize=6,
            alpha=0.35,
            color='white')
        
    # add neurons
    nx.draw(Network, pos, node_size=500/N_tot, node_color='blue', with_labels=False, font_color='whitesmoke', alpha=0.7, ax=ax)
    if inhibitory_pop is not None:
        nx.draw(Network2, pos2, node_size=500/N_tot , node_color='chocolate', with_labels=False, alpha=0.7, ax=ax)
    _ = ax.axis('off')
    

    ax.set_aspect('equal')
    ax.set_title('Placing of excitatory neurons (blue), inhibitory neurons (orange), \n and electrodes (black) with measurement range (grey)', fontsize=12)
    
    if network_dict["neuron_position_type"] == "random": 
        # span_width_elec = np.max(x_electrodes / umeter) - np.min(x_electrodes / umeter)
        radius_neuron_placement =  runsettings_dict["radius_MEA_well"] / umeter 
        circ_placement = Circle((0, 0), radius_neuron_placement, fill=False, color='k', alpha=0.8, linewidth=2, linestyle='--')
        ax.add_patch(circ_placement)
        
    if runsettings_dict.get("save_figs", False):
        outputdir = runsettings_dict.get("output_dir", "/home/")
        noise_idx = runsettings_dict["noise_seed"]
        simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
        plt.savefig(outputdir + simname + '_TopologyPlot.png')
    
    plt.show()


        
def onechannel_plot(runsettings_dict, t, voltagetrace):
    for i in range(12):
        
        plt.figure(dpi=300)
        plt.plot(t, voltagetrace[i,:], 'k', linewidth=0.75)
        plt.xlabel('time (s)')
        plt.ylabel('Voltage (mV)')
        plt.xlim([5, 10])
        plt.ylim([-15, 30])
        plt.title(f"Voltage measured by electrode {i+1}")
        ax = plt.gca()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        
        if runsettings_dict.get("save_figs", False):
            outputdir = runsettings_dict.get("output_dir", "/home/")
            noise_idx = runsettings_dict["noise_seed"]
            simname = runsettings_dict.get("sim_name", f"simulation_{noise_idx}")
            plt.savefig(outputdir + simname + 'OneChannelPlot.png')       
        plt.show()



def STD_plot(runsettings_dict, output_monitors, model_components_dict, synapse_dict, N_E=0):
    
    spikes_E = output_monitors.get("spikes_E", None)
    spikes_I = output_monitors.get("spikes_I", None)
    
    traceconn_EE = output_monitors.get("traceconn_EE", None)
    
    fig, axs = plt.subplots(4, 1, sharex='all', dpi=200, figsize=[8, 6.5])
    
    if "stf" in model_components_dict["E_mechanisms"]:
        utrace = np.sum(traceconn_EE.u_d, axis = 0) / len(traceconn_EE.u_d)
    else:
        U_EE = synapse_dict["strength_depression_e2e"]
        utrace = np.full((len(traceconn_EE.x_d[1,:]),), U_EE)
        
    axs[0].plot(spikes_E.t / second, spikes_E.i, '.k', ms=1)
    axs[0].set_ylabel('Neuron index')
    plt.title('Short term plasticity in the E-E synapses')

    axs[1].plot(traceconn_EE.t / second, np.sum(traceconn_EE.x_d, axis=0) / len(traceconn_EE.x_d), 'k')
    axs[1].set_ylabel('x_d')
    axs[1].grid(visible=True, linestyle=':')

    axs[2].plot(traceconn_EE.t / second, utrace, 'k')
    axs[2].set_ylabel('u_d')
    axs[2].grid(visible=True, linestyle=':')

    axs[3].plot(traceconn_EE.t / second, (np.sum(traceconn_EE.x_d, axis=0) / len(traceconn_EE.x_d))*utrace, 'k')
    axs[3].set_ylabel('x*u')
    axs[3].set_xlabel('time (s)')
    axs[3].grid(visible=True, linestyle=':')
        
    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    fig.tight_layout()
    # savefig(outputdir + simname + 'STDplotEE.png', dpi=300)
    plt.show()

    if model_components_dict["inhibitory_neurons"]:
        traceconn_EI = output_monitors.get("traceconn_EI", None)
        traceconn_IE = output_monitors.get("traceconn_IE", None)
        traceconn_II = output_monitors.get("traceconn_II", None)
        U_EI = synapse_dict["strength_depression_e2i"]
        U_IE = synapse_dict["strength_depression_i2e"]
        U_II = synapse_dict["strength_depression_i2i"]
        
        fig, axs = plt.subplots(4, 1, sharex='all', dpi=200, figsize=[6, 4])
        axs[0].plot(spikes_E.t / second, spikes_E.i, '.k', ms=1)
        axs[0].plot(spikes_I.t / second, spikes_I.i + N_E, '.r', ms=1)
        if "stf" in model_components_dict["E_mechanisms"]:
            utrace = np.sum(traceconn_EI.u_d, axis = 0) / len(traceconn_EI.u_d)
        else:            
            utrace = np.full((len(traceconn_EI.x_d[1,:]),), U_EI)
            
        axs[0].set_ylabel('Neuron index')
        plt.title('Short term plasticity in the E-I synapses')

        axs[1].plot(traceconn_EI.t / second, np.sum(traceconn_EI.x_d, axis=0) / len(traceconn_EI.x_d), 'k')
        axs[1].set_ylabel('x_d')
        axs[1].grid(visible=True, linestyle=':')

        axs[2].plot(traceconn_EI.t / second, utrace, 'k')
        axs[2].set_ylabel('u_d')
        axs[2].grid(visible=True, linestyle=':')

        axs[3].plot(traceconn_EI.t / second, (np.sum(traceconn_EI.x_d, axis=0) / len(traceconn_EI.x_d)) * utrace, 'k')
        axs[3].set_ylabel('x*u')
        axs[3].set_xlabel('time (s)')
        axs[3].grid(visible=True, linestyle=':')
        fig.tight_layout()
        # savefig(outputdir + simname + 'STDplotEI.png', dpi=300)
        plt.show()

        fig, axs = plt.subplots(4, 1, sharex='all', dpi=200, figsize=[6,4])
        axs[0].plot(spikes_E.t / second, spikes_E.i, '.k', ms=1)
        axs[0].plot(spikes_I.t / second, spikes_I.i + N_E, '.r', ms=1)
        if "stf" in model_components_dict["I_mechanisms"]:
            utrace = np.sum(traceconn_IE.u_d, axis = 0) / len(traceconn_IE.u_d)
        else:
            utrace = np.full((len(traceconn_IE.x_d[1,:]),), U_IE)
        axs[0].set_ylabel('Neuron index')
        plt.title('Short term plasticity in the I-E synapses')

        axs[1].plot(traceconn_IE.t / second, np.sum(traceconn_IE.x_d, axis=0) / len(traceconn_IE.x_d), 'k')
        axs[1].set_ylabel('x_d')
        axs[1].grid(visible=True, linestyle=':')

        axs[2].plot(traceconn_IE.t / second, utrace, 'k')
        axs[2].set_ylabel('u_d')
        axs[2].grid(visible=True, linestyle=':')

        axs[3].plot(traceconn_IE.t / second, (np.sum(traceconn_IE.x_d, axis=0) / len(traceconn_IE.x_d)) * utrace, 'k')
        axs[3].set_ylabel('x*u')
        axs[3].set_xlabel('time (s)')
        axs[3].grid(visible=True, linestyle=':')
        fig.tight_layout()
        # savefig(outputdir + simname + 'STDplotIE.png', dpi=300)
        plt.show()

        fig, axs = plt.subplots(4, 1, sharex='all', dpi=200, figsize=[6,4])
        axs[0].plot(spikes_E.t / second, spikes_E.i, '.k', ms=1)
        axs[0].plot(spikes_I.t / second, spikes_I.i + N_E, '.r', ms=1)
        if "stf" in model_components_dict["I_mechanisms"]:
            utrace = np.sum(traceconn_II.u_d, axis = 0) / len(traceconn_II.u_d)
        else:
            utrace = np.full((len(traceconn_II.x_d[1,:]),), U_II)
        axs[0].set_ylabel('Neuron index')
        plt.title('Short term plasticity in the I-I synapses')

        axs[1].plot(traceconn_II.t / second, np.sum(traceconn_II.x_d, axis=0) / len(traceconn_II.x_d), 'k')
        axs[1].set_ylabel('x_d')
        axs[1].grid(visible=True, linestyle=':')

        axs[2].plot(traceconn_II.t / second, utrace, 'k')
        axs[2].set_ylabel('u_d')
        axs[2].grid(visible=True, linestyle=':')

        axs[3].plot(traceconn_II.t / second, (np.sum(traceconn_II.x_d, axis=0) / len(traceconn_II.x_d)) * utrace, 'k')
        axs[3].set_ylabel('x*u')
        axs[3].set_xlabel('time (s)')
        axs[3].grid(visible=True, linestyle=':')
        fig.tight_layout()
        # savefig(outputdir + simname + 'STDplotII.png', dpi=300)
        plt.show()



def synapse_plot(runsettings_dict, output_monitors, model_components_dict):
    trace_E = output_monitors.get("trace_E", None)
    trace_I = output_monitors.get("trace_I", None)
    E_mechanisms = model_components_dict.get("E_mechanisms", [])
    I_mechanisms = model_components_dict.get("I_mechanisms", [])
    
    currents = []

    if "nmda" in E_mechanisms and hasattr(trace_E, "I_nmda"):
        currents.append(("NMDA", trace_E.I_nmda))    
    if "ampa" in E_mechanisms and hasattr(trace_E, "I_ampa"):
        currents.append(("AMPA", trace_E.I_ampa))    
    if "gaba" in E_mechanisms and hasattr(trace_E, "I_gaba"):
        currents.append(("GABA", trace_E.I_gaba))
    currents.append(("Total", trace_E.I_syn))
    
  
    nplots = len(currents)
    fig, axs = plt.subplots(nplots, 1, sharex=True, dpi=200, figsize=[6, 1*nplots])
    if nplots == 1:
        axs = [axs]
    t = trace_E.t / second
    computed_currents = []
    
    for i, (label, data) in enumerate(currents):           
        current = np.mean(data, axis=0) * -1 / pA
        computed_currents.append(current)
        axs[i].plot(t, current, 'k')
        axs[i].set_ylabel(label)
        axs[i].grid(True, linestyle=':')

    if trace_I is not None:
        currentsI = []
        if "nmda" in I_mechanisms and hasattr(trace_I, "I_nmda"):
            currentsI.append(("NMDA", trace_I.I_nmda))    
        if "ampa" in E_mechanisms and hasattr(trace_E, "I_ampa"):
            currentsI.append(("AMPA", trace_I.I_ampa))    
        if "gaba" in E_mechanisms and hasattr(trace_E, "I_gaba"):
            currentsI.append(("GABA", trace_I.I_gaba))
        currentsI.append(("Total", trace_I.I_syn))        
        tI = trace_I.t / second


        for i, (label, data) in enumerate(currents):
     
            currentI = np.mean(data, axis=0) * -1 / pA
            axs[i].plot(tI, currentI, 'r')     
        
        axs[-1].legend(['Excitatory', 'Inhibitory'])
    
    axs[-1].set_xlabel("time (s)")
    plt.title("Average Synaptic currents")
    axs.spines["top"].set_visible(False)
    axs.spines["right"].set_visible(False)
    fig.tight_layout()
    # savefig(outputdir + simname + 'SynapticCurrents.png', dpi=300)
    plt.show()

def adaptation_plot(runsettings_dict, output_monitors):
    trace_E = output_monitors.get("trace_E", None)
    trace_I = output_monitors.get("trace_I", None)
    plt.figure()
    plt.plot(trace_E.t / ms, np.sum(trace_E.I_AHP, axis=0) / len(trace_E.I_AHP)/pA, 'k')
    if trace_I is not None:
        plt.plot(trace_I.t / ms, np.sum(trace_I.I_AHP, axis=0) / len(trace_I.I_AHP) / pA, 'r')
        plt.legend(['Excitatory', 'Inhibitory'])
    plt.title('Average after-hyperpolarizing current')
    plt.xlabel('t (ms)')
    plt.ylabel('I_AHP (pA)')
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # savefig(outputdir + simname + 'Adaptation.png', dpi=300)
    plt.show()



def mechanism_plot(runsettings_dict, output_monitors, N_E=0, N_I=0):
            
    trace_E = output_monitors.get("trace_E", None)
    trace_I = output_monitors.get("trace_I", None)
    spikes_E = output_monitors.get("spikes_E", None)
    spikes_I = output_monitors.get("spikes_I", None)
    
    transient = runsettings_dict.get("sim_transient", 0 * second)
    dt2 = runsettings_dict["time_step"]    
    simtime = runsettings_dict["sim_time"]
    fs = 1 / (dt2 / second)
    spikerateE, _ = Calciumtrace(np.array(spikes_E.t), simtime, transient, fs)
    
    
    fig, axs = plt.subplots(4, 1, sharex='all', dpi=200, figsize=[5, 4])
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    if trace_I is not None:
        I_syn_I_mean = np.mean(trace_I.I_syn/ pA, axis=0)
        axs[0].plot(trace_I.t / second, I_syn_I_mean, 'slategray')
    I_syn_E_mean = np.mean(trace_E.I_syn/ pA, axis=0)
    axs[0].plot(trace_E.t / second, I_syn_E_mean, '#D2691E')
    axs[0].set_ylabel('PSC (pA)')
    axs[0].set_ylim([-1, 300])

    if trace_I is not None:
        spikerateI, _ = Calciumtrace(np.array(spikes_I.t), simtime, transient, fs)
        axs[1].plot(np.arange(len(spikerateI))*25*ms, spikerateI/N_I, 'slategray')
    axs[1].plot(np.arange(len(spikerateE))*25*ms, spikerateE/N_E, '#D2691E')
    axs[1].set_ylabel('FR (Hz)')
    axs[1].set_ylim([-1, 100])
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
       
    if trace_I is not None:
        traceconn_EI = output_monitors.get("traceconn_EI", None)
        axs[2].plot(traceconn_EI.t / second, np.sum(traceconn_EI.x_d, axis=0) / len(traceconn_EI.x_d), 'slategray')
    traceconn_EE = output_monitors.get("traceconn_EE", None)
    axs[2].plot(traceconn_EE.t / second, np.sum(traceconn_EE.x_d, axis=0) / len(traceconn_EE.x_d), '#D2691E')
    axs[2].set_ylabel('STD')
    axs[2].set_ylim([0.4, 1.1])
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
       
    if trace_I is not None:
        axs[3].plot(trace_I.t / second, np.sum(trace_I.I_AHP, axis=0) / len(trace_I.I_AHP) * -1 / pA, 'slategray')
    axs[3].plot(trace_E.t / second, np.sum(trace_E.I_AHP, axis=0) / len(trace_E.I_AHP) * -1 / pA, '#D2691E')
    axs[3].set_ylabel('AHP (PA)')
    axs[3].set_xlabel('time (s)')
    axs[3].set_ylim([-1, 35])
    axs[3].spines['top'].set_visible(False)
    axs[3].spines['right'].set_visible(False)
    if trace_I is not None:
        axs[0].legend(['I', 'E'])
    fig.tight_layout()
    plt.xlim([transient / second, (simtime - transient) / second])
    plt.title("Trans-membrane I, firing rate, STD, AHP Plot", fontsize=12)
    # savefig(outputdir + simname + 'mechplotzoom.png', dpi=300)
    plt.show()
    




def spikerate_plot(runsettings_dict, output_monitors):
    
    spikes_E = output_monitors.get("spikes_E", None)
    spikes_I = output_monitors.get("spikes_I", None)
    
    transient = runsettings_dict.get("sim_transient", 0 * second)
    dt2 = runsettings_dict["time_step"]    
    simtime = runsettings_dict["sim_time"]
    fs = 1 / (dt2 / second)
    
    if spikes_I is not None:
        APs = np.append(spikes_E.t, spikes_I.t)
    else:
        APs = np.array(spikes_E.t)
    spikerate, _ = Calciumtrace(APs, simtime, transient, fs)
    
    plt.figure(dpi=300, figsize=(6, 3))
    plt.plot(np.arange(len(spikerate))*25*ms, spikerate, 'k')
    plt.xlabel('time (s)', fontsize=12)
    plt.ylabel('Spikerate (spikes/second)', fontsize=12)
    
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #plt.xlim([11,20])
    # savefig(outputdir + simname + 'spikerate.png')
    plt.show()

# def calciumplot(runsettings_dict, output_monitors, N_E=0, N_I=0):
     
#     spikes_E = output_monitors.get("spikes_E", None)
#     spikes_I = output_monitors.get("spikes_I", None)
    
#     transient = runsettings_dict.get("sim_transient", 0 * second)
#     dt2 = runsettings_dict["time_step"]    
#     simtime = runsettings_dict["sim_time"]
#     fs = 1 / (dt2 / second)
    
#     if spikes_I is not None:
#         APs = np.append(spikes_E.t, spikes_I.t)
#     else:
#         APs = np.array(spikes_E.t)
#     spikerate, _ = Calciumtrace(APs, simtime, transient, fs)
    
#     fig, axs = plt.subplots(2, 1, sharex='all', dpi=200, figsize=[8, 6.5])
#     axs[0].plot(spikes_E.t / second, spikes_E.i, '.k', ms=0.7)
#     if spikes_I is not None:
#         axs[0].plot(spikes_I.t / second, spikes_I.i + N_E, '.r', ms=0.7)
#     axs[0].set_ylabel('Neuron index')
#     axs[0].set_xlim([transient/second, simtime/second])
#     catime = (np.arange(len(spikerate))*25*ms+transient)/second
#     axs[1].plot(catime, Catrace/(N_E+N_I)/1000, color='maroon')
#     #axs[1].scatter(catime[peakinds], Catrace[peakinds]/10000)
#     axs[1].set_ylabel('Calcium signal intensity')
#     axs[1].set_xlabel('Time (s)')
#     #axs[1].grid(visible=True, linestyle=':')
#     axs[1].set_ylim([-0.1, 1.5])
#     axs[1].set_xlim([transient/second, simtime/second])
#     fig.tight_layout()
#     # savefig(outputdir + simname + 'rasterenCa.png', dpi=300)
#     plt.show()

def calcium_plot(runsettings_dict, output_monitors, N_E=0, N_I=0):
     
    spikes_E = output_monitors.get("spikes_E", None)
    spikes_I = output_monitors.get("spikes_I", None)
    
    transient = runsettings_dict.get("sim_transient", 0 * second)
    dt2 = runsettings_dict["time_step"]    
    simtime = runsettings_dict["sim_time"]
    fs = 1 / (dt2 / second)
    
    if spikes_I is not None:
        APs = np.append(spikes_E.t, spikes_I.t)
    else:
        APs = np.array(spikes_E.t)
    spikerate, Catrace = Calciumtrace(APs, simtime, transient, fs)
    fig, axs = plt.subplots(2, 1, sharex='all', dpi=200, figsize=[8, 6.5])
    axs[0].plot(spikes_E.t / second, spikes_E.i, '.', mfc='#D2691E', mec='#D2691E', ms=1)
    
    if N_I > 0:
        axs[0].plot(spikes_I.t / second, spikes_I.i + N_E, '.', mfc='slategray', mec='slategray', ms=1)
    axs[0].set_ylabel('Neuron index')
    axs[0].set_xlim([transient / second, simtime / second])
    axs[0].axis('off')
    catime = (np.arange(len(spikerate)) * 25 * ms + transient) / second
    axs[1].plot(catime, Catrace / (N_E + N_I) / 1000, color='#69140E', linewidth=2.5)
    axs[1].axis('off')
    # axs[1].scatter(catime[peakinds], Catrace[peakinds]/10000)
    axs[1].set_ylabel('Calcium signal intensity')
    axs[1].set_xlabel('Time (s)')
    # axs[1].grid(visible=True, linestyle=':')
    axs[1].set_ylim([-0.1, 1.5])
    axs[1].set_xlim([transient / second, simtime / second])
    fig.tight_layout()
    # savefig(outputdir + simname + 'rasterenCa.png', dpi=300)
    plt.show()

def patchplot(runsettings_dict, output_monitors, N_E=0, N_I=0):
     
    trace_E = output_monitors.get("trace_E", None)
    
    I_syn = trace_E[1].I_syn * -1 / pA
    
    plt.figure(dpi=200)
    plt.plot(trace_E[1].t/second, I_syn, 'k')
    plt.xlabel('time (s)')
    plt.ylabel('Synaptic current')
    plt.xlim([10.5, 13])
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    plt.show()


#%% ### Select Plots and construct them: ###


def get_plots_HH(simulation_dict, output_monitors, model_components_dict, 
                 electrode_info_dict, excitatory_pop, inhibitory_pop): 
    
    runsettings_dict = simulation_dict["runsettings_dict"]
    network_dict = simulation_dict["network_dict"]    
    plot_dict = simulation_dict["plot_dict"]
    analysis_dict = simulation_dict.get("analysis_dict", None)
    synapse_dict = simulation_dict["synapse_dict"]
    
    N_E=network_dict.get("n_e_neurons", 0)    
    N_I=network_dict.get("n_i_neurons", 0)
    
    
    
    
    if plot_dict.get("neuron_voltplot"): 
        volt_plot(runsettings_dict, output_monitors)
    if plot_dict.get("neuron_rasterplot"):
        raster_plot(runsettings_dict, output_monitors, N_E=N_E)   
    if plot_dict.get("STDplot", None):
        STD_plot(runsettings_dict, output_monitors, model_components_dict, synapse_dict, N_E=0)        
    if plot_dict.get("synapseplot", None):
        synapse_plot(runsettings_dict, output_monitors, model_components_dict)
    if plot_dict.get("adaptationplot", None):
        if "AHP" in model_components_dict["E_mechanisms"]:
            adaptation_plot(runsettings_dict, output_monitors)
    if plot_dict.get("summary_mechanismplot", None):
        if "AHP" in model_components_dict["E_mechanisms"]:
            mechanism_plot(runsettings_dict, output_monitors, N_E=N_E, N_I=N_I)
    if plot_dict.get("neuron_spikerateplot", None):
        spikerate_plot(runsettings_dict, output_monitors)
    if plot_dict.get("calciumplot", None):
        calcium_plot(runsettings_dict, output_monitors, N_E=N_E, N_I=N_I)
    if plot_dict.get("patchplot", None): 
        patchplot(runsettings_dict, output_monitors, N_E=N_E, N_I=N_I)
        
         
    electrodeplot= plot_dict.get("electrode_voltplot",False)
    rasterlecplot = plot_dict.get("electrode_rasterplot", False)
    topologyplot = plot_dict.get("topologyplot", False)
    onechannelplot = plot_dict.get("one_electrode_voltplot", False)
     
    if electrodeplot or rasterlecplot or topologyplot or onechannelplot:
         
        t = output_monitors["trace_E"].t / second
        x_electrodes = electrode_info_dict["x_electrodes"]
        y_electrodes = electrode_info_dict["y_electrodes"]
        APs = electrode_info_dict["APs"]
        voltagetraces = electrode_info_dict["voltagetraces"]
        NBs = electrode_info_dict.get("NBs", None)        
        if electrodeplot: 
            electrode_plot(runsettings_dict, t, voltagetraces)
        if rasterlecplot: 
            time_bin = analysis_dict["analysis_args"].get("time_bin", 25 * ms)
            rasterlec_plot(runsettings_dict, t, APs, NBs=NBs, time_bin=time_bin)
        if topologyplot: 
            topology_plot(runsettings_dict, network_dict, x_electrodes, y_electrodes, excitatory_pop=excitatory_pop, 
            inhibitory_pop=inhibitory_pop, N_E=N_E, N_I=N_I)
        if onechannelplot: 
            onechannel_plot(runsettings_dict, t, voltagetraces)
            
    return














