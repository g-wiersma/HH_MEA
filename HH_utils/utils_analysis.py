#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:56:22 2026

Utilities file for MEA data analysis. 

@author: Guido
"""

import numpy as np
import pandas as pd
from brian2 import ms, second, umeter, mV, meter
import scipy
from scipy.signal import find_peaks
from scipy.fft import fft
from scipy.stats import norm
from itertools import combinations
import matplotlib.pyplot as plt



def get_electrode_info(runsettings_dict, output_monitors, excitatory_pop=None, inhibitory_pop=None, N_E=0, N_I=0): 
    trace_E = output_monitors.get("trace_E", None)
    trace_I = output_monitors.get("trace_I", None)
    
    dt2 = runsettings_dict["time_step"]
    # set up a filter to filter the voltage signal
    fs = 1 / (dt2 / second)
    fc = 100  # Cut-off frequency of the filter
    w = fc / (fs / 2)  # Normalize the frequency
    b, a = scipy.signal.butter(2, w, 'high')
    voltagetraces = np.zeros((12, len(trace_E.t)))

    # determine from which neurons the electrodes measure a signal (faster than measuring everything)
    elec_grid_dist = runsettings_dict["electrode_grid_distance"]    # electrode grid size (there are 12 electrodes)
    elec_range = 3 * elec_grid_dist                                      # total width/heigth that the electrode grid spans
    elec_offset = elec_range / 2
   

    # Find the excitatory neurons that need to be measures by each electrode
    elecrangesE = np.full((16, 40), np.nan)
    elecrangesI = np.full((16, 40), np.nan)
    
    r_detect_neuron = runsettings_dict["radius_detect_neurons"] # the radius in which an electrode measures neurons
    r2 = r_detect_neuron ** 2 
    
    xE = excitatory_pop.x
    yE = excitatory_pop.y
    
    if inhibitory_pop:
        xI = inhibitory_pop.x
        yI = inhibitory_pop.y
        
    x_electrodes = []
    y_electrodes = []
    
    for i in [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14]: # note: of the 4x4 grid, the corner electrodes are skipped.
    
        x_electrode = (i % 4) * elec_grid_dist  - elec_offset
        y_electrode = (i // 4) * elec_grid_dist - elec_offset
        x_electrodes.append(x_electrode)
        y_electrodes.append(y_electrode)
        
        dx = xE - x_electrode
        dy = yE - y_electrode
        mask = dx**2 + dy**2 < r2
    
        idx = np.where(mask)[0] # indeces of which neurons are measured by the electrode
        n = min(len(idx), 40) # make sure these are not more than 40 
        elecrangesE[i, :n] = idx[:n]
    
        if inhibitory_pop is not None:
            dx = xI - x_electrode
            dy = yI - y_electrode
            mask = dx**2 + dy**2 < r2
    
            idx = np.where(mask)[0]
            n = min(len(idx), 40) # make sure these are not more than 40 
            elecrangesI[i, :n] = idx[:n]
   

    # compute the signal every electrodes measures and filter it
    k = 0
    APs = []
    
    # convert neuron positions to NumPy arrays for efficiency
    xE_positions = excitatory_pop.x[:]
    yE_positions = excitatory_pop.y[:]
    if inhibitory_pop is not None:
        xI_positions = inhibitory_pop.x[:]
        yI_positions = inhibitory_pop.y[:]
    
    for k,i in enumerate([1,2,4,5,6,7,8,9,10,11,13,14]):  # skip corners
        x_electrode = x_electrodes[k]
        y_electrode = y_electrodes[k]
    
        # get valid neuron indices for this electrode
        templist = elecrangesE[i, ~np.isnan(elecrangesE[i, :])].astype(int)
        Voltage = np.zeros(len(trace_E.V[0]))
                
        # excitatory contribution
        dx = xE_positions[templist] - x_electrode
        dy = yE_positions[templist] - y_electrode
        distances = np.sqrt(dx**2 + dy**2) / meter

        d_cutoff = np.min([r_detect_neuron * 0.2, 5* umeter])
        scaling_weights = np.ones_like(distances)        
        far_field_mask = distances > d_cutoff
        scaling_weights[far_field_mask] = d_cutoff / distances[far_field_mask]
        
        Voltage += np.sum((trace_E[templist].V / mV) * scaling_weights[:, np.newaxis],  axis=0)
        
        # inhibitory contribution
        if inhibitory_pop is not None:
            templist2 = elecrangesI[i, ~np.isnan(elecrangesI[i, :])].astype(int)
            dx = xI_positions[templist2] - x_electrode
            dy = yI_positions[templist2] - y_electrode
            distances = np.sqrt(dx**2 + dy**2) / meter
            scaling_weights = np.ones_like(distances)        
            far_field_mask = distances > d_cutoff
            scaling_weights[far_field_mask] = d_cutoff / distances[far_field_mask]
            
            Voltage += np.sum((trace_E[templist2].V / mV) * scaling_weights[:, np.newaxis],  axis=0)
        
    
        # high-pass filter
        Voltage = Voltage - np.mean(Voltage)
        Voltagefilt = scipy.signal.filtfilt(b, a, Voltage)
        voltagetraces[k, :] = Voltagefilt
    
        # detect APs
        threshold = 4 * np.sqrt(np.mean(Voltagefilt**2))
        APstemp, _ = find_peaks(np.abs(Voltagefilt), height=threshold)
        APs.extend([(k, t / second) for t in APstemp])
    
        k += 1
    
    APs = np.array(APs)
    electrode_info_dict = {
        "x_electrodes": x_electrodes, 
        "y_electrodes": y_electrodes, 
        "voltagetraces": voltagetraces, 
        "APs": APs, 
        "elecrangesE": elecrangesE
            }
    
    if inhibitory_pop is not None:
        electrode_info_dict["elecrangesI"] = elecrangesI
       
    return electrode_info_dict

def Calciumtrace(APs, simtime, transient, fs, timeBin = 25 * ms, tau = 2 * second):
   
    #initialize
    APsinbin = np.zeros((int(np.floor((simtime - transient) / timeBin))))
    spikerate = APsinbin

    # delete transient
    APs_wot = APs[APs > transient / second]
    APs_wot = APs_wot - transient / second

    APt = APs_wot
    binsize = timeBin / second
    for k in range(int(np.floor((simtime - transient) / timeBin))):
        APsinbin[k] = sum((APt > (k * binsize)) & (APt < ((k + 1) * binsize)))
        spikerate[k] = APsinbin[k] * second / timeBin

    tau = tau / timeBin
    Catrace = np.convolve(spikerate, np.exp(-np.arange((simtime-transient)/timeBin)/tau),'full')
    Catrace = Catrace[:len(spikerate)]
    
    return spikerate, Catrace


#%% ### Feature Functions: ###

def compute_spikerate(analysis_args, APs, rec_time, transient, fs):
    """Compute the spike rate by binning APs"""
    # delete transient
    APs_wot = APs[APs[:, 1] > transient* fs, :]
    APs_wot[:, 1] = APs_wot[:, 1] - transient * fs

    # bin spikes
    numelectrodes = analysis_args.get("numelectrodes", 12)
    time_bin = analysis_args["time_bin"] / second
    bins = int((rec_time - transient) / time_bin)
    binsize = time_bin * fs

    APsinbin = np.zeros((numelectrodes, int(np.floor((rec_time - transient) / time_bin))))
    spikerate = np.zeros((numelectrodes, int(np.floor((rec_time - transient) / time_bin))))

    for l in range(numelectrodes):
        APt = APs_wot[APs_wot[:, 0] == l, :]  # go through one electrode first
        APt = APt[:, 1]                      # take only the timestamps
        for k in range(int(np.floor((rec_time - transient) / time_bin))):
            APsinbin[l, k] = sum((APt > (k * binsize)) & (APt < ((k + 1) * binsize)))
            spikerate[l, k] = APsinbin[l, k] / time_bin

    APsinbintot = np.sum(APsinbin, axis=0)
    spikeratetot = np.sum(spikerate, axis=0)

    return APsinbin, spikerate, APsinbintot, spikeratetot

def smooth_spikerate(analysis_args, spikeratetot):
    """Apply Gaussian smoothing to the spike rate."""
    width, sigma = analysis_args["smoothing"].values()
    x = np.arange(width) - width // 2
    kernel = norm.pdf(x, scale=sigma)
    kernel /= kernel.sum()
    return np.convolve(spikeratetot, kernel, mode='same')

def detect_bursts(analysis_args, spikeratesmooth, APsinbin, spikerate, APsinbintot):
    """Detect network bursts."""
    args_burst_detec = analysis_args["burst_detection"]
    numelectrodes = analysis_args.get("numelectrodes", 12)

    start_th = args_burst_detec["Startth_factor"] * max(spikeratesmooth)
    time_th = int((args_burst_detec["time_th"] / second) / (analysis_args["time_bin"] / second) )
    actelec = np.sum(np.mean(spikerate, axis=1) > args_burst_detec["Act_thres"])    # calculate the number of active electrodes
    elec_th = args_burst_detec["Eth_factor"] * actelec                        # determine the amount of electrodes that need to be active within a burst
    stop_th = args_burst_detec["Stopth_factor"] * max(spikeratesmooth)

    # detect fragmented bursts as peaks of the smoothed spikerate
    MBth = args_burst_detec["FBth_factor"] * max(spikeratesmooth)
    peaks, ph = find_peaks(spikeratesmooth, height=MBth, prominence=args_burst_detec["FBprom_factor"] * max(spikeratesmooth))

    # initialize
    i = 0
    NBcount = 0
    maxNB = 1000
    NBs = np.zeros((maxNB, 4))

    # detect
    while (i + time_th) < len(spikeratesmooth):
        if (all(spikeratesmooth[i:i + time_th] > start_th)) \
                & (sum(np.sum(APsinbin[:, i:i + time_th], axis=1) > time_th) > elec_th):
            NBs[NBcount, 2] = NBs[NBcount, 2] + sum(APsinbintot[i:i + time_th])
            NBs[NBcount, 0] = i
            i = i + time_th
            while any(spikeratesmooth[i:i + 2] > stop_th):
                NBs[NBcount, 2] = NBs[NBcount, 2] + APsinbintot[i]
                i = i + 1
            NBs[NBcount, 3] = sum((peaks > NBs[NBcount, 0]) & (peaks < i))
            NBs[NBcount, 1] = i
            NBcount = NBcount + 1
        else:
            i = i + 1

    NBs = NBs[0:NBcount, :]

    return NBs

def partial_smooth_histogram(hist, window_size=3, boundary_bins=5):
    """Partially smooth the inter-burst-interval histogram while preserving the first few bins."""
    smoothed_hist = hist.copy()
    smoothed_hist[boundary_bins:] = np.convolve(hist[boundary_bins:], np.ones(window_size) / window_size,
                                                mode='same')
    return smoothed_hist

def detect_megabursts(analysis_args, NBs, plot=False):
    """Detect clusters of NBs and group them based on the histogram of inter-burst-intervals"""
    # Step 1: Calculate inter-NB-intervals in seconds
    lowerbound, upperbound, resetvalue = analysis_args["burst_detection"]["megaburst"].values()
    time_bin = analysis_args["time_bin"] / second
    start_times = NBs[:, 0] * time_bin  # Convert indices to time in seconds
    inter_NB_intervals = (np.array(NBs[1:, 0]) - np.array(NBs[0:-1, 1])) * time_bin

    # get histogram of inter-NB-intervals to detect grouped NBs
    hist, bin_edges = np.histogram(inter_NB_intervals, bins=50)
    hist = np.concatenate((np.array([0]), hist))
    smoothed_hist = partial_smooth_histogram(hist)

    if plot:
        plt.figure(dpi=300)
        plt.plot(bin_edges, hist, label='Original histogram', color='gray', alpha=0.5)
        plt.plot(bin_edges, smoothed_hist, label='Smoothed histogram', color='blue')
        plt.xlabel('Inter-NB interval (s)')
        plt.ylabel('Frequency')
        plt.title('Histogram of Inter-NB intervals')
        

    # Find a threshold for clustering bursts based on the peaks in the INBI histogram
    min_peak_height = np.max(smoothed_hist) * 0.1  # Set the minimum peak height as a percentage of the max
    peaks, _ = find_peaks(smoothed_hist, height=min_peak_height)

    if len(peaks) >= 2:
        # Find the index range between the first two significant peaks
        peak1_idx = peaks[0]
        peak2_idx = peaks[1]

        if bin_edges[peak1_idx-1] > 0.8:
            Tth = resetvalue
        else:
            # Slice the histogram between these two peaks and find the minimum
            valley_region = smoothed_hist[peak1_idx:peak2_idx]
            min_value = np.min(valley_region)
            minimal_indices = np.where(valley_region == min_value)[0] + peak1_idx

            # Choose the middle of the range of minimal values as Tth
            middle_min_idx = minimal_indices[len(minimal_indices) // 2]
            Tth = bin_edges[middle_min_idx]

            # Ensure Tth is within bounds
            Tth = min(max(Tth, lowerbound), upperbound)
    else:
        # If no clear bimodal distribution, use resetvalue
        Tth = resetvalue

    if plot:
        
        plt.axvline(Tth, color='red', linestyle='--', label=f'Threshold = {Tth:.3f} s')
        plt.legend()
        plt.show()

    # Step 4: Cluster NBs into megabursts
    megaNBs = []
    current_start = NBs[0, 0]  # Start time of the first megaburst
    current_maxima = NBs[0, 3]
    num_APs = NBs[0, 2]

    for i in range(1, len(NBs)):
        interval = (NBs[i, 0] - NBs[i - 1, 1]) * time_bin  # Calculate the interval in seconds
        if interval <= Tth:
            # Add to current megaburst
            current_maxima += NBs[i, 3]
            num_APs += NBs[i, 2]
        else:
            # End the current megaburst and start a new one
            current_stop = NBs[i - 1, 1]
            megaNBs.append([current_start, current_stop, num_APs, current_maxima])
            current_start = NBs[i, 0]
            current_maxima = NBs[i, 3]
            num_APs = NBs[i, 2]

    # Add the last megaburst
    current_stop = NBs[-1, 1]
    megaNBs.append([current_start, current_stop, num_APs, current_maxima])

    return np.array(megaNBs)

def compute_MAC(spikeratesmooth):
    """Compute MAC as defined by Maheswaranathan."""
    yf = fft(spikeratesmooth)
    return max(abs(yf[1:])) / abs(yf[0])

def rasterlecplot(analysis_args, APs, NBs, dt, plotstart, plotstop, color1, color2):
    """ Generate a raster plot with shaded regions indicating detected bursts """
    time_bin = analysis_args["time_bin"] / second
    fig, ax = plt.subplots(figsize=(8.5, 2.2), dpi=400)

    # Plot spikes as raster plot
    ax.plot(APs[:, 1] * dt, APs[:, 0], '|', markeredgecolor=color1, ms=6, alpha=.3)

    # Shade network bursts
    for nb in NBs:
        ax.axvspan(nb[0] * time_bin, nb[1] * time_bin, color=color2, alpha=0.2,
                   label="Network Burst" if 'NB' not in locals() else "")

    # Add grid lines to separate electrodes
    for i in range(12):
        ax.hlines(i + 0.5, -3, plotstop+5, colors='black', linewidths=0.7)

    ax.set_xlim([plotstart, plotstop])
    ax.set_ylim([-0.5, 12])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Electrode")
    plt.tight_layout()
    plt.show()


#%% ### Combined ISI and bursting features: ###

#%% ### Combined functions: ##
def compute_ISI_measures(analysis_dict, APs_wot, fs, rec_time, transient):
    """Compute ISI-based measures."""
    numelectrodes = analysis_dict.get("numelectrodes", 12)

    # Compute continuous ISI arrays
    time_vector = np.arange(0, (rec_time - transient), 1/fs)
    isi_arrays = np.zeros((numelectrodes, len(time_vector)))

    for electrode in range(numelectrodes):
        # Extract spike times for the current electrode
        electrode_spike_times = APs_wot[APs_wot[:, 0] == electrode, 1].astype(int)

        for i in range(len(electrode_spike_times) - 1):
            spike1 = electrode_spike_times[i]
            spike2 = electrode_spike_times[i + 1]
            tisi = (spike2 - spike1) / fs

            # Fill ISI values in the appropriate range
            if i == 0:
                isi_arrays[electrode, 0:spike1] = np.nan
                
            isi_arrays[electrode, spike1:spike2] = tisi
            
            if (i + 1) == (len(electrode_spike_times) - 1):
                isi_arrays[electrode, spike2:] = np.nan

    # Compute ISI measures
    meanisi_array = np.nanmean(isi_arrays, axis=0)
    meanISI = np.nanmean(meanisi_array)
    sdmeanISI = np.nanstd(meanisi_array)
    sdtimeISI = np.nanmean(np.nanstd(isi_arrays, axis=0))

    # Calculate the ISI-distance and ISI correlations
    all_combinations = list(combinations(list(np.arange(numelectrodes)), 2))
    isi_distances = np.zeros(len(all_combinations))
    isicoefficients = np.zeros(len(all_combinations))
    N = len(isi_arrays[0, :])
    j = 0
    # iterate through the electrode combinations
    for electrode1_key, electrode2_key in all_combinations:
        # Get the ISI arrays for the selected electrodes
        isit1_wn = isi_arrays[electrode1_key, :]
        isit2_wn = isi_arrays[electrode2_key, :]
        isit1 = isit1_wn[~np.isnan(isit1_wn)]
        isit2_wn = isit2_wn[~np.isnan(isit1_wn)]
        isit2 = isit2_wn[~np.isnan(isit2_wn)]
        isit1 = isit1[~np.isnan(isit2_wn)]

        isi_diff = isit1 / isit2
        isi_distances[j] = np.mean(np.where(isi_diff <= 1, abs(isi_diff - 1), -1 * (1 / isi_diff - 1)))

        if (i != j) & (not list(isit1) == list(isit2)):
            isicoefficients[j] = ((N * sum(isit1 * isit2) - sum(isit1) * (sum(isit2)))
                                  * ((N * sum(isit1 ** 2) - sum(isit1) ** 2) ** (-0.5))
                                  * ((N * sum(isit2 ** 2) - sum(isit2) ** 2) ** (-0.5)))

        j += 1

    meanisicorr = np.mean(isicoefficients)
    sdisicorr = np.std(isicoefficients)
    isi_distance = np.mean(isi_distances)
    return meanisicorr, sdisicorr, isi_distance, meanISI, sdmeanISI, sdtimeISI

def compute_bursting_features(analysis_dict, NBs, numAPs, rec_time, transient):
    time_bin = analysis_dict["analysis_args"]["time_bin"] / second
    numelectrodes = analysis_dict.get("numelectrodes", 12)
    NBcount = len(NBs[:,0])
    MNBR = NBcount * 60 / (rec_time - transient)        # Mean Network Burst Rate
    NBdurations = (np.array(NBs[:, 1]) - np.array(NBs[:, 0])) * time_bin
    MNBD = np.mean(NBdurations)                            # Mean Network Burst Duration
    PSIB = sum(NBs[:, 2] / numAPs) * 100                # Percentage of spike in Network Burst
    MFR = numAPs / (rec_time - transient)  / numelectrodes
    IBI = (np.array(NBs[1:, 0]) - np.array(NBs[0:-1, 1])) * time_bin
    CVIBI = np.std(IBI) / np.mean(IBI)                  # Coefficient of Variation of the inter-burst-intervals
    if NBcount == 0:
        MNBD = 0.0
        MNMBs = 0.0
        NFBs = 0
    else:
        NFBs = sum(NBs[:, 3]) / NBcount

    if NBcount < 2:
        CVIBI = 0.0

    return MFR, MNBR, MNBD, PSIB, NFBs, CVIBI


#%% ### Collect features for HH simulator: ###

def get_features_HH(analysis_dict, APs, rec_time, transient, fs, return_NBs=False, plot=False, plotstart=0, plotstop=50):
    # Compute the spikerate and smooth it
    analysis_args = analysis_dict["analysis_args"]
    APsinbin, spikerate, APsinbintot, spikeratetot = compute_spikerate(analysis_args, APs, rec_time, transient, fs)
    spikeratesmooth =  smooth_spikerate(analysis_args, spikeratetot)

    # Detect network bursts and group in megaburts (this replaces the original bursts with the megabursts if megaNB=true)
    NBs = detect_bursts(analysis_args, spikeratesmooth, APsinbin, spikerate, APsinbintot)
    if analysis_dict.get("calc_burst_features", False):
        if analysis_dict["analysis_args"].get("burst_detection", False) is not False:
            do_detect_megabursts = analysis_dict["analysis_args"]["burst_detection"].get("detect_mega_network_bursts", True)
    if len(NBs[:,1]) > 0 and do_detect_megabursts:
        NBs = detect_megabursts(analysis_args, NBs, plot=False)
    
    features = {}
    # Compute the network burst features based on either the detected bursts or detected megabursts (analysis_dict, NBs, numAPs, rec_time, transient)
    if analysis_dict.get("calc_burst_features", True):
        MFR, MNBR, MNBD, PSIB, NFBs, CVIBI = compute_bursting_features(analysis_dict, NBs, sum(APs[:, 1] > transient * fs), rec_time, transient)
        features.update({
            "MFR": MFR,
            "MNBR": MNBR,
            "MNBD": MNBD,
            "PSIB": PSIB,
            "#FBs": NFBs,
            "CV_INBI": CVIBI
        })
        
    # Compute ISI related measured using all APs after the transient 
    if analysis_dict.get("calc_spike_features", True):
        meanisicorr, sdisicorr, isi_distance, meanISI, sdmeanISI, sdtimeISI = compute_ISI_measures(analysis_dict, APs[APs[:, 1] > transient * fs, :], fs, rec_time, transient)
        MAC = compute_MAC(spikeratesmooth)      # Maximum autocorrelation component.
        features.update({
           "mean_ISI_corr": meanisicorr,
           "std_ISI_corr": sdisicorr,
           "ISI_distance": isi_distance,
           "mean_ISI": meanISI,
           "std_mean_ISI": sdmeanISI,
           "std_ISI_times": sdtimeISI,
           "MAC": MAC
       })


    features_df = pd.DataFrame([features])
    if return_NBs: 
        return NBs, features_df
    else:
        return features_df













