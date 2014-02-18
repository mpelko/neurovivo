import quantities as pq
import numpy as np
from neurovivo.common import Trace
import neurovivo.common as cmn

# -----------------------------------------------------------------------------
def spike_times(trace, treshold=0*pq.mV):
# -----------------------------------------------------------------------------
    """
    Checks for the spike times. Spike time is the time, when the voltage
    crosses over the treshold value.
    This is the slow but sure version.
    """
    data = trace._data
    assert data.units == pq.mV,\
           "trace for counting spikes must be in mV units"
    if len(data) <= 0:
        print "WARNING: Counting spikes in trace of lenght 0."
        return []*pq.ms 
    spike_times = []
    spike_flag = data[0] < treshold
    for i, voltage in enumerate(data):
        if spike_flag and voltage > treshold:
            spike_times.append(trace._time[i]) 
            spike_flag = False
        if not spike_flag and voltage < treshold:
            spike_flag = True

    return spike_times

# -----------------------------------------------------------------------------
def spike_times_fast(trace, treshold=0*pq.mV):
# -----------------------------------------------------------------------------
    """
    Checks for the spike times. Spike time is the time, when the voltage
    crosses over the treshold value.
    This is the fast version.
    """
    data = trace._data
    assert data.units == pq.mV,\
           "trace for counting spikes must be in mV units"
    if len(data) <= 0:
        print "WARNING: Counting spikes in trace of lenght 0."
        return []*pq.ms 
    data_bool = data > treshold
    spike_times = trace._time[data_bool - np.roll(data_bool,1)]
    useless_first = 0
    if data_bool[0]: 
        useless_first += 1
    if data_bool[0] != data_bool[-1]:
        useless_first += 1
    return spike_times[useless_first::2]

# -----------------------------------------------------------------------------
def spike_count(trace, treshold=0*pq.mV):
# -----------------------------------------------------------------------------
    """
    Returns the number of spikes in the trace.
    """
    return len(spike_times_fast(trace, treshold))

# -----------------------------------------------------------------------------
def was_spike_elicited(trace, window_range=None, treshold=0*pq.mV):
# -----------------------------------------------------------------------------
    """
    Checks if the spike was elicited during the time in the  window range (in 
    the form of (start_time, end_time]). Treshold dictates the definition of
    the spike.
    """
    if len(trace._data) == 0:
        print "WARNING: the checking if spike elicited on the empty trace."
        return False
    
    # If no window is given, we check the whole trace
    if not window_range:
        window_range = (trace._time[0], trace._time[-1])
    
    spikes = spike_times_fast(trace, treshold)
    
    sp_in_window = spikes[np.all([spikes>window_range[0],
                                  spikes<=window_range[1]], axis=0)]

    return len(sp_in_window) > 0

# -----------------------------------------------------------------------------
def probability_profile(trace, bins=20):
# -----------------------------------------------------------------------------
    """
    Returns the voltage probability profile.
    """
    data = trace._data
    return np.histogram(data,bins=bins,density=True)

# -----------------------------------------------------------------------------
def spike_times_simple(strace, threshold=0):
# -----------------------------------------------------------------------------
    """
    Checks for the spike times. Spike time is the time, when the voltage
    crosses over the treshold value.
    This is the fast version.
    """
    data = strace._data
    if len(data) <= 0:
        print "WARNING: Counting spikes in trace of lenght 0."
        return [] 
    data_bool = data > threshold
    time = np.arange(len(data))*strace._dt
    spike_times = time[data_bool - np.roll(data_bool,1)]
    useless_first = 0
    if data_bool[0]: 
        useless_first += 1
    if data_bool[0] != data_bool[-1]:
        useless_first += 1
    return spike_times[useless_first::2]

# -----------------------------------------------------------------------------
def spike_count_simple(strace, treshold=0):
# -----------------------------------------------------------------------------
    """
    Returns the number of spikes in the trace.
    """
    return len(spike_times_simple(strace, treshold))

# -----------------------------------------------------------------------------
def trace_moments_simple(strace):
# -----------------------------------------------------------------------------
    """
    Returns the first 4 moments of voltage trace distribution.
    The spikes are not cut out from the trace.
    """
    from scipy import stats
    data = strace._data
    return [np.mean(data), np.std(data)**2, stats.skew(data), stats.kurtosis(data)]

# -----------------------------------------------------------------------------
def remove_spikes_from_trace(trace, threshold=-40*pq.mV):
# -----------------------------------------------------------------------------
    """
    Returns the trace object whith cut spikes replaced by threshold values.
    """
    data = trace._data[:]
    data[data>threshold]=threshold
    trace._data = data

# -----------------------------------------------------------------------------
def get_big_events(trace, threshold=None, min_amplitude=2*pq.mV, time_limit=3*pq.ms):
# -----------------------------------------------------------------------------
    """
    Returns the big events according to the Paolo analysis.
    """
    import scipy.io as sio
    import os
    
    data = trace._data[:]
    time = trace._time[:]
    res = {"voltage":data,"times":time, "min_amplitude":min_amplitude, "timelimit":time_limit, "analysis_folder":cmn.HOME + "/ws/neurovivo/analysis/"}
    temp_dir = cmn.HOME + "/tmp/"
    cmn.mkdir_p(temp_dir)
    full_path = temp_dir+"tmp_matlab_data_{}.mat".format(np.random.randint(0,100000))
    sio.savemat(full_path, res)
    os.system('octave -f -q ' + cmn.HOME + "/ws/neurovivo/analysis/get_big_events.m " + full_path)
    big_event = sio.loadmat(full_path)
    print big_event
    big_event_amplitudes = big_event["result"][0]
    big_event_start_bins = big_event["result"][1]
    cmn.rm(full_path)
    return big_event_amplitudes * pq.mV, big_event_start_bins 

# -----------------------------------------------------------------------------
def power_spectrum_old(trace):
# -----------------------------------------------------------------------------
    '''
    Returns the power spectrum of trace
    '''
    trace = trace-trace.mean()
    times = trace._time
    times = times.rescale(pq.s)
    transform = np.abs(trace.fft())**2
    transform = transform[:len(times)/2]
    freqs = np.fft.fftfreq(len(trace._data), times[1]-times[0])[:len(times)/2]
    return freqs, transform
    
# -----------------------------------------------------------------------------
def power_spectrum(trace, time_window=1000, overlap_window=None, enforce_pow2=True, freq_range=None):
# -----------------------------------------------------------------------------
    '''
    Returns the power spectrum of trace. Parameters set to fit with what Paolo 
    does within his project.
    '''
    from matplotlib.mlab import psd
    from matplotlib.pylab import detrend_linear, detrend_mean, detrend_none
    from neurovivo.common import nextpow2
        
    dt = (trace._time[1]-trace._time[0]).rescale(pq.ms).item()
    Fs = 1000./dt
    NFFT=int(time_window/dt)
    if enforce_pow2:
        NFFT=nextpow2(NFFT)
        time_window = 1.*time_window*NFFT/int(time_window/dt)
    if overlap_window == None:
        overlap_window = 0.75 * time_window
    noverlap=int(overlap_window/dt)
    window = np.bartlett(NFFT)
    #print dt, Fs, NFFT, noverlap
    pows, freqs = psd(detrend_mean(trace._data), NFFT=NFFT, Fs=Fs, noverlap=noverlap, window=window, detrend=detrend_none)
    if freq_range == None:
        return pows, freqs
    else:
        freqs, pows = cmn.splice_two_vectors(freqs, pows, freq_range)
        return pows, freqs

# -----------------------------------------------------------------------------
def specgram(trace, NFFT=256, noverlap=128):
# -----------------------------------------------------------------------------
    '''
    Returns the spectrogram of the trace
    '''
    from matplotlib.pylab import specgram
    from matplotlib.pylab import detrend_mean
    dt = (trace._time[1]-trace._time[0]).item()
    Fs = 1000./dt
    return specgram(trace._data, NFFT=NFFT, Fs=Fs, detrend=detrend_mean, noverlap=noverlap)

# -----------------------------------------------------------------------------
def full_analysis(strace, cutoff=50, threshold=-45,):
# -----------------------------------------------------------------------------
    """
    strace - voltage trace of class SimpleTrace
    cutoff - first part of the trace cut away from the analysis (in ms)
    threshold - cutting voltage
    Returns [spike_rate, mean, var, skew, kurtosis]
    """
    from scipy import stats
    from common import SimpleTrace
    voltage = strace._data
    dt = strace._dt
    cut_trace = voltage[int(cutoff/dt):]
    rate = 1000.*spike_count_simple(SimpleTrace(dt, cut_trace))/(dt*len(cut_trace))
    data = cut_trace[cut_trace<threshold]
    return rate, np.mean(data), np.std(data)**2, stats.skew(data), stats.kurtosis(data)

# -----------------------------------------------------------------------------
def graphical_analysis(strace, comment=None, cutoff=50, threshold=-45, bins=50):
# -----------------------------------------------------------------------------
    """
    Graphical report of the trace.
    strace - voltage trace of class SimpleTrace
    cutoff - first part of the trace cut away from the analysis (in ms)
    threshold - cutting voltage
    """
    import matplotlib.pylab as pyl
    voltage = strace._data
    time = np.arange(len(voltage))*strace._dt
    pyl.figure()
    rate, mean, var, skew, kurt = full_analysis(strace, cutoff, threshold)
    if comment:
        pyl.suptitle(comment)
    sp1 = pyl.subplot(2,1,1)
    sp1.plot(time,voltage)
    sp1.set_title("Spike rate = {0}".format(rate))
    sp1.set_xlabel("time [ms]")
    sp1.set_ylabel("V [mV]")
    sp2 = pyl.subplot(2,1,2)
    cut_trace = voltage[int(cutoff/strace._dt):]
    data = cut_trace[cut_trace<threshold]
    sp2.hist(data, bins=bins, histtype="stepfilled", normed=1)
    xlim = sp2.get_xlim()
    pyl.text(xlim[0]+0.7*(xlim[1]-xlim[0]), 0.6,
             "mean = {0}\nvar={1}\nskew={2}\nkurt={3}".format(mean, var, skew, kurt))
    sp2.plot([mean, mean],[0,1],"r")
    sp2.plot([mean-np.sqrt(var)/2., mean+np.sqrt(var)/2.],[0.1,0.1], "r")
    sp2.set_xlabel("V [mV]")
    sp2.set_ylabel("normalized distribution")
    
