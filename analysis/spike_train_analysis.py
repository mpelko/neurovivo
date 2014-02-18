'''
Created on Dec 29, 2011

@author: mpelko

A collection of methods for the analysis of a spike train.
All the times are in ms.
'''

import numpy as np
import random
import operator
from neurovivo.common import SimpleTrace
from neurovivo.common import mymath

def correlation(st1, st2, dt=1, window=20):
    """
    The pairwise correlation of two spike trains.
    st1, st2 - the two spike trains
    dt - the time resolution used
    window - the window for averaging used
    """

    total_time, start_time = st1._total_time, st1._start_time
    assert (total_time - start_time) > window, "Too big time window."
    for st in [st1, st2]:
        assert st._total_time == total_time and st._start_time == start_time
    
    # take the dynamic rates
    dr1 = st1.dynamic_rate(dt=dt, window=window)
    dr2 = st2.dynamic_rate(dt=dt, window=window)
    
    # normalize
    dr1_norm = dr1 - np.mean(dr1)
    dr2_norm = dr2 - np.mean(dr2)
    
    return np.sum(dr1_norm*dr2_norm)*dt/(total_time-start_time)

def pearson_correlation(st1, st2, dt=1, window=20):
    """
    The pearson correlation of two spike trains.
    st1, st2 - the two spike trains
    dt - the time resolution used
    window - the window for averaging used
    """

    # only take non-empty spike trains
    s_trains = []
    for st in [st1,st2]:
        if len(st._spike_times) > 0:
            s_trains.append(st)
    sts=s_trains
    
    if len(sts) < 2:
        return np.nan

    total_time, start_time = st1._total_time, st1._start_time
    assert (total_time - start_time) > window, "Too big time window."
    for st in [st1, st2]:
        assert st._total_time == total_time and st._start_time == start_time
    
    # take the dynamic rates
    dr1 = st1.dynamic_rate(dt=dt, window=window)
    dr2 = st2.dynamic_rate(dt=dt, window=window)
    
    # normalize
    dr1_norm = dr1 - np.mean(dr1)
    dr2_norm = dr2 - np.mean(dr2)
    
    return np.sum(dr1_norm*dr2_norm)*dt/(total_time-start_time)/(np.std(dr1_norm)*np.std(dr2_norm))

def spike_count_correlation(st1, st2, dt=1, window=20):
    """
    The spike count correlation of two spike trains.
    st1, st2 - the two spike trains
    dt - the time resolution used
    window - the window for averaging used
    """

    # only take non-empty spike trains
    s_trains = []
    for st in [st1,st2]:
        if len(st._spike_times) > 0:
            s_trains.append(st)
    sts=s_trains
    
    if len(sts) < 2:
        return np.nan

    total_time, start_time = st1._total_time, st1._start_time
    assert (total_time - start_time) > window, "Too big time window."
    for st in [st1, st2]:
        assert st._total_time == total_time and st._start_time == start_time
    
    # take the dynamic rates
    dsc1 = st1.dynamic_spike_count(dt=dt, window=window)
    dsc2 = st2.dynamic_spike_count(dt=dt, window=window)
    
    # normalize
    dsc1_norm = dsc1 - np.mean(dsc1)
    dsc2_norm = dsc2 - np.mean(dsc2)
    
    return np.sum(dsc1_norm*dsc2_norm)*dt/(total_time-start_time)

def average_pairwise_correlation_full(sts, dt=1, window=20):
    """
    Returns the average pairwise_correlation among the set of spike trains sts.
    """
    
    # only take non-empty spike trains
    s_trains = []
    for st in sts:
        if len(st._spike_times) > 0:
            s_trains.append(st)
    sts=s_trains
    
    if len(sts) < 2:
        return np.nan
    
    drs = []
    total_time, start_time = sts[0]._total_time, sts[0]._start_time
    assert (total_time - start_time) > window, "Too big time window."
    for st in sts:
        assert st._total_time == total_time and st._start_time == start_time
        dr = st.dynamic_rate(dt=dt, window=window)
        dr_norm = dr - np.mean(dr)
        drs.append(dr_norm)
    N = len(drs)
    total_corr = 0
    for i in xrange(N):
        for j in xrange(i+1,N):
            total_corr += np.sum(drs[i]*drs[j])
            
    return (2*total_corr/(N*(N-1)))*dt/(total_time-start_time)

def average_3rd_order_correlation_full(sts, dt=1, window=20):
    """
    Returns the average 3rd order correlation among the set of spike trains sts.
    """
    # only take non-empty spike trains
    s_trains = []
    for st in sts:
        if len(st._spike_times) > 0:
            s_trains.append(st)
    sts=s_trains

    if len(sts) < 3:
        return np.nan
    
    drs = []
    total_time, start_time = sts[0]._total_time, sts[0]._start_time
    assert (total_time - start_time) > window, "Too big time window."
    for st in sts:
        assert st._total_time == total_time and st._start_time == start_time
        dr = st.dynamic_rate(dt=dt, window=window)
        dr_norm = dr - np.mean(dr)
        drs.append(dr_norm)
    N = len(drs)
    total_corr = 0
    n = 0
    for i in xrange(N):
        for j in xrange(i+1,N):
                for k in xrange(j+1,N):
                    total_corr += np.sum(drs[i]*drs[j]*drs[k])
                    n+=1
    return (total_corr/n)*dt/(total_time-start_time)

def nth_order_correlation(sts, nth_order, dt=1, windows=[20], nr_combinations=1000, statistics=True, corr_type="normalized"):
    """
    Returns the n-th order correlation among the set of spike trains 
    sts using nr_combinations numbers of combinations.
    If statistics = True, it returns the [mean and std] of the results, otherwise it
    returns the results for all nr_combinations.
    The correlation can be either evaluated on dynamic rates or on dynamic spike counts or on z-score (). 
    """
    
    # only take non-empty spike trains
    s_trains = []
    for st in sts:
        if len(st._spike_times) > 0:
            s_trains.append(st)
    sts=s_trains
    
    assert len(sts) >= nth_order, "more spike trains needed for {}th order".format(nth_order)
    
    total_time, start_time = sts[0]._total_time, sts[0]._start_time

    chosen_drs_inds = []
    for comb in xrange(nr_combinations):
        #print comb
        chosen_drs_inds.append(random.sample(xrange(len(sts)), nth_order))
    
    total_corrs = []
    for w in windows:
        assert (total_time - start_time) > w, "Too big time window."    

        values = []
        for i, st in enumerate(sts):
            #print "train {0} of {1}".format(i+1, len(sts))
            assert st._total_time == total_time and st._start_time == start_time
            if corr_type=="dynamic_rate" or corr_type=="normalized":
                value = st.dynamic_rate(dt=dt, window=w)
            elif corr_type=="spike_count":
                value = st.dynamic_spike_count(dt=dt, window=w)
            else:
                assert False, "bad corr_type"
            value_norm = value - np.mean(value)
            if corr_type=="normalized":
                value_norm = value_norm/np.std(value_norm)
            values.append(value_norm)
        total_corr = []
        
        for comb in xrange(nr_combinations):
            #print comb
            total_corr.append(np.sum(reduce(np.multiply, [values[ind] for ind in chosen_drs_inds[comb]])))
        total_corr = np.array(total_corr)*dt/(total_time-start_time)
        total_corrs.append(total_corr)
    if statistics:
        return [np.mean(total_corrs,1), np.std(total_corrs,1)]
    else:
        return total_corrs
    
def convolve_with_expexp_kernel(st, tau1, tau2, dt=0.1, kernel_length=250):
    """
    Going form spike trains to conductance changes in time by convolving the 
    spike train with the exponential-exponential kernel.
    st -            SpikeTrain
    tau1, tau2 -    time constants for the kernel
    dt -            sampling time in ms (should be smaller then tau1 and tau2)
    kernel_length - the length of kernel in ms
    
    Returns a SimpleTrace (so the info on start_time is lost).
    """
    assert dt<tau1 and dt<tau2, "choose smaller dt for the taus."
    spike_times = st._spike_times
    spike_times = spike_times-st._start_time
    discretised_spikes = np.histogram(spike_times, int((st._total_time-st._start_time)/dt), (0, st._total_time-st._start_time))[0]
    kernel = mymath.expexp_fun(np.arange(int(kernel_length/dt)), tau1, tau2)
    data = np.convolve(discretised_spikes, kernel)[:-(int(kernel_length/dt)-1)]
    
    return SimpleTrace(dt, data, name="convolved spike train")