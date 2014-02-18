"""
A set of functions for generating and analysing correlated spike trains based
on the procedure described by Staude et al. in: 

Staude B, Grun S, Rotter S. Higher-order correlations in non-stationary
parallel spike trains: statistical modeling and inference. Front Comput
Neurosci. 2010;4. Available at: http://www.ncbi.nlm.nih.gov/pubmed/20725510.
"""

import numpy as np
import correlated_input as ci

DEBUG = False

# -----------------------------------------------------------------------------
def generate_N_spike_trains(N, aver,  P_probs_m, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains following he Staude process,
    where the "aver" is average spikes per bin (has to be < 1) and P_probs_m is
    the probability distribution of having P spikes in the population at the
    same time-bin (for P > 0).
    """
    
    if not np_rnd:
        np_rnd = np.random

    verify_parameters(N, aver, P_probs_m)
    P_probs_m = np.array(P_probs_m)/sum(P_probs_m)
        
    mean_P = np.sum([P_probs_m[i]*(i+1) for i in xrange(len(P_probs_m))])
    aver_g = N * aver / mean_P
    #print mean_P, aver_g
    
    if aver_g == 0:
        grandfather_st = np.zeros(time)
    else:
        ISI = int(1. / aver_g)
        grandfather_st = ci.poisson_train(ISI, time, np_rnd)
    
    trains = np.zeros((N,time))
    
    # Creating a range from 0 to 1 based on the P_prob.
    borders = np.zeros(len(P_probs_m)+1)
    current_sum = 0.
    for P, P_prob in enumerate(P_probs_m):
        current_sum += P_prob 
        borders[P+1] = current_sum 
    intermediate_randoms = np_rnd.rand(sum(grandfather_st))
    Ps_when_gf_spike = np.digitize(intermediate_randoms, borders)

    gf_spike_index = 0    
    for t, bin in enumerate(grandfather_st):
        if bin == 1:
            P = Ps_when_gf_spike[gf_spike_index]
            spike_pattern = np.concatenate([np.ones(P), np.zeros(N-P)])
            np_rnd.shuffle(spike_pattern)
            trains[:,t] = spike_pattern
            gf_spike_index += 1
    
    return trains

# -----------------------------------------------------------------------------
def average(N,aver,P_probs_m):
# -----------------------------------------------------------------------------
    """
    Average of the spike train - probability that a spike is elicited in the
    time-bin (this related to spiking rate via the time-bin width).
    """
    verify_parameters(N, aver, P_probs_m)
    
    P_probs = get_P_probs(N, aver, P_probs_m)
    
    return sum([(1.0*P/N)*P_probs[P] for P in xrange(len(P_probs))])

# -----------------------------------------------------------------------------
def corr2(N, aver, P_probs_m):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of the 2nd order correlation.
    """
    verify_parameters(N, aver, P_probs_m)
    
    P_probs = get_P_probs(N, aver, P_probs_m)
    
    avg_x = average(N,aver,P_probs_m)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_probs[P]\
                  for P in xrange(2,len(P_probs))])
    result = (avg_xy - avg_x**2)
    return result 

# -----------------------------------------------------------------------------
def corr3(N, aver, P_probs_m):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of correlation of 3rd order.
    """
    verify_parameters(N, aver, P_probs_m)

    P_probs = get_P_probs(N, aver, P_probs_m)

    avg_x = average(N,aver,P_probs_m)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_probs[P]\
                  for P in xrange(2,len(P_probs))])
    avg_xyz = sum([1.0*P*(P-1)*(P-2)/(1.0*N*(N-1)*(N-2))*P_probs[P]\
                   for P in xrange(3,len(P_probs))])
    result = avg_xyz - 3 * avg_xy*avg_x + 2 * avg_x**3
    
    return result 

# -----------------------------------------------------------------------------
def verify_parameters(N, aver,  P_probs_m):
# -----------------------------------------------------------------------------
    """
    Checks if the combination of parameters is realistic.
    """
    
    assert len(P_probs_m) <= N, "You can only specify probabilities for P up to" +\
                            "N in the P_probs_m."
    assert np.all(np.array(P_probs_m)>=0), "P_probs_m positive."
    P_probs_m = np.array(P_probs_m)/sum(P_probs_m)

    mean_P = np.sum([P_probs_m[i]*(i+1) for i in xrange(len(P_probs_m))])
    aver_g = N * aver / mean_P

    if aver_g == 0:
        return      # This would be an empty spike train, so parameters are ok.

    assert aver_g <= 1, "The desired parameters would require too high" +\
                 "grandfather train density. Reduce N, aver or change P_probs_m"

    ISI = int(1. / aver_g)
    aver_g_real = 1. / ISI
    
    assert abs(aver_g_real - aver_g) < 0.7, "aver_g is too destorted. Choose "+\
                                            "different parameters.\n" +\
                                            "Current parameters: " +\
            "N = {0}, aver = {1}, P_probs_m = {2}".format(N, aver, P_probs_m)                                                    

# -----------------------------------------------------------------------------
def get_P_probs(N, aver, P_probs_m):
# -----------------------------------------------------------------------------
    """
    Taking the parameters it returns P_probs - the probability distribution for
    all Ps, including P=0 
    """

    if aver == 0:
        result = np.zeros(len(P_probs_m) + 1)
        result[0] = 1.
        return result

    mean_P = np.sum([P_probs_m[i]*(i+1) for i in xrange(len(P_probs_m))])
    aver_g = N * aver / mean_P
    
    ISI = int(1. / aver_g)

    P_0 = 1.0*(ISI-1)/ISI
    divider = 1./(1-P_0)
    P_probs = np.concatenate([[P_0],np.array(P_probs_m)/divider])
    return P_probs
