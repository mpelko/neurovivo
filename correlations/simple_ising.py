"""
A set of functions for generating and analysing a simple max-ent model (which
can be written in Ising form with fixed J and h parameters for N neurons). 
For further understanding you can look at:

Tkacik, Gasper, Jason S. Prentice, Vijay Balasubramanian, and Elad Schneidman.
2010. "Optimal population coding by noisy spiking neurons." Proceedings of the
National Academy of Sciences 107(32): 14419 -14424.

or:

Macke, Jakob, Philipp Berens, Alexander Ecker, Andreas Tolias, and Matthias
Bethge. 2009. "Generating spike trains with specified correlation
coefficients." Neural computation 21(2): 397-423.
"""

import numpy as np

DEBUG = False

# -----------------------------------------------------------------------------
def generate_N_spike_trains(N, h, J, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains following the minimum model
    (maximum-entropy distribution), which can be written in the Ising form as:
    ...
    """
    
    if not np_rnd:
        np_rnd = np.random
    
    # Calculating the probabilities of P spikes being elicit in a population 
    # within one time step.
    P_probabilities = P_probability(N,h,J)
    if DEBUG:
        assert sum(P_probabilities) == 1.

    # Creating a range from 0 to 1 based on the P_probabilities.
    borders = np.zeros(N+2)
    current_sum = 0.
    for P, P_prob in enumerate(P_probabilities):
        current_sum += P_prob 
        borders[P+1] = current_sum 
    # Randomly assigning the number of Ps for the time bins, based on the
    # previously calculated probability profile.
    intermediate_randoms = np_rnd.rand(time)
    Ps_over_time = np.digitize(intermediate_randoms, borders) - 1
    # Randomly shuffeling those P spikes among the neurons in the population
    # for each time bin.
    result = np.zeros((N, time))
    for t in xrange(time):
        P = Ps_over_time[t]
        spike_pattern = np.concatenate([np.ones(P), np.zeros(N-P)])
        np_rnd.shuffle(spike_pattern)
        result[:,t] = spike_pattern

    return result

# -----------------------------------------------------------------------------
def P_probability(N,h,J):
# -----------------------------------------------------------------------------
    """
    Calculating the probabilities of P spikes being elicit in a population
    within one time step.
    """
    import gmpy2
    def top_f(N, P, h, J):
        h_mul = (2.*P)-N 
        #result = h*h_mul + J*(N*(N-1)/2 - 2*P*(N-P)) + J*N/2.
        result = h*h_mul + J*(N*(N-1)/2. - 2*P*(N-P))
        return result

    Z = sum([gmpy2.bincoef(N,P)*np.exp(top_f(N,P,h,J)) for P in xrange(N+1)])
    P_prob = [(gmpy2.bincoef(N,P)/Z)*np.exp(top_f(N,P,h,J)) for P in xrange(N+1)]
    return P_prob

# -----------------------------------------------------------------------------
def average(N,h,J):
# -----------------------------------------------------------------------------
    """
    Average of the spike train - probability that a spike is elicited in the
    time-bin (this related to spiking rate via the time-bin width).
    """
    P_prob = P_probability(N,h,J)
    return sum([(1.0*P/N)*P_prob[P] for P in xrange(N+1)])

# -----------------------------------------------------------------------------
def corr2(N, h, J):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of the 2nd order correlation.
    """
    P_prob = P_probability(N,h,J)
    avg_x = average(N,h,J)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_prob[P] for P in xrange(2,N+1)])
    result = (avg_xy - avg_x**2)
    return result 

# -----------------------------------------------------------------------------
def corr3(N, h, J):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of correlation of 3rd order.
    """
    P_prob = P_probability(N,h,J)
    avg_x = average(N,h,J)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_prob[P] for P in xrange(2,N+1)])
    avg_xyz = sum([1.0*P*(P-1)*(P-2)/(1.0*N*(N-1)*(N-2))*P_prob[P] for P in xrange(3,N+1)])
    result = avg_xyz - 3 * avg_xy*avg_x + 2 * avg_x**3
    
    return result 

# -----------------------------------------------------------------------------
def find_good_hJ_parameters(N, aver, corr):
# -----------------------------------------------------------------------------
    """
    Returns the h and J parameters for the desired "N" spike trains with
    "average" spiking ratio per time bin and average population pairwise
    correlation corr2
    """
    from scipy.optimize import fmin
    #from scipy.optimize import fmin_tnc
    #from scipy.optimize import anneal   # really buggy. Doesn't produce results.
    
    assert N>0

    av_des = aver     # desired average
    corr_des = corr   # desired correlation
    h_cand = -0.1     # candidate parameter h
    J_cand = 0.01     # candidate parameter J
    
    def opt_func(param):
        h = param[0]
        J = param[1]
        av_curr = average(N, h, J)
        av_diff = (av_curr - av_des) 
        corr_curr = corr2(N, h, J)
        corr_diff = (corr_curr - corr_des) 
        return abs(av_diff) + abs(corr_diff)
        #return abs(av_diff)/abs(av_des + 0.0000001) + abs(corr_diff)/abs(corr_des + 0.0000001)

    result = fmin(opt_func, [h_cand, J_cand])
    #result = fmin_tnc(opt_func, [h_cand, J_cand], approx_grad=True)[0]
    #result = anneal(opt_func, [h_cand, J_cand], maxeval = 200)[0]

    return result
