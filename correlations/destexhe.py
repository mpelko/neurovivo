"""
Generating and investigating the spike trains as suggested by Destexhe in:

Destexhe, Alain, and Denis Pare. 1999. "Impact of Network Activity on the
Integrative Properties of Neocortical Pyramidal Neurons In Vivo." J
Neurophysiol 81(4): 1531-1547.
"""

import numpy as np
from itertools import chain

DEBUG = False
    
# -----------------------------------------------------------------------------
def create_n_trains(n, n2_trains, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Makes n trains of spikes timings based on the subset of n2_trains.
    """
    if not np_rnd:
        np_rnd = np.random 
    n2 = len(n2_trains)
    result=[[n2_trains[np_rnd.randint(0, n2),i] for i in xrange(time)] for j in xrange(n)]
    return np.array(result, dtype=int)

# -----------------------------------------------------------------------------
def generate_N_trains(N, N2, ISI, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains based on N2 uncorrelated possion 
    spike trains with the interspike interval "ISI" and time length "time". 
    "ISI" and "time" are in the units "dt".
    """
    import correlated_input as ci
    N2_trains = np.array([ci.poisson_train(ISI, time, np_rnd=np_rnd) for i in xrange(N2)])
    result = create_n_trains(N, N2_trains, time, np_rnd)
    return [np.array(result, dtype=int), np.array(N2_trains, dtype=int)]

# -----------------------------------------------------------------------------
def create_n_trains_sparse(n, n2_trains, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Makes n trains of spikes timings based on the subset of n2_trains in the 
    sparse form.
    """
    if not np_rnd:
        np_rnd = np.random 
    n2_limit = 1./len(n2_trains)
    result = [[] for _ in xrange(n)]
    spike_times = list(chain(*n2_trains))
    spike_times.sort()
    for spike_time in spike_times:
        copy_ind = np.arange(n)[np_rnd.random(n)<n2_limit]
        [result[i].append(spike_time) for i in copy_ind]
    np_result = [np.array(result[i]) for i in xrange(len(result))]
    return np_result

# -----------------------------------------------------------------------------
def generate_N_trains_sparse(N, N2, ISI, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains based on N2 uncorrelated possion 
    spike trains with the interspike interval "ISI" and time length "time". 
    "ISI" and "time" are in the units "dt".
    """
    import correlated_input as ci
    ci.poisson_train(ISI, time, np_rnd=np_rnd, dense=True) 
    N2_trains = [ci.poisson_train(ISI, time, np_rnd=np_rnd, dense=True) for _ in xrange(N2)]
    result = create_n_trains_sparse(N, N2_trains, time, np_rnd)
    return [result, N2_trains]

# -----------------------------------------------------------------------------
def corr2(p, n):
# -----------------------------------------------------------------------------
    """
    Returns the 2nd correlation value for the spike trains generated with
    average p and n independent cells (N2 in the paper of Destexhe).
    """
    return p*(1-p)/n
   
# -----------------------------------------------------------------------------
def corr3(p, n):
# -----------------------------------------------------------------------------
    """
    Returns the 3rd correlation value for the spike trains generated with
    average p and n independent cells (N2 in the paper of Destexhe).
    """
    return p*(1-p)*(1-2*p)/n**2

# -----------------------------------------------------------------------------
def corr4(p, n):
# -----------------------------------------------------------------------------
    """
    Returns the 4th correlation value for the spike trains generated with
    average p and n independent cells (N2 in the paper of Destexhe)
    """
    return ((p*(1-p)**3-p**2*(1-p)*(1-2*p))/n**3)+3*(n-1)*p**2*(1-p)**2/n**3
