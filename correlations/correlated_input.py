"""
Creating correlated input trains
"""

import numpy as np

DEBUG = False

# -----------------------------------------------------------------------------
def make_array_train(array, time):
# -----------------------------------------------------------------------------
    """
    Makes a train of spikes timings based on the array of inter-spike intervals
    """
    result = [0 for i in xrange(time)]
    
    spike = 0
    for number in array:
        spike += number
        if spike >= time:
            print sum(array), spike
        result[spike] = 1

    return result

# -----------------------------------------------------------------------------
def binary_to_sparse_train(binary, jitter = False, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Takes the spike train in binary (001001001-like) format and returns:
    [sparse_train, length], where length is the length of the original train 
    and sparse_train is a list of spike timings. If jitter = True the timings
    are randomly jittered within the time window.
    """
    if jitter and not np_rnd:
        np_rnd = np.random 
    
    length = len(binary)
    
    sparse_train = np.where(np.array(binary)==1)[0]
    
    if jitter:
        jitt = np_rnd.rand(len(sparse_train))
        sparse_train = sparse_train + jitt

    return [sparse_train, length]

# -----------------------------------------------------------------------------
def sparse_to_binary_train(sparse_train, length):
# -----------------------------------------------------------------------------
    """
    Takes the spike train in sparse format - ISI timings (5 12 2 4 -like) and it's length 
    where length is the length of the original train and sparse_train is a list 
    of spike timings. format and returns a binary representation (01010100010).
    """
    assert sum(sparse_train) < length, "Cannot create a binary_train. The" +\
           " specified length is too small."
    
    binary_train = [0 for i in length]
    spike = 0
    for ISI in sparse_train:
        spike += ISI
        binary_train[spike] = 1
    return binary_train
    
# -----------------------------------------------------------------------------
def create_n_trains(n, n2_trains, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Makes n trains of spikes timings based on the subset of n2_trains.
    """
    if not np_rnd:
        np_rnd = np.random 

    n2 = len(n2_trains)
    result = np.zeros([n,time])
    for i in xrange(time):
        for j in xrange(n):
            result[j][i] = n2_trains[np_rnd.randint(0, n2)][i]
    return result

# -----------------------------------------------------------------------------
def correlated_trains_ising_simple(N, mean, corr, time=1000, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains following the minimum model
    (maximum-entropy distribution), which can be written in the Ising form as:
    ... TODO ...
    """
    
    import correlations.simple_ising as si 
    [h, J] = si.find_good_hJ_parameters(N, mean, corr)
    result = si.generate_N_spike_trains(N, h, J, time, np_rnd)
    
    return result

# -----------------------------------------------------------------------------
def correlated_trains_destexhe(N, N2, ISI, time=1000, np_rnd=None, sparse=False):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains .
    """
    import neurovivo.correlations.destexhe as des
    if sparse:
        result = des.generate_N_trains_sparse(N,N2,ISI,time,np_rnd)
    else:
        result = des.generate_N_trains(N,N2,ISI,time,np_rnd)
    return result
   
# -----------------------------------------------------------------------------   
def perturb_spike_trains(trains, mean, sigma, np_rnd=None):
# -----------------------------------------------------------------------------    
    """
    Takes the spike trains in the binary format and randomly perturbs them by
    timing defined by mean and delta (gaussian distribution).
    """
    perturbed_trains = np.array(trains[:])
    for i,train in enumerate(trains):
        #print i
        perturbed_trains[i] = perturb_spike_train(train, mean, sigma, np_rnd)
    return perturbed_trains

# -----------------------------------------------------------------------------   
def perturb_spike_train(train, mean, sigma, np_rnd):
# -----------------------------------------------------------------------------    
    """
    Takes the spike train in the binary format and randomly perturbs it by
    timing defined by mean and delta (gaussian distribution).
    If the perturbed spike falls out of the length range, it is lost. 
    """

    if not np_rnd:
        np_rnd = np.random

    new_train = [0 for i in train]
    for i, spike in enumerate(train):
        if spike == 1:
            candidate = i + int(np_rnd.normal(mean, sigma) + 0.5)
            flag = True
            while flag and candidate < (len(train) - 1) and candidate >= 0:
                if new_train[candidate] != 1:
                    new_train[candidate] = 1
                    flag = False
                else:
                    candidate = i + int(np_rnd.normal(mean, sigma) + 0.5)
    return new_train
    
# -----------------------------------------------------------------------------
def poisson_train_old(ISI, time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates a spike train with the interspike interval "ISI" and time length 
    "time". "ISI" and "time" are in the units "dt".
    """
    if not np_rnd:
        np_rnd = np.random

    train_length = 0
    while train_length < time:
        candidate = np_rnd.poisson(ISI, (time/ISI)*2)
        train_length = sum(candidate)
    i = 0; train_length = 0
    while train_length < time:
        train_length += candidate[i]
        i += 1
    assert sum(candidate[:(i-1)]) < time, "Spike train is longer then the total time."
    return make_array_train(candidate[:(i-1)], time)

# -----------------------------------------------------------------------------
def poisson_train(ISI, time, bin_size=1, np_rnd=None, dense=False):
# -----------------------------------------------------------------------------
    """
    Creates a spike train with the interspike interval "ISI" and time length 
    "time". "ISI" and "time" are in the units "dt".
    """
    if not np_rnd:
        np_rnd = np.random

    # choose how many spikes there should be
    num_spike = np_rnd.poisson(1.*time/ISI)
    if not dense:
        assert num_spike < (time/bin_size), "The Possion train can not be so dense."
    
    # create those spikes in the continuous space
    continuous_spikes = np_rnd.rand(num_spike)*time
    continuous_spikes.sort()
    
    if dense:
        return continuous_spikes
    # discretise the spikes in 1ms time bin
    discr_spikes = discretize_train(continuous_spikes, time, 1, True, np_rnd)
    
    return discr_spikes

# -----------------------------------------------------------------------------
def discretize_train(spike_times, time, bin_size, redistribute=True, 
                     np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Takes the spike times given as real numbers and returns a spike train in
    bins (1 there is a spike, 0 there is no spike).
    If redistribute = True, the extra spikes (more then 1 in a bin) are 
    randomly redistributed the available bins.   
    """
    #print spike_times
    nr_bins = int(np.round(time/bin_size))
    if redistribute:
        assert nr_bins >= len(spike_times), "too many spikes or too few bins."
    
    if not np_rnd:
        np_rnd = np.random

    if len(spike_times) == 0:
        return np.zeros(nr_bins, dtype=int)

    # discretise the spikes
    discr_spikes = np.histogram(spike_times, 1.0*time*np.arange(0, nr_bins+1)/(nr_bins+1))[0]
    #print discr_spikes
    # cut out the extra spikes
    discr_spikes = np.array(discr_spikes >= 1, dtype=int)
    if redistribute:
        
        # how many spikes need to be redistributed (as there are more then 1 in 
        # the bin)
        redistribute_nr = len(spike_times) - sum(discr_spikes)
        #print redistribute_nr
        while redistribute_nr > 0:
            candidate = np_rnd.randint(time)
            if discr_spikes[candidate] == 0:
                discr_spikes[candidate] = 1
                redistribute_nr -= 1
    
    return discr_spikes

# -----------------------------------------------------------------------------
def alternative_poisson(rate, time, dt=0.0001, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates the posson train by advancing in time by dt and setting the spike
    if the gods of randomness decide to do so. The rate parameter should give a
    proper rate of the train (typically a bit less), the time gives the length 
    of the spike train in seconds.
    In order for this to work, rate*dt should be smaller then THRESHOLD 
    (= probability that we would have 2 spikes in that time window).
    the output are the times when the train spikes (in seconds). 
    """
    THRESHOLD = 0.01
    assert rate*dt < THRESHOLD, "use a smaller rate or smaller dt"
    
    if not np_rnd:
        np_rnd = np.random
    randoms = np_rnd.rand(int(time/dt))
    spikes = []
    limit = rate*np.exp(-rate*dt/2.)*dt
    for i,rand in enumerate(randoms):
        if rand < limit:
            spikes.append(i*dt)
    return spikes

# -----------------------------------------------------------------------------
def modulated_poisson_old(R0, amp, freq, time, dt=0.0001, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates the posson train by advancing in time by dt and setting the spike
    if the gods of randomness decide to do so. 
    The rate of the train is dynamic and changes with time (t) as:
    rate = R0*(1*amp*sin(freq*t))  
    
    The time gives the length of the spike train in seconds.
    In order for this to work, (R0*(1+amp))*dt should be smaller then THRESHOLD 
    (= probability that we would have 2 spikes in that time window).
    the output are the times when the train spikes (in seconds). 
    """
    THRESHOLD = 0.01
    assert (R0*(1+amp))*dt < THRESHOLD, "use a smaller mean rate or smaller dt"
    
    if not np_rnd:
        np_rnd = np.random
    
    spikes = np.arange(0,time,dt)
    randoms = np_rnd.rand(len(spikes))
    limits = R0*(1+amp*np.sin(freq*spikes))*np.exp(-R0*dt/2.)*dt
    return spikes[randoms<limits] 

# -----------------------------------------------------------------------------
def modulated_poisson(modulation_function, spike_tries, total_time, np_rnd=None):
# -----------------------------------------------------------------------------
    """
    Creates the modulated Poisson train.
    
    modulation_function - Probability function of time (argument is 
                          time in ms).
    spike_tries         - The numer of tries for a spike. Together with 
                          modulation function this determines the mean rate.
    total_time          - total length of the spike train (in ms)
    
    The output are the times when the train spikes (in ms). 
    """
    if not np_rnd:
        np_rnd = np.random
    
    rand_tries = np_rnd.rand(spike_tries,2)
        
    try_values = modulation_function(rand_tries[:,0]*total_time)
    
    good_idx = try_values > rand_tries[:,1]
    
    spike_times = rand_tries[good_idx,0]*total_time

    spike_times.sort()
            
    return spike_times

# -----------------------------------------------------------------------------
def get_sin_mod_parameters(freq, total_time, mean_rate, modulation_strength):
# -----------------------------------------------------------------------------
    """
    Returns the paramateres (modulation_function, spike_tries) for the 
    sinusoidal modulation of a certain modulation strength for the creation of
    rate-modulated spike trains.
    freq -                the frequency of modulation in Hz
    total_time -          the total time of spike train in ms
    mean_rate -           in Hz
    modulation_strength - a number in the interval [0,1], where 1 means 
                          complete modulation and 0 no modulation
    
    The output is in a form (modulation_function, spike_tries)
    """
    import scipy.integrate
    
    def mod_function(time):
        return (1-modulation_strength/2.) + modulation_strength*np.sin(np.pi*2*freq*time/1000)/2.

    prob_lost_spikes = 1. - (1.*scipy.integrate.quad(mod_function,0,total_time)[0]/total_time)
    
    spike_tries = int(round((mean_rate * total_time/1000.)/(1-prob_lost_spikes)))
    
    return (mod_function, spike_tries) 

# -----------------------------------------------------------------------------
def get_sin_mod_OU_parameters(freq, theta, sigma, total_time, mean_rate, modulation_strength, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Returns the paramateres (modulation_function, spike_tries) for the 
    sinusoidal modulation with the OU noise in the phase noise 
    of a certain modulation strength for the creation of
    rate-modulated spike trains.
    freq -                the frequency of modulation in Hz
    theta -               the OU noise parameter
    sigma -               the phase noise parameter - how much the phase oscillates around 0
    total_time -          the total time of spike train in ms
    mean_rate -           in Hz
    modulation_strength - a number in the interval [0,1], where 1 means 
                          complete modulation and 0 no modulation
    
    The output is in a form (modulation_function, spike_tries)
    """
    import scipy.integrate
    from scipy.interpolate.interpolate import interp1d
    import neurovivo.common as cmn

    times = precision*np.arange(int(total_time/precision))
    #theta = 1.*theta/precision
    OU_noise = cmn.OU_process_basic(theta, times)
    OU_interpolator = interp1d(times, OU_noise, "linear", bounds_error=False)    
    def mod_function(time):
        return (1-modulation_strength/2.) + modulation_strength*np.sin(np.pi*2*freq*time/1000+2*sigma*np.pi*OU_interpolator(time))/2.

#    discrete_function = [mod_function(time) for time in times]
    discrete_function = mod_function(times)

    prob_lost_spikes = 1. - (1.*scipy.integrate.trapz(discrete_function, dx=precision)/total_time)
    
    spike_tries = int(round((mean_rate * total_time/1000.)/(1-prob_lost_spikes)))
    
    return (mod_function, spike_tries) 

# -----------------------------------------------------------------------------
def get_basic_OU_function(theta=0.001, total_time=1000, precision=0.1):
# -----------------------------------------------------------------------------
    from scipy.interpolate.interpolate import interp1d
    import neurovivo.common as cmn

    times = precision*np.arange(int(total_time/precision))
    #theta = 1.*theta/precision
    OU_noise = cmn.OU_process_basic(theta, times)
    
    OU_interpolator = interp1d(times, OU_noise, "linear", bounds_error=False)
    
    return OU_interpolator

# -----------------------------------------------------------------------------
def get_sin_mod_OU_function(freq, theta, sigma, total_time, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Returns the modulation_function for the sinusoidal modulation with the 
    OU noise in the phase noise for the creation of rate-modulated spike trains.
    freq -                the frequency of modulation in Hz
    theta -               the OU noise parameter
    sigma -               the phase noise parameter - how much the phase oscillates around 0
    total_time -          the total time of spike train in ms
    mean_rate -           in Hz
    
    The output is in a form (modulation_function, spike_tries)
    """
    from scipy.interpolate.interpolate import interp1d
    import neurovivo.common as cmn

    times = precision*np.arange(int(total_time/precision))
    OU_noise = cmn.OU_process_basic(theta, times)
    OU_interpolator = interp1d(times, OU_noise, "linear", bounds_error=False)    
    def mod_function(time):
        return np.sin(np.pi*2*(freq*time/1000+sigma*OU_interpolator(time)))
    
    return mod_function

# -----------------------------------------------------------------------------
def get_OU_and_sin_modulated_function(OU_func, sinOU_func, OU_factor=1, sin_factor=1, total_time=1000, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Returns the modulation_function for the OU proces including the
    sinusoidal modulation (OU modulated in phase) of a certain modulation 
    strength for the creation of rate-modulated spike trains.
    OU_func  -            the OU function
    sinOU_func -          the OU modulated sine function
    OU_factor -           factor for OU function
    sin_factor -          factor for OU function
    total_time -          the total time of spike train in ms
    precision  -          the discretization parameter in ms
    
    The output is in a modulation_function of time
    """

    def tmp_function(time):
        return 2 + OU_factor*OU_func(time) + sin_factor*sinOU_func(time)
    
    times = precision*np.arange(int(total_time/precision))
    #sample_max_inv = 1./np.max([tmp_function(time) for time in times])
    sample_max_inv = 1./np.max(tmp_function(times))
    
    # making sure function returns a value between 0 and 1
    def mod_function(time):
        tmp = sample_max_inv*tmp_function(time)
        tmp[tmp<0]=0
        return tmp
    
    return mod_function 


# -----------------------------------------------------------------------------
def get_OU_and_sin_modulated_function_medium_old(OU_func, OU_factor=1, sin_freq=2.5, sin_factor=1, total_time=1000, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Returns the modulation_function for the OU proces including the
    sinusoidal modulation of a certain modulation strength for the creation of
    rate-modulated spike trains.
    theta -               the OU noise parameter
    OU_factor -           the OU process multiplier
    sin_freq -            the frequency of modulation in Hz
    sin_factor -          the sin process multiplier
    total_time -          the total time of spike train in ms
    precision  -          the discretisation parameter in ms
    
    The output is in a modulation_function of time
    """
    import neurovivo.common as cmn
        
    def tmp_function(time):
        return 1 + OU_factor*OU_func(time) + 1 + sin_factor*np.sin(np.pi*2*sin_freq*time/1000)
    
    times = precision*np.arange(int(total_time/precision))
    sample_max_inv = 1./np.max([tmp_function(time) for time in times])
    
    # making sure function returns a value between 0 and 1
    def mod_function(time):
        tmp = sample_max_inv*tmp_function(time)
        if tmp<0:
            tmp = 0
        return tmp
    
    return mod_function 

# -----------------------------------------------------------------------------
def get_OU_and_sin_modulated_function_old(theta=0.001, OU_factor=1, sin_freq=2.5, sin_factor=1, total_time=1000, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Returns the modulation_function for the OU proces including the
    sinusoidal modulation of a certain modulation strength for the creation of
    rate-modulated spike trains.
    theta -               the OU noise parameter
    OU_factor -           the OU process multiplier
    sin_freq -            the frequency of modulation in Hz
    sin_factor -          the sin process multiplier
    total_time -          the total time of spike train in ms
    precision  -          the discretisation parameter in ms
    
    The output is in a modulation_function of time
    """
    from scipy.interpolate.interpolate import interp1d
    import neurovivo.common as cmn

    times = precision*np.arange(int(total_time/precision))
    #theta = 1.*theta/precision
    OU_noise = cmn.OU_process_basic(theta, times)
    
    OU_interpolator = interp1d(times, OU_noise, "linear", bounds_error=False)    
        
    def tmp_function(time):
        return 1 + OU_factor*OU_interpolator(time) + 1 + sin_factor*np.sin(np.pi*2*sin_freq*time/1000)
    
    
    sample_max_inv = 1./np.max([tmp_function(time) for time in times])
    
    # making sure function returns a value between 0 and 1
    def mod_function(time):
        tmp = sample_max_inv*tmp_function(time)
        if tmp<0:
            tmp = 0
        return tmp
    
    return mod_function 

# -----------------------------------------------------------------------------
def get_number_of_tries(func, mean_rate=10, total_time=1000, precision=0.1):
# -----------------------------------------------------------------------------
    """
    Given a dynamic rate function (func), total time (total_time) and the 
    desired mean rate (mean_rate) it calculates the required number of tries
    for spike train in order to produce a spike train with the the desired mean
    rate.
    func -                the function (must be defined for the whole time 
                          period and should not be bigger then 1 in this period)
    rate -                desired mean rate
    total_time -          the total time of spike train in ms
    mean_rate -           in Hz
    precision -           a discretization parameter in ms.
    
    The output is in a form (modulation_function, spike_tries)
    """
    import scipy.integrate

    times = precision*np.arange(int(total_time/precision))
#    samples = [func(time) for time in times]
    #samples = func(times)
    #assert np.max(samples) <=1 and np.min(samples)>=0, "{}, {}".format(np.max(samples), np.min(samples))

#    discrete_function = [func(time) for time in times]
    discrete_function = func(times)
    
    prob_lost_spikes = 1. - (1.*scipy.integrate.trapz(discrete_function, dx=precision)/total_time)
    
    spike_tries = int(round((mean_rate * total_time/1000.)/(1-prob_lost_spikes)))
    
    return spike_tries