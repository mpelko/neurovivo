"""
Creating correlated input trains
"""

import numpy as np
from pylab import *
import random
import cProfile

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
def binary_to_sparse_train(binary):
# -----------------------------------------------------------------------------
    """
    Takes the spike train in binary (001001001-like) format and returns:
    [sparse_train, length], where length is the length of the original train 
    and sparse_train is a list of spike timings.
    """
    sparse_train = []
    length = len(binary)
    for i,bit in enumerate(binary):
        if bit == 1:
            sparse_train.append(i)
    return [sparse_train, length]

# -----------------------------------------------------------------------------
def sparse_to_binary_train(sparse_train, length):
# -----------------------------------------------------------------------------
    """
    Takes the spike train in sparse format (5 12 2 4 -like) and it's length 
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
def create_n_trains(n, n2_trains, time):
# -----------------------------------------------------------------------------
    """
    Makes n trains of spikes timings based on the subset of n2_trains.
    """
    n2 = len(n2_trains)
    result = np.zeros([n,time])
    for i in xrange(time):
        for j in xrange(n):
            result[j][i] = n2_trains[random.randint(0, n2-1)][i]
    return result

# -----------------------------------------------------------------------------
def visualize_spike_train(trains):
# -----------------------------------------------------------------------------
    """
    Visulizes the spike trains. trains are given as a list of spike trains in 
    binary code.
    """
    #TODO: fix it!!!
    length = len(trains)
    y_width = 0.1
    figure
    for i, train in enumerate(trains):
        for spike in train:
            Line2D([spike, i-0.4],[spike, i+0.4], c = "blue")
    show()

# -----------------------------------------------------------------------------
def visualize_spike_train_fast(trains):
# -----------------------------------------------------------------------------
    """
    Visulizes the spike trains. trains are given as a list of spike trains in 
    binary code.
    """
    length = len(trains)
    y_width = 1
    for i, train in enumerate(trains):
        train_to_plot = np.array(binary_to_sparse_train(train)[0])
        y_position = np.zeros(len(train_to_plot)) + i * y_width
        plot(train_to_plot, y_position,".", color = "blue")
    xlabel("dt")
    ylabel("N (neuron spike trains)")
    ylim((-y_width,(i + 1) * y_width))

# -----------------------------------------------------------------------------
def plot_train_correlation(train1, train2, spike_width = 1):
# -----------------------------------------------------------------------------
    correlation = train_correlation(train1, train2, spike_width)
    d_length = (len(correlation)/2)
    plot(range(-d_length, d_length+1), correlation)
    xlabel("$\Delta$dt")
    ylabel("corr")
    
# -----------------------------------------------------------------------------
def correlation_sum(train1, train2, delta=0):
# -----------------------------------------------------------------------------
    """
    Calculates the correlation sum from two spike trains, where the 2nd is 
    shifted by delta.
    Spike trains should be np.arrays for efficiency.
    """
    # Making sure trains are Numpy arrays.
    np_train1 = np.array(train1)
    np_train2 = np.array(train2)
    result = 0
    if delta <= 0:
        result = sum(np_train1[:len(np_train1) + delta] * np_train2[-delta:])
    else:
        result = sum(np_train1[delta:] * np_train2[:len(np_train1) - delta])
    return result
   
# -----------------------------------------------------------------------------
def population_pairwise_correlation(trains):
# -----------------------------------------------------------------------------
    """
    Calculates the mean pairwise correlation and its std deviation in the
    population of spike trains (trains).
    """
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    n = len(np_trains)
    correlations = np.zeros((n,n))
    for i, train1 in enumerate(np_trains):
        for j, train2 in enumerate(np_trains[i+1:]):
            normalisation_sum = 1.*(sum(train1) + sum(train2))
            assert normalisation_sum != 0, "two trains withot spikes"
            correlations[i,j] = 2 * correlation_sum(train1, train2) /\
                                normalisation_sum
    correlations = correlations.flatten()
    correlations = correlations[correlations!=0]
    return [np.average(correlations), np.std(correlations)]

# -----------------------------------------------------------------------------
def train_correlation(train1, train2, spike_width = 1):
# -----------------------------------------------------------------------------
    """
    Returns the cross-correlation function of the two spike trains.
    The width of specific spike can be set (default = 1).
    """
    assert len(train1) == len(train2), "The lengths of the trains need to be"+\
                                       " the same to calc. the correlation."
    # TODO: Change the calculations of the wide spike train correlations!!!
    total_time = len(train1)
    corr_time = int(total_time / 1000)
    result = []
    delta = -corr_time
    np_train1 = np.array(train1)
    np_train2 = np.array(train2)
    
    if spike_width > 1:
        # Fixing the train1, so that it widens the spikes.
        for i,spike in enumerate(train1):
            if spike == 1:
                for j in xrange(spike_width-2):
                    if i+j+1 < total_time:
                        np_train1[i+j+1] = 1
        # Fixing the train2, so that it widens the spikes.
        for i,spike in enumerate(train2):
            if spike == 1:
                for j in xrange(spike_width-2):
                    if i+j+1 < total_time:
                        np_train2[i+j+1] = 1
                    
    for i in xrange((2 * corr_time) + 1):
        result.append(correlation_sum(np_train1, np_train2, delta))
        delta += 1
    result = 1.0*np.array(result) / sum(np_train1)
    return result
    
# -----------------------------------------------------------------------------
def common_elements(list1, list2):
# -----------------------------------------------------------------------------
    """
    Calculates the number of common elements of 2 lists.
    """
    result = 0;
    for element in list1:
        if list2.__contains__(element):
            result += 1
    return result

# -----------------------------------------------------------------------------
def correlated_trains_destexhe(N, N2, ISI, time):
# -----------------------------------------------------------------------------
    """
    Creates N correlated possion spike trains based on N2 uncorrelated possion 
    spike trains with the interspike interval "ISI" and time length "time". 
    "ISI" and "time" are in the units "dt".
    """
    N2_trains = [poisson_train(ISI, time) for i in xrange(N2)]
    result = create_n_trains(N, N2_trains, time)
    return result
   
# -----------------------------------------------------------------------------   
def perturb_spike_trains(trains, mean, sigma):
# -----------------------------------------------------------------------------    
    """
    Takes the spike trains in the binary format and randomly perturbs them by
    timing defined by mean and delta (gaussian distribution).
    """
    for i,train in enumerate(trains):
        print i
        trains[i] = perturb_spike_train(train, mean, sigma)
    return trains

# -----------------------------------------------------------------------------   
def perturb_spike_train(train, mean, sigma):
# -----------------------------------------------------------------------------    
    """
    Takes the spike train in the binary format and randomly perturbs it by
    timing defined by mean and delta (gaussian distribution).
    If the perturbed spike falls out of the length range, it is lost. 
    """
    new_train = [0 for i in train]
    for i, spike in enumerate(train):
        if spike == 1:
            candidate = i + int(random.gauss(mean, sigma) + 0.5)
            flag = True
            while flag and candidate < (len(train) - 1) and candidate >= 0:
                if new_train[candidate] != 1:
                    new_train[candidate] = 1
                    flag = False
                else:
                    candidate = i + int(random.gauss(mean, sigma) + 0.5)
    return new_train
    
# -----------------------------------------------------------------------------
def poisson_train(ISI, time):
# -----------------------------------------------------------------------------
    """
    Creates a spike train with the interspike interval "ISI" and time length 
    "time". "ISI" and "time" are in the units "dt".
    """
    train_length = 0
    while train_length < time:
        candidate = np.random.poisson(ISI, (time/ISI)*2)
        train_length = sum(candidate)
    i = 0; train_length = 0
    while train_length < time:
        train_length += candidate[i]
        i += 1
    assert sum(candidate[:(i-1)]) < time, "Spike train is longer then the total time."
    return make_array_train(candidate[:(i-1)], time)
    
# -----------------------------------------------------------------------------
def correlation_matrix(trains):
# -----------------------------------------------------------------------------
    """
    Creates a correlation matrix on the spike trains "trains".
    """
    result = [[0 for i in xrange(len(trains))] for j in xrange(len(trains))]
    for i in xrange(len(trains)):
        print i
        for j in xrange(len(trains)):
            if i>=j:
                corr = sum(np.array(trains[i])*np.array(trains[j]))*2.0/(sum(np.array(trains[i]))+sum(np.array(trains[j])))
                result[i][j] = corr
                result[j][i] = corr
    return result
    
# -----------------------------------------------------------------------------
def main():
# -----------------------------------------------------------------------------
    """
    Main function
    """
    time = 500000
    ISI = 100
    
    trains = correlated_trains_destexhe(10, 10, ISI, time)
    trains_perturbed = perturb_spike_trains(trains[:], 3, 5)
    trains_perturbed2 = perturb_spike_trains(trains[:], 0, 5)
    
    #visualize_spike_train_fast(trains)
    
    fig = figure(1)
    #ax = fig.add_subplot(212)
    plot_train_correlation(trains[7], trains[8], 20)
    plot_train_correlation(trains_perturbed[7], trains_perturbed[8],20)
    plot_train_correlation(trains_perturbed[7], trains_perturbed2[8],20)
    legend(["unperturbed","perturbed: mean = 0", "perturbed: mean = 3"])
    show()
    
    #imshow(correlation_matrix(trains), interpolation='nearest')
    #colorbar()
    
    binary_train = binary_to_sparse_train(trains[0])[0]
    ISIs = [binary_train[i+1] - binary_train[i] for i in xrange(len(binary_train)-1)]
    print min(ISIs), max(ISIs)
    hist(ISIs, 20)
    show()
    
    #plot(np.transpose([trains[7]), trains[8]]),np.transpose([np.zeros(len(train)), np.ones(len(train))]),"|")
    #ylim((-1,3))
    
if __name__ == "__main__":
    main()
    #cProfile.run("main()")
