"""
Analysis of the correlations.
"""

import numpy as np
import correlations.correlated_input as ci

# -----------------------------------------------------------------------------
def correlation_order4(x, y, z, q):
# -----------------------------------------------------------------------------
    """
    Calculates the pearson-like correlation between the vectors x,y,z,q.
    Difference to the order2 case is that here we don't normalize.
    x and y have to be numpy-array-supported types.
    """
    ### make the functions 0-mean and 1-std.###
    std_x = np.std(x); std_y = np.std(y); std_z = np.std(z); std_q = np.std(q)
        # for constant functions, the coefficint is 0 
    if std_x == 0 or std_y == 0 or std_z == 0 or std_q == 0:
        return 0
    x = x - np.mean(x)
    y = y - np.mean(y)
    z = z - np.mean(z)
    q = q - np.mean(q)
   
    ### calculate and return the result ###
    return sum(x*y*z*q)/len(x)

# -----------------------------------------------------------------------------
def correlation_order3(x, y, z):
# -----------------------------------------------------------------------------
    """
    Calculates the pearson-like correlation between the vectors x,y,z.
    Difference to the order2 case is that here we don't normalize.
    x and y have to be numpy-array-supported types.
    """
    ### make the functions 0-mean and 1-std.###
    std_x = np.std(x); std_y = np.std(y); std_z = np.std(z)
        # for constant functions, the coefficint is 0 
    if std_x == 0 or std_y == 0 or std_z == 0:
        return 0
    x = x - np.mean(x)
    y = y - np.mean(y)
    z = z - np.mean(z)
   
    ### calculate and return the result ###
    return np.sum(x*y*z)/len(x) 

# -----------------------------------------------------------------------------
def correlation_order2(x, y):
# -----------------------------------------------------------------------------
    """
    Calculates the pearson-like correlation between the vectors x,y.
    Difference to the order2 case is that here we don't normalize.
    x and y have to be numpy-array-supported types.
    """
    ### make the functions 0-mean and 1-std.###
    std_x = np.std(x); std_y = np.std(y)
        # for constant functions, the coefficint is 0 
    if std_x == 0 or std_y == 0:
        return 0
    x = x - np.mean(x)
    y = y - np.mean(y)
   
    ### calculate and return the result ###
    return np.sum(x*y)/len(x) 

# -----------------------------------------------------------------------------
def pearson_correlation(x, y):
# -----------------------------------------------------------------------------
    """
    Calculates the pearson correlation between the vectors train1, train2.
    http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
    x and y have to be numpy-array-supported types.
    """
    #import scipy.stats.stats as stats
    #return stats.pearsonr(x,y)[0]
    
    ### make the functions 0-mean and 1-std. ###
    std_x = np.std(x); std_y = np.std(y)
        # for constant functions, the coefficint is 0 
    if std_x == 0 or std_y == 0:
        return 0
    x = (x - np.mean(x))/std_x
    y = (y - np.mean(y))/std_y

    ### calculate and return the result ###
    return sum(x*y)/len(x) 
    
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
   
    # cutting out the overlaping parts
    if delta <= 0:
        train1_part = np_train1[:len(np_train1) + delta]
        train2_part = np_train2[-delta:]
    else:
        train1_part = np_train1[delta:]
        train2_part = np_train2[:len(np_train1) - delta]
   
    # setting mean to 0 and variation to 1
    result = pearson_correlation(train1_part, train2_part)
    
    return result
   
# -----------------------------------------------------------------------------
def population_pairwise_pearson_correlation(trains):
# -----------------------------------------------------------------------------
    """
    Calculates the mean pairwise correlation (pearson) and its std deviation in
    the population of spike trains (trains).
    """
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    n = len(np_trains)
    correlations = np.zeros(n*(n-1)/2)
    corr_calc = 0
    for i, train1 in enumerate(np_trains):
        for train2 in np_trains[i+1:]:
            correlations[corr_calc] = pearson_correlation(train1, train2)
            corr_calc += 1
    return [np.average(correlations), np.std(correlations)]

# -----------------------------------------------------------------------------
def population_correlation_order2(trains):
# -----------------------------------------------------------------------------
    """
    Calculates the mean correlation of order 2 and its std deviation in the
    population of spike trains (trains).
    """
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    n = len(np_trains)
    correlations = np.zeros(n*(n-1)/2)
    corr_calc = 0
    for i, train1 in enumerate(np_trains):
        for train2 in np_trains[i+1:]:
            correlations[corr_calc] = correlation_order2(train1, train2)
            corr_calc += 1
    return [np.average(correlations), np.std(correlations)]

# -----------------------------------------------------------------------------
def population_correlation_order3(trains):
# -----------------------------------------------------------------------------
    """
    Calculates the mean correlation of order 3 and its std deviation in the
    population of spike trains (trains).
    """
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    n = len(np_trains)
    correlations = np.zeros(n*(n-1)*(n-2)/6)
    corr_calc = 0
    for i, train1 in enumerate(np_trains):
        for j, train2 in enumerate(np_trains[i+1:]):
            for k, train3 in enumerate(np_trains[i+j+2:]):
                correlations[corr_calc] =\
                                    correlation_order3(train1, train2, train3)
                corr_calc += 1
    return [np.average(correlations), np.std(correlations)]

# -----------------------------------------------------------------------------
def population_correlation_order4(trains):
# -----------------------------------------------------------------------------
    """
    Calculates the mean correlation of order 4 and its std deviation in the
    population of spike trains (trains).
    """
    assert len(trains) > 3;
    
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    n = len(np_trains)
    correlations = np.zeros(n*(n-1)*(n-2)*(n-3)/24)
    corr_calc = 0
    for i, train1 in enumerate(np_trains):
        for j, train2 in enumerate(np_trains[i+1:]):
            for k, train3 in enumerate(np_trains[i+j+2:]):
                for train4 in np_trains[i+j+k+3:]:
                    correlations[corr_calc] =\
                                    correlation_order4(train1, train2, train3, train4)
                corr_calc += 1
    return [np.average(correlations), np.std(correlations)]

# -----------------------------------------------------------------------------
def population_correlation_estimate(trains, order = 2, combinations = 10000):
# -----------------------------------------------------------------------------
    """
    Estimates mean correlation of specified order and its std deviation in the
    population of spike trains (trains).
    It does so by calculating from a specified number of combinations.
    """
    
    n = len(trains)
    time = len(trains[0])
    assert n >= order,\
            "Too few spike trains for order {0}".format(order)
    assert order >= 2, "Only correlations of order 2 or higher are sensible"
    
    # Making sure trains are Numpy arrays.
    np_trains = np.array(trains)
    #max_combinations = gmpy.bincoef(n,order)
    #combinations = min(max_combinations, combinations)
    correlations = np.zeros(combinations)
    for i in xrange(combinations):
        np.random.shuffle(np_trains)       # to get new batch of trains
        trains_to_multiply = np_trains[:order,:]
        res = trains_to_multiply[0][:]-np.mean(trains_to_multiply[0])
        #print "****"
        #print trains_to_multiply
        #print "########"
        for j in range(1,order):
            res *= (trains_to_multiply[j]-np.mean(trains_to_multiply[j]))
            #print res
        #print 1.0*sum(res)/time
        correlations[i] = 1.0*sum(res)/time

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
    result = 1.0*np.array(result) 
    return result
    
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
                norm = (sum(np.array(trains[i]))+sum(np.array(trains[j])))
                assert norm > 0, "at least two spike trains are empty!"
                corr = sum(np.array(trains[i])*np.array(trains[j]))*2.0/norm
                result[i][j] = corr
                result[j][i] = corr
    return result

# -----------------------------------------------------------------------------
def spike_probability_histogram(trains):
# -----------------------------------------------------------------------------
    """
    Only here due to backward compatibility in the code. To disapear once I 
    refractore it. The name was misleading.
    TODO: refractore
    """
    
    return spike_count_histogram(trains)

# -----------------------------------------------------------------------------
def spike_count_histogram(trains):
# -----------------------------------------------------------------------------
    """
    Given a spike train it spits out the counts of occurances of patterens with
    a number of spikes in patteren.
    """
    N = len(trains)

    spikes_per_pateren = np.sum(trains, axis=0)
    counts = np.zeros(N+1)
    for i in xrange(N+1):
        counts[i] = len(spikes_per_pateren[spikes_per_pateren[:] == i]) 
    
    return counts

# -----------------------------------------------------------------------------
def full_visual_analysis(trains, title=None):
# -----------------------------------------------------------------------------
    """
    Visually plots teh spike analysis.
    """
    # TODO - for a spike train, plot the train, plot the P_prob, plot the 
    # averages per each train, calculate higher order distributions
    import matplotlib.pylab as pyl
    from plotting import plotting
    
    pyl.figure(2351, figsize=(11,9))
    if title:
        pyl.suptitle("Viusal train analysis: {0}". title)
    else:
        pyl.suptitle("Viusal train analysis")
    sp1 = pyl.subplot(2,1,1)
    plotting.plot_spike_trains(trains)
    sp2 = pyl.subplot(2,2,3)
    P_prob_num = spike_probability_histogram(trains)
    pyl.bar(range(1,len(P_prob_num)), P_prob_num[1:], align="center")
    sp2.set_xlabel("spike # per time-bin")
    sp3 = pyl.subplot(2,2,4)
    means = np.mean(trains, 1)
    if len(trains) < 100:
        pyl.bar(range(len(trains)), means, align="center")
        sp3.set_ylabel("mean")
        sp3.set_xlabel("trains")
    else:
        pyl.hist(means)
        sp3.set_ylabel("# of trains")
        sp3.set_xlabel("mean")
    
    
    pyl.figure(251, figsize=(10,8))
    sp1 = pyl.subplot(1,1,1)
    limit_order = min(10,len(trains))
    [corrs, corr_err] = correlations_hist(trains, limit_order=limit_order)
    sp1.bar(range(2,limit_order+1), corrs, yerr=corr_err, align="center")
    

# -----------------------------------------------------------------------------
def correlations_hist(trains, combinations = 50, limit_order = None):
# -----------------------------------------------------------------------------
    """
    Returns the correlations (orders up to limit_order). If limit_order is None, 
    then it plots all higher order correlations (can be lenghty with many 
    neurons).
    """
    
    if not limit_order:
        limit_order = len(trains)
    corrs = np.zeros(limit_order-1)
    corr_err = np.zeros(limit_order-1)
    for order in xrange(2,limit_order+1):
        [corrs[order-2], corr_err[order-2]] =\
                   population_correlation_estimate(trains, order, combinations)
    
    return [corrs, corr_err]
    
# -----------------------------------------------------------------------------
def mean_from_P_probs(N, P_probs):
# -----------------------------------------------------------------------------
    """
    Average of the spike train - probability that a spike is elicited in the
    time-bin (this related to spiking rate via the time-bin width).
    N - number of neurons
    P_probs - array of len(P_probs) < N elements, where the i-th element gives 
    the probability of i spikes being elicited in a population at a certain 
    time-bin
    """
    __verify_parameters__(N, P_probs)
    
    return sum([(1.0*P/N)*P_probs[P] for P in xrange(len(P_probs))])

# -----------------------------------------------------------------------------
def corr2_from_P_probs(N, P_probs):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of the 2nd order correlation.
    N - number of neurons
    P_probs - array of len(P_probs) < N elements, where the i-th element gives 
    the probability of i spikes being elicited in a population at a certain 
    time-bin

    """
    __verify_parameters__(N, P_probs)
    
    avg_x = mean_from_P_probs(N, P_probs)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_probs[P]\
                  for P in xrange(2,len(P_probs))])
    result = (avg_xy - avg_x**2)
    return result 

# -----------------------------------------------------------------------------
def corr3_from_P_probs(N, P_probs):
# -----------------------------------------------------------------------------
    """
    Theoretical calculation of correlation of 3rd order.
    N - number of neurons
    P_probs - array of len(P_probs) < N elements, where the i-th element gives 
    the probability of i spikes being elicited in a population at a certain 
    time-bin
    """
    __verify_parameters__(N, P_probs)

    avg_x = mean_from_P_probs(N, P_probs)
    avg_xy = sum([1.0*P*(P-1)/(1.0*N*(N-1))*P_probs[P]\
                  for P in xrange(2,len(P_probs))])
    avg_xyz = sum([1.0*P*(P-1)*(P-2)/(1.0*N*(N-1)*(N-2))*P_probs[P]\
                   for P in xrange(3,len(P_probs))])
    result = avg_xyz - 3 * avg_xy*avg_x + 2 * avg_x**3
    
    return result 

# -----------------------------------------------------------------------------
def __verify_parameters__(N, P_probs):
# -----------------------------------------------------------------------------
    """
    Checks the parameters are valid.
    N - number of neurons
    P_probs - array of len(P_probs) < N elements, where the i-th element gives 
    the probability of i spikes being elicited in a population at a certain 
    time-bin
    """
    assert N > 0
    assert abs(np.sum(P_probs)-1) < 0.00001, "sum of probabilities should be"+\
                                    " 1 and not {0}".format(np.sum(P_probs))
    
# -----------------------------------------------------------------------------
def find_P_prob_params(N, mean, corr2):
# -----------------------------------------------------------------------------
    """
    Taking the parameters it returns P_probs - the probability distribution for
    all Ps, including P=0 
    """
    #from scipy.optimize import fmin_l_bfgs_b
    from scipy.optimize import fmin_tnc # seems better

    assert N>2
    assert mean < 1
    # TODO - assert corr2 < f(N,mean) 

    mean_des = mean           # desired average
    corr_des = corr2          # desired correlation
    P_probs_m = np.zeros(N)   # candidate P_probs_m
    P_probs_m[-1] = 1.        # candidate parameter J
    
    if mean == 0:
        result = np.zeros(N+1)
        result[0] = 1.
        return result
    
    def opt_func(P_probs_m, N, mean_des, corr_des):
        P_probs = __get_P_probs__(N, mean_des, P_probs_m)
        corr_curr = corr2_from_P_probs(N, P_probs)
        corr_diff = (corr_curr - corr_des)
        #print P_probs_m, P_probs, corr_curr 
        return abs(corr_diff)/abs(corr_des + 0.0000001)
    
    bounds = [(0,1) for i in xrange(N)]
    opt_result = fmin_tnc(opt_func, P_probs_m, args=(N,mean_des,corr_des),\
                              approx_grad = True, bounds=bounds)
    
    return __get_P_probs__(N, mean_des, opt_result[0])

# -----------------------------------------------------------------------------
def find_P_prob_params_corr3(N, mean, corr2, corr3, method="tnc"):
# -----------------------------------------------------------------------------
    """
    Taking the parameters it returns P_probs - the probability distribution for
    all Ps, including P=0 
    """
    #from scipy.optimize import fmin_l_bfgs_b
    from scipy.optimize import fmin_tnc # seems better
    from scipy.optimize import anneal
    from scipy.optimize import brute

    assert N > 2
    assert mean < 1
    # TODO - assert corr2 < f(N,mean)
    # TODO - assert corr3 < f(N,mean,corr2) 

    mean_des = mean           # desired average
    corr_des = corr2          # desired correlation
    corr3_des = corr3          # desired correlation
    P_probs_m = np.zeros(N)   # candidate P_probs_m
    P_probs_m[-1] = 1.        # candidate parameter J
    
    if mean == 0:
        result = np.zeros(N+1)
        result[0] = 1.
        return result
    
    def opt_func(P_probs_m, N, mean_des, corr_des):
        P_probs = __get_P_probs__(N, mean_des, P_probs_m)
        corr_curr = corr2_from_P_probs(N, P_probs)
        corr_diff = (corr_curr - corr_des)
        corr3_curr = corr3_from_P_probs(N, P_probs)
        corr3_diff = (corr3_curr - corr3_des)
        #print corr_curr, corr3_curr, P_probs_m
        #return (abs(corr_diff)/abs(corr_des + 0.0000001) + abs(corr3_diff)/abs(corr3_des + 0.0000001))**.1
        return 2*abs(corr_diff) + abs(corr3_diff)
    
    if method == "tnc":
        bounds = [(0,1) for i in xrange(N)]
        opt_output = fmin_tnc(opt_func, P_probs_m, args=(N,mean_des,corr_des),\
                            approx_grad = True, bounds=bounds)
        print opt_output
        opt_result = opt_output[0]
        
    if method == "anneal":
        opt_output = anneal(opt_func, P_probs_m, args=(N,mean_des,corr_des),\
                        lower = np.zeros(N), upper = np.ones(N))
        print opt_output
        opt_result = opt_output[0]
    
    if method == "brute":
        opt_result = brute(opt_func, args=(N,mean_des,corr_des),\
                        ranges = [(0.00000000001,0.99999999999) for i in xrange(N)], Ns=20)
    
    return __get_P_probs__(N, mean_des, opt_result)

#------------------------------------------------------------------------------ 
def __get_P_probs__(N, mean, P_probs_m):
#------------------------------------------------------------------------------ 
    P_probs_mm = 1.0*np.array(P_probs_m)/sum(P_probs_m)
    mean_P = np.sum([P_probs_mm[i]*(i+1) for i in xrange(len(P_probs_mm))])
    if mean_P == 0:
        return np.concatenate([[1],P_probs_mm])
    P_probs = np.array(P_probs_mm) * mean * N / mean_P
    return np.concatenate([[1-sum(P_probs)],P_probs]) 

#------------------------------------------------------------------------------ 
def __get_P_probs_m__(P_probs):
#------------------------------------------------------------------------------ 
    P_probss = np.array(P_probs)
    return P_probss[1:]/sum(P_probss[1:]) 

#------------------------------------------------------------------------------ 
def print_info_from_P_probs(N, P_probs):
#------------------------------------------------------------------------------ 
    print "mean = {0}".format(mean_from_P_probs(N, P_probs))
    print "corr2 = {0}".format(corr2_from_P_probs(N, P_probs))
    if N>=3:
        print "corr3 = {0}".format(corr3_from_P_probs(N, P_probs))
       
#------------------------------------------------------------------------------ 
def ISI_histogram(trains, bins=10):
#------------------------------------------------------------------------------ 
    """
    Returns the hsitogram of the ISI intervals along with the std deviation and
    border values.
    The bins argument follows the rules of the bins argument in 
    numpy.histogram: 
    http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
    """
    import itertools
    
    if isinstance(bins, int):
        nr_bins = bins
    else:
        nr_bins = len(bins)-1
    N = len(trains)
        
    ISI_hists = np.zeros((N, nr_bins))
    trains_ISIs = []
    
    for train in trains:
        #print (train)
        sparse_train = np.array(ci.binary_to_sparse_train(train)[0])
        trains_ISIs.append(sparse_train[1:] - sparse_train[:-1])
    
    
    if isinstance(bins, int):
        all_ISIs = list(itertools.chain.from_iterable(trains_ISIs))
        bins = np.linspace(min(all_ISIs), max(all_ISIs), nr_bins+1)
    
    for i, train_ISIs in enumerate(trains_ISIs):
        ISI_hists[i] = np.histogram(train_ISIs, bins)[0]
        
    return [np.mean(ISI_hists,0), np.std(ISI_hists, 0), bins]

#------------------------------------------------------------------------------ 
def cov_from_P_prob(k, N, P_probs):
#------------------------------------------------------------------------------ 
    """
    Returns the k-th order of covariance among N spike trains with given P_prob
    """
    import gmpy
    
    assert k <= N, "you can only calculate up to N-th order of covariace"
    assert len(P_probs) <= N+1, "Invalid P_prob"
    
    P_probs = np.array(P_probs)/sum(P_probs)
    mean = sum([1.0*i*P_prob for i,P_prob in enumerate(P_probs)])/N
    
    total_cov = 0
    for i in xrange(k+1):
        raw_corr = raw_corr_from_P_prob(k-i, N, P_probs)
        contribution = ((-mean)**i)*gmpy.bincoef(k,i)*raw_corr 
        total_cov += contribution 
        #print raw_corr, total_cov, contribution
    
    return float(total_cov) 
    
#------------------------------------------------------------------------------ 
def raw_corr_from_P_prob(k, N, P_probs):
#------------------------------------------------------------------------------ 
    """
    Returns the k-th order of raw correlation among N spike trains with given 
    P_prob. <x_1, x_2, ..., x_k>, k<=N.
    """
        
    assert len(P_probs) <= N+1, "Invalid P_prob"
    
    if k >= len(P_probs):
        return 0.0
    
    raw_corr = 0.0
    for P in xrange(k, len(P_probs)):
        multiplicator = 1.
        for i in xrange(k):
            #print "III"
            multiplicator *= 1.0*(P-i)/(N-i)
        #print P, P_probs[P], multiplicator
        raw_corr += P_probs[P] * multiplicator
    
    return raw_corr
    
