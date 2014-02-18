"""
Generating MIP (Kuhn, Aertsen & Rotter, 2003):
1 A poisson process is generated as the mother process. Note that r is the 
  total rate of the MIP, and the rate of the mother process is r/(n*corr2). 
  All temporal dimensions are in ms.
2 The cluster size is drawn from a binomial distribution.
3 spk_times stores the spike times in an array.

Code by Yim, Man Yi
Refractored to fit into NeuroVivo environment by Miha Pelko
"""

import numpy as np
import correlations.correlated_input as ci

# -----------------------------------------------------------------------------
def create_n_trains(n, mean, corr2, time, np_rnd=None):
# -----------------------------------------------------------------------------
	"""
	Creates the trains according to MIP algorithm and then discretizes them in 
	the bins.
	"""
	spk_times = create_n_trains_cont(n, mean, corr2, time, np_rnd)
	return [ci.discretize_train(spk_times[i], time, 1) for i in xrange(n)]


# -----------------------------------------------------------------------------
def create_n_trains_cont(n, mean, corr2, time, np_rnd=None):
# -----------------------------------------------------------------------------
	if not np_rnd:
		np_rnd = np.random 
	
	assert corr2>0, "MIP can only create correlated trains."
	
	num_spike = np_rnd.poisson(int(round(mean*time/corr2)))

	proc_gm = np_rnd.rand(num_spike)*time # Mother process
	proc_gm.sort()
	spk_times = [[] for i in xrange(n)]
	for spike_time in proc_gm:
		copy_to_how_many = np_rnd.binomial(n, corr2)
		if copy_to_how_many:
			train_ids = range(n)
			np_rnd.shuffle(train_ids)
			trains_getting_spike = train_ids[:copy_to_how_many]
			for train in trains_getting_spike:
				spk_times[train].append(spike_time)
	return spk_times