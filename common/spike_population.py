'''
Class SpikePopulation for methods on a population of neurons with the spike 
trains defined by the spike times and optionally by the total time. 
'''
import numpy as np
import neurovivo.analysis.spike_train_analysis as sta
import random as rnd
from neurovivo.common.spike_train import SpikeTrain

class SpikePopulation(object):
  
    def __init__(self, spike_trains, total_time=None, population_size=None, name = "Unnamed"):
        """
        Class SpikePopulation containing the spike trains in the population. 
        All times are presumed to be in ms.
        """
        self._name = name
        self.spike_trains=spike_trains

        # Setting the total_time
        times = [train._total_time for train in self.spike_trains]
        if not len(times) == 0:
            time = np.max(times)
        else:
            time = 0
            
        if not total_time == None:
            assert total_time >= time, "the specified total_time ({0}) seems to be too low. Found a spike at time = {1}".format(total_time, time)
        else: 
            total_time = time
        self._total_time = total_time

        # Setting the population_size
        if not population_size == None:
            assert population_size >= len(self.spike_trains)
        else:
            self.population_size = len(self.spike_trains)

    def rates(self):
        """
        Returns the average rates for the spike trains in the population.
        If there are no spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 1:
            return np.NaN
        else:
            return [train.average_rate() for train in self.spike_trains]
            
    def rate(self):
        """
        Returns the mean of the average rate for the spike trains in the population.
        If there are no spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 1:
            return np.NaN
        else:
            return np.mean([train.average_rate() for train in self.spike_trains])

    def rate_std(self):
        """
        Returns the std of the average rate for the spike trains in the population.
        If there are less thent 2 spike trains present, it returns NaN.
        """
        if len(self.spike_trains) < 2:
            return np.nan
        else:
            return np.std([train.average_rate() for train in self.spike_trains])

    def CV_ISIs(self):
        """
        Returns the coefficients of variation of interspike intervals
        for the spike trains in the population.
        If there are no spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 1:
            return np.NaN
        else:
            return [train.CV_ISI() for train in self.spike_trains]
                    
    def CV_ISI(self):
        """
        Returns the average coefficient of variation of interspike intervals
        for the spike trains in the population.
        If there are no spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 1:
            return np.NaN
        else:
            CV_ISIs = np.array(self.CV_ISIs())
            # ignore the trains with too few spikes (SpikeTrain.CV_ISI method returns -1 for those)
            CV_ISIs = CV_ISIs[CV_ISIs >= 0]
            return np.mean(CV_ISIs)
   
    def CV_ISI_std(self):
        """
        Returns the standard deviation of the coefficient of variation of 
        interspike intervals for the spike trains in the population.
        If there are less then 2 trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 2:
            return np.NaN
        else:
            CV_ISIs = np.array(self.CV_ISIs())
            # ignore the trains with too few spikes (SpikeTrain.CV_ISI method returns -1 for those)
            CV_ISIs = CV_ISIs[CV_ISIs >= 0]
            return np.std(CV_ISIs)

    def local_variations(self):
        """
        Returns the local variations of interspike intervals
        for the spike trains in the population.
            Shigeru Shinomoto et al., "Relating Neuronal Firing Patterns 
            to Functional Differentiation of Cerebral Cortex," 
            PLoS Comput Biol 5, no. 7 (July 10, 2009): e1000433.
        If there are less then 2 spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 2:
            return np.NaN
        else:
            return [train.local_variation() for train in self.spike_trains]

    def local_variation(self):
        """
        Returns the average local variation of interspike intervals
        for the spike trains in the population.
            Shigeru Shinomoto et al., "Relating Neuronal Firing Patterns 
            to Functional Differentiation of Cerebral Cortex," 
            PLoS Comput Biol 5, no. 7 (July 10, 2009): e1000433.
        If there are less then 2 spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 2:
            return np.NaN
        else:
            lvs = np.array(self.local_variations())
            # ignore the trains with too few spikes (SpikeTrain.local_variation method returns -1 for those)
            lvs = lvs[lvs >= 0]
            return np.mean(lvs)

    def local_variation_std(self):
        """
        Returns the std of the local variation of interspike intervals
        for the spike trains in the population.
            Shigeru Shinomoto et al., "Relating Neuronal Firing Patterns 
            to Functional Differentiation of Cerebral Cortex," 
            PLoS Comput Biol 5, no. 7 (July 10, 2009): e1000433.
        If there are less then 2 spike trains present, it returns np.NaN.
        """
        if len(self.spike_trains) < 2:
            return np.NaN
        else:
            lvs = np.array(self.local_variations())
            # ignore the trains with too few spikes (SpikeTrain.local_variation method returns -1 for those)
            lvs = lvs[lvs >= 0]
            return np.std(lvs)
        
    def pairwise_correlations(self, windows, n=50, corr_type="dynamic_rate", dt=1):
        """
        Returns an array of n arrays of len(windows) results. windows is a parameter
        of window widths used for calculating correlation, n is the number of pairs
        used for calculation, corr_type is string for one of the following options:
        - dynamic_rate
        - pearson
        - spike_count
        """
        corr_functions={"dynamic_rate":sta.correlation, 
                        "pearson":sta.pearson_correlation, 
                        "spike_count":sta.spike_count_correlation}
        
        #print "starting pairwise " + corr_type
        correlations = [[] for _ in xrange(n)]
        for i in xrange(n):
            [st1, st2] = rnd.sample(self.spike_trains,2)
            for w in windows:
                correlations[i].append(corr_functions[corr_type](st1, st2, dt=dt, window=w))
        return correlations

    def nth_correlation(self, windows, nth_corr, nr_combinations=50, corr_type="dynamic_rate", dt=1, statistics=True):
        """
        Returns the nth order of correlation.
        Windows is a parameter of window widths used for calculating correlation,
        n is the number of combinations used for calculation, corr_type is string 
        for one of the following options:
        - dynamic_rate
        - spike_count
        - normalized
        If statistics = True, it returns the [mean and std] of the results, otherwise it
        returns the results for all nr_combinations.
        """

        print "starting nth_correlation: where n = " + str(nth_corr) + "; type = " + corr_type

        return sta.nth_order_correlation(self.spike_trains, nth_corr, dt=dt, windows=windows, nr_combinations=nr_combinations, corr_type=corr_type, statistics=statistics)
        
    def as_single_spike_train(self):
        """
        Returns a single SpikeTrain containing all the spikes in the population.
        """
        from neurovivo.common.spike_train import SpikeTrain
        spike_times = []
        for train in self.spike_trains:
            spike_times.extend(train._spike_times)
        spike_times.sort()
        return SpikeTrain(spike_times,total_time=self._total_time)
    
    def population_part(self, t_min, t_max):
        """
        Returns a SpikePopulation object including only the spikes between t_min and t_max.
        The spike times are rebased (new_time = old_time - t_min), total_time of the population
        is t_max-t_min.
        """
        assert t_max > t_min
        total_time = t_max-t_min
        new_trains = []
        for train in self.spike_trains:
            times = train._spike_times
            times = times - t_min
            times = times[times>0]
            times = times[times<total_time]
            new_trains.append(SpikeTrain(times,total_time=total_time))
        return SpikePopulation(new_trains)