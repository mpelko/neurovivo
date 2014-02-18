'''
Class SpikeTrain for methods on a single spike train defined by the spike times
and optionally by the total time. 
'''
import numpy as np
import neurovivo.analysis.spike_train_analysis as sta

class SpikeTrain(object):
  
    def __init__(self, spike_times, start_time=0, total_time=None, name="Unnamed"):
        """
        Class SpikeTrain containing the spike timings of the train. All times
        are presumed to be in ms.
        """
        if not len(spike_times) == 0:
            assert start_time <= np.min(spike_times), "{0}".format(np.min(spike_times))
        
        self._name = name
        #just in case we also sort them. Possible to optimise later ...
        self._spike_times = np.sort(spike_times)
        if total_time == None:
            if len(spike_times) == 0:
                total_time = 0
            else:
                total_time = np.max(spike_times) + start_time
        self._total_time = total_time
        self._start_time = start_time
   
    def average_rate(self):
        return 1000.*len(self._spike_times)/self._total_time
    
    def spike_times(self):
        return self._spike_times
    
    # legacy code
    def averageRate(self):
        return self.average_rate()
    
    def ISIs(self):
        if len(self._spike_times) < 2:
            return []
        else:
            return self._spike_times[1:] - self._spike_times[:-1]
    
    def nr_spikes_between(self,t_min,t_max):
        '''
        Returns the number of spikes in the interval [t_min, t_max)
        '''
        assert t_min < t_max
        A = self._spike_times>=t_min
        B = self._spike_times<t_max
        return np.sum(A*B)

    def spike_count(self, time, window = 20):
        """
        returns the number of spikes within the desired window at
        a specific time.
        """
        assert time >=self._start_time and time < self._total_time # making sure
        half_win = window/2.
        return self.nr_spikes_between(time-half_win, time+half_win)

    def moving_average(self, time, window = 20):
        """
        returns the moving average window value for the desired window size at
        a specific time.
        """
        #assert time >=self._start_time and time <= self._total_time, "{}, {}".format(time,self._total_time) # making sure
        half_win = window/2.
        divider = 1. * window
        #if time < (self._start_time + half_win):
        #    divider = time + half_win
        #elif time > (self._total_time-half_win):
        #    divider = (self._start_time+self._total_time)-time+half_win
        return 1000*self.nr_spikes_between(time-half_win, time+half_win)/divider
    
    def CV_ISI(self):
        """
        Rreturns the CV (coefficient of variation) of ISI (inter-spike
        interval) for the spike train. 
        If there are less then 3 spikes present (less then 2 ISIs) it returns None.
        """
        ISIs = self.ISIs()
        if len(ISIs) < 2:
            return np.nan
        else:
            return np.std(ISIs)/np.mean(ISIs)

    def local_variation(self):
        """
        Rreturns the coefficient of local variation for the spike train. The idea taken from the paper:
        Shigeru Shinomoto et al., "Relating Neuronal Firing Patterns to Functional Differentiation of Cerebral Cortex," PLoS Comput Biol 5, no. 7 (July 10, 2009): e1000433.
        If there are less then 3 spikes present (less then 2 ISIs) it returns -1.
        """
        ISIs = self.ISIs()
        if len(ISIs) < 2:
            return np.nan
        else:
            return 3.*np.sum((1.*(ISIs[:-1]-ISIs[1:])/(ISIs[:-1]+ISIs[1:]))**2)/(len(ISIs)-1)
    
    # TODO: Merge the dynamic_rate and dynamic_spike_count functions
    # they are essentially the same, only normalised by window
    # All the fast methods are now redundant as they have been 
    # substituted by the dynamic_calculation method.
    def dynamic_calculations(self, dt=1, window=20, calculation="spike_count"):
        """
        returns the dynamic rate (calucaltion="spike_count") or dynamic spike
        count (calculation="rate") for the spike train using a given window 
        with the chosen dt step.
        """
        assert dt<=window, "Having smaller window then dt. Possibly losing spikes. Not good." 
        assert window%dt==0, "Not permitted. Write your own routine."
        if not window % (2*dt) == 0:
            bins = np.arange(self._start_time-dt/2., self._start_time + self._total_time + dt, dt)
        else:
            bins = np.arange(self._start_time, self._start_time + self._total_time + dt, dt)
        discrete=np.histogram(self._spike_times, bins)[0]
        result = 0
        if dt==window:
            result=discrete
        else:
            kernel = np.ones(window/dt)
            cutoff = (int(window/dt)-1)/2
            result = np.convolve(discrete,kernel)
            result = result[cutoff:(len(result)-cutoff)] 
        if calculation=="rate":
            result = 1000.*result/window
        return result
    
    def dynamic_rate_slow(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with the chosen dt step.
        """
        time_points = np.arange(((self._total_time-self._start_time) / dt) + 1) * dt
        return np.array([self.moving_average(time,window) for time in time_points])
    
    def dynamic_rate_fast(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with 3the chosen dt step.
        """
        assert window < 250 # for bigger windows are slower then the slow method
        assert (self._total_time-self._start_time) % dt == 0, "Bad dt."
        assert dt<=window, "Having smaller window then dt. Possibly losing spikes. Not good." 
        assert (window/2.) % dt==0
        time_points = np.arange(((self._total_time) / dt) + 1) * dt + self._start_time
        binned_counts = np.histogram(self.spike_times(), time_points)[0]
        binned_counts_pp = np.r_[[0]*(window/2),binned_counts,[0]*(window/2)]
        result = np.sum([binned_counts_pp[i:i+len(time_points)] for i in xrange(window)],0)
        return 1000.*result/window
        
    def dynamic_rate(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with the chosen dt step.
        """
        try:
            return self.dynamic_calculations(dt, window, calculation="rate")
        except:
            try: 
                return self.dynamic_rate_fast(dt=dt, window=window)
            except:
                return self.dynamic_rate_slow(dt=dt, window=window)
            
    def dynamic_spike_count_slow(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with the chosen dt step.
        """
        time_points = np.arange(self._start_time,self._total_time,dt)
        return np.array([self.spike_count(time,window) for time in time_points])

    def dynamic_spike_count_fast(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with the chosen dt step.
        """
        assert window < 250 # because bigger windows are slower then the slow method
        assert (self._total_time-self._start_time) % dt == 0, "Bad dt."
        assert dt<=window, "Having smaller window then dt. Possibly losing spikes. Not good." 
        assert (window/2.) % dt==0
        time_points = np.arange(((self._total_time) / dt) + 1) * dt + self._start_time
        binned_counts = np.histogram(self.spike_times(), time_points)[0]
        binned_counts_pp = np.r_[[0]*(window/2),binned_counts,[0]*(window/2)]
        result = np.sum([binned_counts_pp[i:i+len(time_points)] for i in xrange(window)],0)
        return result    

    def dynamic_spike_count(self, dt=1, window=20):
        """
        returns the dynamic rate for the spike train using a given window with the chosen dt step.
        """
        try:
            return self.dynamic_calculations(dt, window, "spike_count")
        except:
            try: 
                return self.dynamic_spike_count_fast(dt=dt, window=window)
            except:
                return self.dynamic_spike_count_slow(dt=dt, window=window)

    def conductance_trace(self, tau1, tau2, dt=0.1, kernel_length=250):
        """
        Going form spike trains to conductance changes in time by convolving the 
        spike train with the exponential-exponential kernel.
        st -            SpikeTrain
        tau1, tau2 -    time constants for the kernel
        dt -            sampling time in ms (should be smaller then tau1 and tau2)
        kernel_length - the length of kernel in ms
        
        Returns a SimpleTrace (so the info on start_time is lost).
        """
        return sta.convolve_with_expexp_kernel(self, tau1, tau2, dt, kernel_length)

    def __str__(self):
        return "SpikeTrain object: " + self._name
