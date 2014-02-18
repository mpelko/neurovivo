'''
Created on Jun 17, 2011

@author: mpelko
'''
import numpy as np

class Epsp(object):
    NOT_AVAILABLE = np.NaN
  
    def __init__(self, trace, dt, stimulus_onset_ind):
        
        assert dt*stimulus_onset_ind > 0.044,\
        "The epsp trace needs to start at least 44 ms before the stimulus_onset."
        
        self.mean_win_ind = int(0.04/dt) # the time window before stimulus used to rebase the trace
        
        start_trace = trace[stimulus_onset_ind-self.mean_win_ind:stimulus_onset_ind]
        self._trace = trace - np.mean(start_trace)
        self._dt = dt
        self._stimulus_onset_ind = stimulus_onset_ind
        self.reset_features()

# Getter public methods ------------------------------------
    def rebase(self, delta):
        """
        adds delta to the trace. 
        """
        self._trace = self._trace + delta
        self.reset_features()

    def reset_features(self):
        """
        Resets all the current feature evaluations.
        """
        self._onset_ind = None
        self._offset_ind = None
        self._peak = None
        self._peak_ind = None
        self._integral = None
        self._rise_slope = None
        self._noise_std = None
        self._rise_time = None

    def onset(self):
        """
        Returns the start time of the epsp rise.
        The threshold is based on the epsp peak value.
        """
        if self._onset_ind == None:
            self._set_onset_offset()
        if self._onset_ind == Epsp.NOT_AVAILABLE or np.isnan(self._onset_ind):
            return self._onset_ind
        else:
            return self._onset_ind * self._dt

    def offset(self):
        """
        Returns the end time of the epsp.
        The threshold is based on the epsp peak value.
        """
        if self._offset_ind == None:
            self._set_onset_offset()
        if self._offset_ind == Epsp.NOT_AVAILABLE or np.isnan(self._offset_ind):
            return self._offset_ind
        else:
            return self._offset_ind * self._dt

    def peak(self):
        """
        Returns the peak voltage of epsp.
        """
        if self._peak == None:
            self._set_peak()
        return self._peak

    def peak_time(self):
        """
        Returns the timing of the peak of the epsp.
        """
        if self._peak_ind == None:
            self._set_peak()
        if self._peak_ind == Epsp.NOT_AVAILABLE or np.isnan(self._peak_ind):
            return self._peak_ind
        else:
            return self._peak_ind * self._dt

    def rise_slope(self):
        """
        Returns the slope of the epsp rise.
        Defined simply by the slope of the line connecting the points at 1/5 and 
        at 4/5 of the epsp peak voltage.
        """
        if self._rise_slope == None:
            self._set_rise_slope()
        return self._rise_slope

    def rise_time(self):
        """
        Returns the time of epsp rising.
        Defined simply by the time between the 1/5 and 
        at 4/5 of the epsp peak voltage.
        """
        if self._rise_time == None:
            self._set_rise_time()
        return self._rise_time
        
    def integral(self):
        """
        Returns the integral of EPSP between the onset and offset time.
        """
        if self._integral == None:
            self._set_integral()
        return self._integral

    def noise_std(self):
        """
        Returns the std of the signal (noise level), based on the 
        std in the time before the stimulation.
        """
        if self._noise_std == None:
            self._set_noise_std()
        return self._noise_std

    def is_good_EPSP(self):
        """
        Checks if the trace really is an EPSP. Returns the quality 
        (Good/Bad/Maybe) and the message with the explanation.
        """
        # peak is too small
        if self.peak() == 0:
            return ["Bad", "Peak is too small"]
        if self.integral() == Epsp.NOT_AVAILABLE or np.isnan(self.integral()):
            return ["Bad", "Could not calculate the integral."]
        # very small EPSP length
        if self.offset()-self.onset() < 0.002:
            return ["Bad", "Too small EPSP length"]
        # very small peak
        if self.peak() <= 5*self.noise_std():
            return ["Maybe", "Relatively small peak"]
        # small EPSP length 
        if self.offset()-self.onset() < 0.01:
            return ["Maybe", "EPSP very short"]
        # too long EPSP
        if self.offset()-self.onset() >= 0.1:
            return ["Maybe", "EPSP very long"] 
        else:
            return ["Good", "No complaints"]

# Helper methods (calculating epsp features) ----------------------------
    def _check_if_good_EPSP(self, feature):
        if self.peak() == 0:
            self.feature = Epsp.NOT_AVAILABLE 
            return False
        return True
            
    def _set_onset_offset(self):
        """
        Determening the start time of the epsp rise.
        The threshold is based on the epsp peak value.
        """
        if not self._check_if_good_EPSP(self._onset_ind) and not self._check_if_good_EPSP(self._offset_ind):
            return
        threshold = max(self.peak()/5.,0.0000001)

        after_stimulus_trace = self._trace[self._stimulus_onset_ind:self._peak_ind]
        #print 
        if len(np.nonzero(after_stimulus_trace < threshold)[0]) == 0:
            onset_ind = 0
        else:
            #print np.nonzero(after_stimulus_trace < threshold)[0]
            onset_ind = np.nonzero(after_stimulus_trace < threshold)[0][-1] + 1
        self._onset_ind = self._stimulus_onset_ind + onset_ind
       
        # the offset has to occur after the peak
        if self._peak_ind == None:
            self._set_peak()
        after_peak_trace = self._trace[self._peak_ind:]
        # sometimes the trace never goes bellow the threshold
        if np.min(after_peak_trace >= threshold):
            self._offset_ind = Epsp.NOT_AVAILABLE
            self._integral = Epsp.NOT_AVAILABLE
            return

        offset_ind = np.nonzero(after_peak_trace <
                                threshold)[0][0]
        self._offset_ind = self._peak_ind + offset_ind
    
    def _set_noise_std(self):
        """
        Determines the std of the signal (noise level), based on the 
        std in the time before the stimulation.
        """
        time_before_stim = int(self.mean_win_ind/self._dt)
        std_trace = self._trace[self._stimulus_onset_ind - time_before_stim:self._stimulus_onset_ind]
        
        self._noise_std = np.std(std_trace)
    
    def _set_peak(self):
        """
        Determening the peak value of the epsp. A peak is the first significant 
        voltage maximum after the simulation onset.
        """
        noise_limit = 3    # how many std over noise do we still consider it a peak
        time_window = 0.07 # time window of interest after the stimulus onset
        
        ### now limiting the expected time for the peak to some window after onset 
        tw_end_ind = self._stimulus_onset_ind + int(time_window/self._dt) 
        after_stimulus_trace = self._trace[self._stimulus_onset_ind:tw_end_ind]
        peak = np.max(after_stimulus_trace)
        peak_ind = np.nonzero(after_stimulus_trace == peak)[0][0]
        
        if peak < noise_limit * self.noise_std():
            peak = 0
            peak_ind = Epsp.NOT_AVAILABLE
            self._onset_ind = Epsp.NOT_AVAILABLE
            self._offset_ind = Epsp.NOT_AVAILABLE
            self._integral = Epsp.NOT_AVAILABLE
            self._rise_slope = Epsp.NOT_AVAILABLE
            self._peak = peak
            self._peak_ind = peak_ind
            return 
        
        else:
            ### FIND THE FIRST peak that exceeds the noise_limit 
            
            # TODO: think how to optimize this:
            #     - idea: deducting shifts to both sides and comparing. 
            for time_point in np.arange(1,len(after_stimulus_trace)-2):
                if after_stimulus_trace[time_point+1] < after_stimulus_trace[time_point] and\
                   after_stimulus_trace[time_point-1] < after_stimulus_trace[time_point]:
                    peak = after_stimulus_trace[time_point]
                    if peak > noise_limit * self.noise_std():
                        self._peak = peak
                        self._peak_ind = time_point + self._stimulus_onset_ind
                        return
        # Finally we are left with ever rising EPSP.
        self._peak = after_stimulus_trace[-1]
        self._peak_ind = tw_end_ind
            
    def _set_rise_slope(self):
        """
        Defined simply by the slope of the line connecting the points at 1/5 and 
        at 4/5 of the epsp peak voltage.
        """
        if not self._check_if_good_EPSP(self._rise_slope):
            return

        peak_value = self.peak()
        start_value = peak_value/5.
        end_value = start_value*4

        after_stimulus_trace = self._trace[self._stimulus_onset_ind:]
        start_index = np.nonzero(after_stimulus_trace > start_value)[0][0]
        start_time = start_index * self._dt
        end_index = np.nonzero(after_stimulus_trace > end_value)[0][0]
        end_time = end_index * self._dt
        
        self._rise_slope = (end_value-start_value)/(end_time-start_time)

    def _set_rise_time(self):
        """
        Defined simply by the time between 1/5 and 
        4/5 of the epsp peak voltage.
        """
        if not self._check_if_good_EPSP(self._rise_time):
            return

        peak_value = self.peak()
        start_value = peak_value/5.
        end_value = start_value*4

        after_stimulus_trace = self._trace[self._stimulus_onset_ind:]
        start_index = np.nonzero(after_stimulus_trace > start_value)[0][0]
        start_time = start_index * self._dt
        end_index = np.nonzero(after_stimulus_trace > end_value)[0][0]
        end_time = end_index * self._dt
        
        self._rise_time = end_time-start_time
        
    def _set_integral(self):
        import scipy.integrate as sint
        
        if not self._check_if_good_EPSP(self._integral):
            return
       
        if self._onset_ind == None or self._offset_ind == None:
            self.onset()
        
        if self.offset() == Epsp.NOT_AVAILABLE or np.isnan(self.offset()):
            self._integral = Epsp.NOT_AVAILABLE
            return
        
        integration_trace = self._trace[self._onset_ind:self._offset_ind]
        #print self._onset_ind, self._offset_ind
        integral = sint.simps(y=integration_trace, dx=self._dt)
        self._integral = integral

# OLD: ----------------------------------------------------
    def _set_onset_old(self):
        """
        Determening the start time of the epsp rise.
        The threshold is based on the variance of the signal before the stimulus.
        """
        
        assert False, "OLD method"
        
        # what is the minimum noise std (mainly for surrogate data)
        minimal_std = 0.02
        std_times_criteria = 4  

        noise_std = self.noise_std()
        std = max(minimal_std, noise_std)
        
        after_stimulus_trace = self._trace[self._stimulus_onset_ind:]
        onset_ind = np.nonzero(after_stimulus_trace >
                                std_times_criteria * std)[0][0]
        self._onset_ind = onset_ind

    def _set_peak_old(self):
        """
        Determening the peak value of the epsp.
        """
        noise_limit = 3    # how many std over noise do we still consider it a peak
        time_window = 0.05 # time window of interest after the stimulus onset
        
        ### now limiting the expected time for the peak to some window after onset 
        tw_end_ind = self._stimulus_onset + int(time_window/self._dt) 
         
        after_stimulus_trace = self._trace[self._stimulus_onset_ind:tw_end_ind] 
        peak = np.max(after_stimulus_trace)
        peak_ind = np.nonzero(after_stimulus_trace == peak)[0][0]
        if peak < noise_limit * self.noise_std():
            peak = 0
            peak_ind = Epsp.NOT_AVAILABLE
            self._onset_ind = Epsp.NOT_AVAILABLE
            self._offset_ind = Epsp.NOT_AVAILABLE
            self._integral = Epsp.NOT_AVAILABLE
            self._rise_slope = Epsp.NOT_AVAILABLE
        self._peak = peak
        self._peak_ind = peak_ind 
