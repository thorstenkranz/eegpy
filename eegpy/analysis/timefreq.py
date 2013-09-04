#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Classes and funtions for calculation of EventRelated Potentials."""

import os.path

import eegpy
from eegpy.misc import FATALERROR, EegpyBase, debug
from eegpy.analysis.wavelet import wt_analyze
from eegpy.filter.freqfilt import filtfilt_low
from phases import phase_coherence
import logging
logger = logging.getLogger(__name__)

try:
    import numpy as np
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


###########
# Classes #
###########

class EventRelatedTimeFrequencyAnalyzerBase(EegpyBase):
    """Base class for time-frequency analysis. Implements basic functionality 
    used for both ERP and Wavelet-Analysis."""
    
    def __init__(self, f=None, start_end=None, baseline=None):
        self._datadict = {} #Keys are conditions, values numpy-arrays
        self._times = {} #Keys are conditions, values are list of timepoints
        self._is_calculated = {} #Keys are conditions, values are bools: Where the data already calculated? 
        self._is_artifact = {}
        
        self._set_eeg_input(f)
        self._set_start_end(start_end)
        self._set_baseline(baseline)
    
    def _set_eeg_input(self, f):
        self._eeg_was_opened_here = False
        self._dfile_name = None
        if f!=None:
            if isinstance(f, basestring):
                if os.path.exists(f):
                    self._dfile_name=f
                    self._eeg = None
                else:
                    raise IOError("File not found")
            else:
                try:
                    assert(f.numDatapoints >= 0)
                    self._eeg = f
                except Exception:
                    raise ValueError("f is neither a valid filename nor a eeg-file-object")
         
    def _set_start_end(self, start_end):
        if start_end is None:
            start_end = [-500,1500] 
        else:
            try:
                se = start_end
                assert_is_tuple_of_int_boundaries(se)
                self._start_end = (se[0],se[1])
            except AssertionError:
                raise ValueError("Error while setting start/end during initialisation. " + 
                                 "start: %s, end: %s" %(str(start_end[0]), str(start_end[1])))
                                 
    def _set_baseline(self, baseline):
        if baseline is None:
            baseline = [-200,0]
        try:
            bl = baseline
            assert_is_tuple_of_int_boundaries(bl)
            assert bl[0] >= self._start_end[0]
            assert bl[1] <= self._start_end[1]
            self._baseline = bl
        except Exception:
            start_end = self._start_end
            raise ValueError("Error while setting baseline during initialisation." +
                             "start: %s, end: %s " % (str(start_end[0]), str(start_end[1])) +
                             "baseline[0]: %s, baseline[1]: %s " % (str(bl[0]), str(bl[1])) )

    def __del__(self):
        self.close()   
    
    def __getitem__(self,key):
        #if debug:
        #    print key
        if self._is_calculated.has_key(key):
            if self._is_calculated[key] and self._datadict.has_key(key):
                return self._datadict[key]
            else:
                self.calculate(key)
                if debug:
                    print self._datadict
                return self._datadict[key]
    
    def add_condition(self,cond_name,ts):
        """Adds a named condition for ERP-analysis. If condition already exists, it's updated."""
        self._is_calculated[cond_name] = False
        if self._times.has_key(cond_name) and debug:
            print "Condition %s is being updated"%cond_name
        try:
            self._times[cond_name] = [int(x) for x in ts]
            self._is_artifact[cond_name] = np.zeros((len(self._times[cond_name])),np.bool)
        except Exception, e:
            del self._is_calculated[cond_name]
            if debug:
                print "Error while assigning timepoints for condition."
            raise e
        
    def check_timepoints(self, remove_false=False):
        """checks all timepoints of all conditions if they make sense.
        I.e. if they are within the interval [0,numdatapoints[
        Parameters:
        remove_false: Triggers wether false timepoints are removed from the conditions."""
        f = self.eeg
        ndp = f.numDatapoints
        all_ok = True
        for k in self.keys():
#            for i in range(len(self._times[k])):
#                if self._times[k][i]+self._start_end[0]<0 or self._times[k][i]+self._start_end[1]>=ndp:
#                    all_ok = False
#                    if remove_false:
#                        print "eegpy.analysis.timefreq: Removal of invalid time points not implemented yet!"
                        #del self._times[k][i]
             for i in range(len(self._times[k])-1,-1,-1): #Go through indices in inverse order
                if self._times[k][i]+self._start_end[0]<0 or self._times[k][i]+self._start_end[1]>=ndp:
                     all_ok = False
                     if remove_false:    
                         del self._times[k][i]
             self.add_condition(k, self._times[k])
        return all_ok
    
    def check_for_artifacts(self, thres1=100.0, thres2=50.0, exclude_channels=[]):
        f = self.eeg
        channels = range(f.num_channels)
        for c in exclude_channels:
            channels.remove(c)
        for k in self.keys():
            data = f.get_data_for_events(self._times[k],self._start_end)
            for tr in range(data.shape[2]):
                for ch in channels:#.shape[1]):
                    if abs(data[:,ch,tr]).max() > thres1 or abs(np.diff(data[:,ch,tr])).max() > thres2:
                        if debug:
                            print "Channel %i, Trial %i,condition %s identified as artifact" % (ch,tr,k)
                        self._is_artifact[k][tr] = True
                        break
                    
    def _get_baseline_data(self, data):
        return data[-self._start_end[0]+self._baseline[0]:-self._start_end[0]+self._baseline[1],...]
    
    def export_data(self,basename="exported_data",conditions=None):
        if conditions == None:
            conditions = self.keys()
        if self.check_timepoints():
            f = self.eeg
            for c in conditions:
                if c in self.keys():
                    data = f.get_data_for_events(self._times[c],self._start_end)
                    for tr in range(data.shape[2]):
                        fn = "%s_%s_%i.dat"%(basename,c,tr)
                        data[:,:,tr].tofile(fn)
    
    def keys(self):
        return self._times.keys()
    
    def get_channel_names(self):
        return self.eeg.getChannelNames()
    
    def get_num_valid_tps(self):
        rv = {}
        for cond in self.keys():
            rv[cond] = int(len(self._is_artifact[cond])-self._is_artifact[cond].sum())
        return rv
    
    channel_names = property(get_channel_names)
    num_valid_tps = property(get_num_valid_tps)
    
    @property
    def event_times(self):
        """Dictionary condition -> event-times"""
        return self._times
    
    @property
    def times_to_trigger(self):
        Fs = self.eeg.Fs
        start, end = self._start_end
        start_time = float(start) / Fs
        end_time = float(end) / Fs
        return np.linspace(start_time, end_time, end-start)
        
    
    @property
    def eeg(self):
        if self._eeg is None:
            if len(self._dfile_name)>0:
                self._eeg = eegpy.open_eeg(self._dfile_name)
                self._eeg_was_opened_here = True
            else:
                return None
        if not self._eeg.is_open:
            self.eeg = eegpy.open_eeg(self._dfile_name)
            self._eeg_was_opened_here = True
        return self._eeg
        
    @property
    def Fs(self):
        return self.eeg.Fs

    def calculate(self, condition):
        pass
    
    def calculate_all(self):
        for c in self.keys():
            self.calculate(c)
    
    def close(self):
        if self._eeg_was_opened_here:
            self._eeg.close()
        

class ERPAnalyzer(EventRelatedTimeFrequencyAnalyzerBase):
    """Used for the creation of ERP-Data."""
    
    def __init__(self, f=None, start_end=None, baseline=None, lp_freq=None):
        EventRelatedTimeFrequencyAnalyzerBase.__init__(self, f, start_end, baseline)
        self._lp_freq = lp_freq
        
    def calculate(self, cond_name=None):
        """After doing the setup, this method actually calculates the ERPs. 
        If no argument is supplied, then the data for all conditions are calculated."""
        if cond_name == None:
            cnames = self._times.keys()
        else:
            cnames = [cond_name]
        
        #self.check_for_artifacts()
        
        for c in cnames:
            logger.debug("Condition %s" % c)
            #assert self._datadict.has_key(c) and self._is_calculated.has_key(c), "Key in dictionary missing!"
            f = self.eeg
            #TODO: check timepoints
            ndp = f.numDatapoints
            tps_ok = True
            logger.debug("datapoints: %s, input: %s" % (ndp, str(f)))
            for t in self._times[c]:
                #if debug:
                #    print "t = ", t
                if t+self._start_end[0]<0 or t+self._start_end[1]>ndp:
                    tps_ok=False
                    break
            if tps_ok:
                #getData
                data = f.get_data_for_events(self._times[c],self._start_end)
                if self._lp_freq!=None:
                    data = filtfilt_low(self._lp_freq,data,Fs=f.Fs)
                #Baseline-Korrektur
                data -= self._get_baseline_data(data).mean(axis=0)
                #calculate average
                self._datadict[c] = data[:,:,~self._is_artifact[c]].mean(axis=2)
                #print "Shape datadict:", self._datadict[c].shape
                self._is_calculated[c] = True
            else:
                raise ValueError("Some problem with the timepoints occured." +
                                 "Number of datapoints: %i, Max. timepoint: %i" %(ndp, max(self._times[c])))

class PhaseCoherenceAnalyzer(EventRelatedTimeFrequencyAnalyzerBase):
    """Used for the creation of Phase-Coherence-Data.
       The parameter "average_over" decides wther averaging is 
       performed over trials (standard) or over a given interval 
       of each trial. In the latter case, "average_over" has to be a 
       [s,e] in samples relative to tp 0.
    """
    
    def __init__(self, fn=None, se=None, average_over="trials"):
        EventRelatedTimeFrequencyAnalyzerBase.__init__(self, fn, se)
        if not average_over=="trials":
            try:
                s,e = average_over
            except Exception:
                raise ValueError("Parameter 'average_over' not valid")
        self.average_over = average_over
        
    def calculate(self, cond_name=None):
        """After doing the setup, this method actually calculates the Phase-Coherence. 
        If no argument is supplied, then the data for all conditions are calculated."""
        if cond_name == None:
            cnames = self._times.keys()
        else:
            cnames = [cond_name]
        
        #open file
        f = eegpy.open_eeg(self._dfile_name)
        #check timepoints
        tps_ok = self.check_timepoints()
        if tps_ok:
            for c in cnames:
                if debug:
                    print "Condition %s" % c
                #assert self._datadict.has_key(c) and self._is_calculated.has_key(c), "Key in dictionary missing!"
                #getData
                data = f.get_data_for_events(self._times[c],self._start_end)
                if self.average_over == "trials":
                    #prepare _datadict[c]
                    self._datadict[c] = np.ones((data.shape[0],data.shape[1],data.shape[1]),"d")
                    #calculate phase coherence
                    for ch1 in range(data.shape[1]):
                        for ch2 in range(ch1):
                            pc = phase_coherence(data[:,ch1,~self._is_artifact[c]],data[:,ch2,~self._is_artifact[c]])
                            self._datadict[c][:,ch1,ch2] = pc[:]
                            self._datadict[c][:,ch2,ch1] = pc[:]
                else:
                    s,e = self.average_over
                    self._datadict[c] = np.ones((data.shape[1],data.shape[1],data.shape[2]),"d")
                    for ch1 in range(data.shape[1]):
                        for ch2 in range(ch1):
                            for tr in range(data.shape[2]):
                                if self._is_artifact[c][tr]:
                                    self._datadict[c][ch1,ch2,tr] = np.nan
                                    self._datadict[c][ch1,ch2,tr] = np.nan
                                else:
                                    pc = phase_coherence(data[s-self._start_end[0]:e-self._start_end[0],ch1,tr],data[s-self._start_end[0]:e-self._start_end[0],ch2,tr])
                                    #print pc
                                    #print s,e,self._start_end
                                    self._datadict[c][ch1,ch2,tr] = self._datadict[c][ch2,ch1,tr] = pc
                    #print "Shape datadict:", self._datadict[c].shape
                self._is_calculated[c] = True
        else:
            #TODO: some kind of error message
            print "Some problem with the timepoints occured."
            
class WaveletAnalyzer(EventRelatedTimeFrequencyAnalyzerBase):
    """Used for time-frequency-analysis with Morlet-Wavelets.
    Certain measures can be calculated from the wavelets, e.g.
    mean_power: The frequency-resolved power averaged over trials
    itpc: Inter-trial phase coherence
    st_power: 
    st_power_nb: No baseline-correction"""
    possible_measures = ["mean_power", "itpc", "st_power","st_power_nb"] # ,"pairwise_pc"]
    
    def __init__(self, f, frequencies=None, measures = None, start_end=None, baseline=None):
        EventRelatedTimeFrequencyAnalyzerBase.__init__(self, f, start_end, baseline)
        self._set_frequencies(frequencies)
        self._set_measures(measures)
        self.st_power_start_end = (0,self._start_end[1]/2)
        
    def _set_frequencies(self, freqs):
        if freqs == None:
            self._freqs = np.array([2.0**x for x in range(8)])
        else:
            self._freqs = np.array(freqs)

    def _set_measures(self, measures):
        if measures==None:
            self._measures == ["mean_power", "itpc", "st_power"]
        else:
            self._measures = list(set([m for m in measures if m in self.possible_measures]))
            assert len(self._measures) > 0, "No valid measure was given"
            
    def calculate(self, cond_name=None):
        """After doing the setup, this method actually calculates the ERPs. 
        If no argument is supplied, then the data for all conditions are calculated."""
        if cond_name == None:
            cnames = self._times.keys()
        else:
            cnames = [cond_name]
        #Check if settings for single-trial-power are o.k.
        if "st_power" in self._measures or "st_power_nb" in self._measures:
            if self.st_power_start_end[0]>self.st_power_start_end[1] or self.st_power_start_end[0]<self._start_end[0] or self.st_power_start_end[1]>self._start_end[1]:
                raise ValueError("The settings for the interval for the single-trial power are incorrect: %s"%(str(self.st_power_start_end)))
        
        tps_ok = self.check_timepoints()
        f = self.eeg
        #ndp = f.num_datapoints
        start, end = self._start_end
        ndp_data = end-start
        n_f = len(self._freqs)
        if tps_ok:
            for c in cnames:
                self._datadict[c] = {}
                if "mean_power" in self._measures:
                    self._datadict[c]["mean_power"] = np.zeros((ndp_data, f.numChannels, n_f),"d")
                if "itpc" in self._measures:
                    self._datadict[c]["itpc"] = np.zeros((ndp_data, f.numChannels, n_f),"d")
                if "st_power" in self._measures:
                    self._datadict[c]["st_power"] = np.zeros((f.numChannels,len(self._times[c]),n_f),"d")
                if "st_power_nb" in self._measures:
                    self._datadict[c]["st_power_nb"] = np.zeros((f.numChannels,len(self._times[c]),n_f),"d")
                    
                for ch in range(f.numChannels):
                    logger.debug("Condition %s, Channel %s" % (c, ch))
                    #getData
                    data = f.get_data_for_events(self._times[c],self._start_end,[ch])
                    #calculate wavelets
                    wts = self._calculate_wavelets(data)
                    if "mean_power" in self._measures:
                        self._datadict[c]["mean_power"][:,ch,:] = self._calculate_mean_power(wts[:,~self._is_artifact[c],:])
                    if "itpc" in self._measures:
                        self._datadict[c]["itpc"][:,ch,:] = self._calculate_itpc(wts[:,~self._is_artifact[c],:])
                    if "st_power" in self._measures:
                        self._datadict[c]["st_power"][ch,:,:] = self._calculate_st_power(wts)
                    if "st_power_nb" in self._measures:
                        self._datadict[c]["st_power_nb"][ch,:,:] = self._calculate_st_power_nb(wts)
                #print "Shape datadict:", self._datadict[c].shape
                self._is_calculated[c] = True
        else:
            print "Some problem with timepoints."
            
    def _calculate_wavelets(self, data):
        rv = np.zeros((data.shape[0],data.shape[2],len(self._freqs)),"D")
        for j in range(data.shape[2]):
            rv[:,j,:] = wt_analyze(data[:,0,j],self._freqs, self.Fs)
        return rv
    
    def _calculate_mean_power(self, good_wts):
        tmp = (abs(good_wts)**2).mean(axis=1)
        for fr in range(tmp.shape[1]):
            unnormalized_power = 10*( np.log10(tmp[:,fr]))
            baseline_data = self._get_baseline_data(unnormalized_power)
            #print baseline_data.shape, baseline_data.min(), baseline_data.mean(), baseline_data.max()
            tmp[:,fr] = unnormalized_power - baseline_data.mean()
            #print self._get_baseline_data(tmp[:,fr]).mean()
        return tmp
            
    def _calculate_itpc(self, good_wts):
        return abs((good_wts[:,:,:] / abs(good_wts[:,:,:])).mean(axis=1))
    
    def _calculate_st_power(self, wts):
        tmp = abs(wts[:,:,:])**2
        for tr in range(tmp.shape[1]):
            for fr in range(tmp.shape[2]):
                unnormalized_power = 10*( np.log10(tmp[:,tr,fr]))
                tmp[:,tr,fr] = unnormalized_power-self._get_baseline_data(unnormalized_power).mean()
        return tmp[-self._start_end[0]+self.st_power_start_end[0]:-self._start_end[0]+self.st_power_start_end[1],:,:].mean(axis=0)
        
    def _calculate_st_power_nb(self, wts):
        tmp = abs(wts[:,:,:])**2
        return tmp[-self._start_end[0]+self.st_power_start_end[0]:-self._start_end[0]+self.st_power_start_end[1],:,:].mean(axis=0)
    
    @property
    def frequencies(self):
        return self._freqs

class WaveletPowerAnalyzer(WaveletAnalyzer):
    """Reduced functionality for easier access to the most relevant analysis of
    WaveletAnalyzer"""
    
    def __init__(self, f, frequencies=None,  start_end=None, baseline=None):
        WaveletAnalyzer.__init__(self, f, frequencies=frequencies,
                                      measures = ["mean_power"],
                                      start_end=start_end, baseline=baseline)
    
    def __getitem__(self, key):
        return WaveletAnalyzer.__getitem__(self, key)["mean_power"]
        
        
def assert_is_tuple_of_int_boundaries(boundaries):
    assert boundaries!=None
    assert len(boundaries) == 2
    assert type(boundaries[0])==type(1000)
    assert type(boundaries[1])==type(1000)
    assert boundaries[0]<boundaries[1]
        
        
        
        
        
