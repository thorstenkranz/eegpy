#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Classes and funtions for calculation of EventRelated Potentials."""

import os.path
import pickle
import time

import eegpy
from eegpy.misc import FATALERROR, EegpyBase, debug
from eegpy.analysis.wavelet import make_wavelets,wt_analyze
from eegpy.filter.freqfilt import filtfilt_low
from phases import phase_coherence

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


###########
# Classes #
###########

class TFBaseAnalyzer(EegpyBase):
    """Base class for time-frequency analysis. Implements basic functionality 
    used for both ERP and Wavelet-Analysis."""
    
    def __init__(self, fn=None, se=None):
        self._datadict = {} #Keys are conditions, values numpy-arrays
        self._times = {} #Keys are conditions, values are list of timepoints
        self._is_calculated = {} #Keys are conditions, values are bools: Where the data already calculated? 
        self._is_artifact = {}
        #if debug:
        #    print "_is_calculated: ", self._is_calculated
        self._dfile_name = None
        self._start_end = [-500,1500] # relative timepoints for ERPs
        if fn!=None:
            fn = str(fn)
            if os.path.exists(fn):
                self._dfile_name=fn 
            else: 
                raise IOError("File not found")
        try:
            if se!=None and type(se[0])==type(1000) and type(se[1])==type(1000) and se[0]<se[1]:
                self._start_end = (se[0],se[1])
        except Exception, e:
            if debug:
                print "Error while setting start/end during initialisation. Using standard values instead."
    
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
        ok = True
        self._is_calculated[cond_name] = False
        if self._times.has_key(cond_name) and debug:
            print "Condition %s is being updated"%cond_name
        try:
            self._times[cond_name] = [int(x) for x in ts]
            self._is_artifact[cond_name] = n.zeros((len(self._times[cond_name])),n.bool)
        except Exception, e:
            ok = False
            del self._is_calculated[cond_name]
            if debug:
                print "Error while assigning timepoints for condition."
        #try:
        #    if ok:
        #        self._datadict[cond_name] = n.zeros((self._start_end[1]-self._start_end[0]),"d")
        #except Exception, e:
        #    del self._times[cond_name]
        #    del self._is_calculated[cond_name]
        #    ok = False
        #    if debug:
        #        print "Error while creating array for data."
        
    def check_timepoints(self, remove_false=False):
        """checks all timepoints of all conditions if they make sense.
        I.e. if they are within the interval [0,numdatapoints[
        Parameters:
        remove_false: Triggers wether false timepoints are removed from the conditions."""
        ndp = 10e50
        try:
            ndp = f.numDatapoints
        except Exception, e:
            try: 
                f.close()
            except Exception,e:
                pass
            f = eegpy.open_eeg(self._dfile_name)
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
        f = eegpy.open_eeg(self._dfile_name)
        channels = range(f.num_channels)
        for c in exclude_channels:
            channels.remove(c)
        for k in self.keys():
            data = f.get_data_for_events(self._times[k],self._start_end)
            for tr in range(data.shape[2]):
                for ch in channels:#.shape[1]):
                    if abs(data[:,ch,tr]).max() > thres1 or abs(n.diff(data[:,ch,tr])).max() > thres2:
                        if debug:
                            print "Channel %i, Trial %i,condition %s identified as artifact" % (ch,tr,k)
                        self._is_artifact[k][tr] = True
                        break
        f.close()
    
    def export_data(self,basename="exported_data",conditions=None):
        if conditions == None:
            conditions = self.keys()
        if self.check_timepoints():
            f = eegpy.open_eeg(self._dfile_name)
            for c in conditions:
                if c in self.keys():
                    data = f.get_data_for_events(self._times[c],self._start_end)
                    for tr in range(data.shape[2]):
                        fn = "%s_%s_%i.dat"%(basename,c,tr)
                        data[:,:,tr].tofile(fn)
            f.close()
    
    def keys(self):
        return self._times.keys()
    
    def get_channel_names(self):
        rv = []
        f = eegpy.open_eeg(self._dfile_name)
        rv = f.getChannelNames()
        f.close()
        return rv
    
    def get_num_valid_tps(self):
        rv = {}
        for cond in self.keys():
            rv[cond] = int(len(self._is_artifact[k])-self._is_artifact[k].sum())
        return rv
    
    channel_names = property(get_channel_names)
    num_valid_tps = property(get_num_valid_tps)
    
    def get_times(self):
        return self._times
    
    ts = property(get_times)

class ERPAnalyzer(TFBaseAnalyzer):
    """Used for the creation of ERP-Data."""
    
    def __init__(self, fn=None, se=None, lp_freq=None):
        TFBaseAnalyzer.__init__(self, fn, se)
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
            if debug:
                print "Condition %s" % c
            #assert self._datadict.has_key(c) and self._is_calculated.has_key(c), "Key in dictionary missing!"
            #TODO: open file
            f = eegpy.open_eeg(self._dfile_name)
            #TODO: check timepoints
            ndp = f.numDatapoints
            tps_ok = True
            if debug:
                print ndp, str(f)
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
                data -= data[-self._start_end[0]-200:-self._start_end[0],...].mean(axis=0)
                #calculate average
                self._datadict[c] = data[:,:,~self._is_artifact[c]].mean(axis=2)
                #print "Shape datadict:", self._datadict[c].shape
                self._is_calculated[c] = True
            else:
                #TODO: some kind of error message
                print "Some problem with the timepoints occured. ndp:", ndp, "max_ts:", max(self._times[c])
                pass 

class PhaseCoherenceAnalyzer(TFBaseAnalyzer):
    """Used for the creation of Phase-Coherence-Data.
       The parameter "average_over" decides wther averaging is 
       performed over trials (standard) or over a given interval 
       of each trial. In the latter case, "average_over" has to be a 
       [s,e] in samples relative to tp 0.
    """
    
    def __init__(self, fn=None, se=None, average_over="trials"):
        TFBaseAnalyzer.__init__(self, fn, se)
        if not average_over=="trials":
            try:
                s,e = average_over
            except Exception,e:
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
                    self._datadict[c] = n.ones((data.shape[0],data.shape[1],data.shape[1]),"d")
                    #calculate phase coherence
                    for ch1 in range(data.shape[1]):
                        for ch2 in range(ch1):
                            pc = phase_coherence(data[:,ch1,~self._is_artifact[c]],data[:,ch2,~self._is_artifact[c]])
                            self._datadict[c][:,ch1,ch2] = pc[:]
                            self._datadict[c][:,ch2,ch1] = pc[:]
                else:
                    s,e = self.average_over
                    self._datadict[c] = n.ones((data.shape[1],data.shape[1],data.shape[2]),"d")
                    for ch1 in range(data.shape[1]):
                        for ch2 in range(ch1):
                            for tr in range(data.shape[2]):
                                if self._is_artifact[c][tr]:
                                    self._datadict[c][ch1,ch2,tr] = n.nan
                                    self._datadict[c][ch1,ch2,tr] = n.nan
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
            
class WaveletAnalyzer(TFBaseAnalyzer):
    """Used for time-frequency-analysis with Morlet-Wavelets.
    Certain measures can be calculated from the wavelets, e.g.
    mean_power: The frequency-resolved power averaged over trials
    itpc: Inter-trial phase coherence"""
#    _wavelets = None
    
    def __init__(self, fn, se=None, freqs=None, Fs=1000.0, measures = None):
        TFBaseAnalyzer.__init__(self, fn, se)
        self._freqs = None
        if freqs == None:
            self._freqs = n.logspace(0.4,2.1,20)
        else:
            self._freqs = freqs
        self._Fs = Fs
        if measures==None:
            self._measures == ["mean_power", "itpc", "st_power","pairwise_pc"]
        else:
            self._measures = []
            for m in measures:
                if m in ["mean_power", "itpc", "st_power","st_power_nb","pairwise_pc"]:
                    if not m in self._measures:
                        self._measures.append(m)
                else:
                    raise ValueError("Unknown measure, %s" %str(m))
        self.st_power_start_end = [0,1000]
        
        #self._wavelets = make_wavelets(self._freqs)
        
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
        
        for c in cnames:
            #assert self._datadict.has_key(c) and self._is_calculated.has_key(c), "Key in dictionary missing!"
            #TODO: open file
            f = eegpy.open_eeg(self._dfile_name)
            #TODO: check timepoints
            ndp = f.numDatapoints
            tps_ok = True
            for t in self._times[c]:
                if t+self._start_end[0]<0 or t+self._start_end[1]>ndp:
                    tps_ok=False
                    break
            if tps_ok:
                self._datadict[c] = {}
                if "mean_power" in self._measures:
                    self._datadict[c]["mean_power"] = n.zeros((self._start_end[1]-self._start_end[0],f.numChannels,len(self._freqs)),"d")
                if "itpc" in self._measures:
                    self._datadict[c]["itpc"] = n.zeros((self._start_end[1]-self._start_end[0],f.numChannels,len(self._freqs)),"d")
                if "st_power" in self._measures:
                    self._datadict[c]["st_power"] = n.zeros((f.numChannels,len(self._times[c]),len(self._freqs)),"d")
                if "st_power_nb" in self._measures:
                    self._datadict[c]["st_power_nb"] = n.zeros((f.numChannels,len(self._times[c]),len(self._freqs)),"d")
                    
                for ch in range(f.numChannels):
                    if debug:
                        print "Condition", c, "Channel", ch
                    #getData
                    data = f.get_data_for_events(self._times[c],self._start_end,[ch])
                    #calculate wavelets
                    wts = n.zeros((data.shape[0],data.shape[2],len(self._freqs)),"D")
                    for j in range(data.shape[2]):
                        wts[:,j,:] = wt_analyze(data[:,0,j],self._freqs, self._Fs)
                        #print wt_analyze(data[:,i,j],self._freqs, self._Fs).shape, self._datadict[c][:,i,j,:].shape
                    if "mean_power" in self._measures:
                        tmp = (abs(wts[:,~self._is_artifact[c],:])**2).mean(axis=1)
                        #print "tmp.shape:", tmp.shape
                        for fr in range(tmp.shape[1]):
                            tmp[:,fr] = 10*( n.log10(tmp[:,fr]) -n.log10 (tmp[-self._start_end[0]-200:-self._start_end[0],fr].mean()))
                        self._datadict[c]["mean_power"][:,ch,:] = tmp
                    if "itpc" in self._measures:
                        tmp = wts[:,:,:] / abs(wts[:,:,:])
                        self._datadict[c]["itpc"][:,ch,:] = abs(tmp.mean(axis=1))
                    if "st_power" in self._measures:
                        tmp = abs(wts[:,:,:])**2
                        for tr in range(tmp.shape[1]):
                            for fr in range(tmp.shape[2]):
                                tmp[:,tr,fr] = 10*( n.log10(tmp[:,tr,fr]) -n.log10 (tmp[-self._start_end[0]-200:-self._start_end[0],tr,fr].mean()))
                        self._datadict[c]["st_power"][ch,:,:] = tmp[-self._start_end[0]+self.st_power_start_end[0]:-self._start_end[0]+self.st_power_start_end[1],:,:].mean(axis=0)
                    if "st_power_nb" in self._measures:
                        tmp = abs(wts[:,:,:])**2
                        #for tr in range(tmp.shape[1]):
                        #    for fr in range(tmp.shape[2]):
                        #        tmp[:,tr,fr] = 10*( n.log10(tmp[:,tr,fr]) -n.log10 (tmp[-self._start_end[0]-200:-self._start_end[0],tr,fr].mean()))
                        self._datadict[c]["st_power_nb"][ch,:,:] = tmp[-self._start_end[0]+self.st_power_start_end[0]:-self._start_end[0]+self.st_power_start_end[1],:,:].mean(axis=0)
                #print "Shape datadict:", self._datadict[c].shape
                self._is_calculated[c] = True
            else:
                #TODO: some kind of error message
                
                pass 
    
#    def mean_power(self, cond, normalize=True):
#        """Calculates the time-frequency-resolved mean power for the condition cond"""
#        rv = abs((self[cond])**2)
#        #Normalisieren
#        if normalize:
#            for i in range(rv.shape[1]):
#                for j in range(rv.shape[2]):
#                    for k in range(rv.shape[3]):
#                        #TODO: more general implementation
#                        rv[:,i,j,k] = rv[:,i,j,k]/rv[(-self._start_end[0])-200:(-self._start_end[0]),i,j,k].mean()
#        rv = rv.mean(axis=2)
#        return rv

        
            
        
    
        
        
        
        
        
