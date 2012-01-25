#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mean filtering as used for simultaneous EEG/fMRI"""

import eegpy
from eegpy.misc import FATALERROR
from eegpy.misc import debug
from eegpy.helper import upsample,downsample,find_max_overlap
#from eegpy.filter.filt_misc import filterRecursively
#if debug:
#    import pylab

#################
# Module-Import #
#################

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

# Global variables 
_t_before = 50 #Time to reduce triggers by

########################
# Function definitions #
########################

class MeanFilter():
    """Implementation for the filtering by subtracting the mean.
    Should work with EEG-files and with arrays, first only with arrays.
    Works in-place!!!"""
    
    _data=None
    _ts = None
    _len = None
    _mean_data = None
    _filter = None#lambda x:x
    _corrts = None
    _t_before = _t_before
    
    def __init__(self, x=None, ts=None, len=None, filter=None, t_before=None):
        if x!=None:
            self.set_data(x)
        if t_before != None:
            self._t_before = t_before
            #print "t_before in meanfilt:"
        if ts!=None:
            self.set_timepoints(ts)
        if len!=None:
            self.set_length(len)
        if filter != None:
            self._filter = filter
        
    def get_data(self):
        return self._data    
                
    def set_data(self, x):
        try:
            a = x.shape
        except Exception, e:
            raise ValueError("Only arrays / memmaps can be used as input for MeanFilter")
        self._data = x
    
    def set_timepoints(self,ts):
        self._ts = [int(x)-self._t_before for x in ts]
        self._ts.sort()
    
    def get_length(self):
        return self._len
    
    def set_length(self,l):
        self._len = int(l)
    
    def set_timepoints_auto(self,start,end=10e50,step=5500,ch_num=0,width=100):
        """If triggers are not good, one can try to automatically find the timepoints"""
        #print "Setting timepoints automatically"
        assert ch_num<self._data.shape[1] and ch_num>=0, "ch_num is not valid"
        ts = []
        t = int(start)
        offset=0
        while t<self._data.shape[0]-step and t<end:
            if t==int(start):
                #template = self._data[t-width/2:t+width/2,ch_num]
                template = self._data[t:t+step,ch_num]
                ts.append(t)
            else:
                #offset = find_max_overlap(template, self._data[t-width/2:t+width/2,ch_num], width/2)
                offset = find_max_overlap(template, self._data[t:t+step,ch_num], width/2)
                #print offset
                ts.append(t+offset)
            if debug:
                print ts[-1],
            t+=step+offset
        self.set_timepoints(ts)
        return ts
                
            
    
    def check_before_filtering(self):
        if self._data == None:
            raise RuntimeError("No data were set yet.")
        if self._ts == None:
            raise RuntimeError("No timepoints were set yet.")
        if self._len == None:
            #Auto calculate _len as the minimal distance between subsequent timepoints
            self._len = n.diff(n.array(self._ts)).min()
            if self._len/n.diff(n.array(self._ts)).mean()<0.1:
                print "Auto-setting of length might have failed..., self._len = ", self._len
        #Now all three values needed should be o.k.
        if debug:
            print "self._data.shape:", self._data.shape
            print "self._ts:", self._ts, n.diff(n.array(self._ts)).mean()
            print "self._len:", self._len
            
    def calc_mean(self):
        #self._corrts = n.zeros((len(self._ts)),"i")
        mean_shape = []
        mean_shape.append(self._len)
        for i in range(1,len(self._data.shape)):
            mean_shape.append(self._data.shape[i])
        #print "mean_shape = ", mean_shape
        #print "Sollte sein: ",self._data[self._ts[0]:self._ts[0]+self._len,...].shape
        self._mean_data = n.zeros(mean_shape,"d")
        #print self._mean_data.shape, self._data[self._ts[0]:self._ts[0]+self._len,...].shape
        for i, t in enumerate(self._ts): #First walk-through, calculating average
            #self._mean_data += self._filter(self._data[t:t+self._len,...])
            try:
                d = self._data[t:t+self._len,...]
                #d_up = upsample(d)
                self._mean_data += d
                #if i==0:
                    #self._mean_data += d#_up
                #else:
                    #self._corrts[i] = find_max_overlap(self._mean_data[:,10],d_up[:,10])
                    #if self._corrts[i]<0:
                    #    self._mean_data[:i,...] += d_up[-i:,...]
                    #elif self._corrts[i]==0:
                    #    self._mean_data += d_up
                    #else:
                    #    self._mean_data[i:,...] += d_up[:-i,...]
                #print "maxOverlap:", self._corrts[i]
                #pylab.clf()
                #pylab.ioff()
                #for i in range(d.shape[1]):
                #    pylab.subplot(8,d.shape[1]/8,i+1)
                #    pylab.plot(d[:,i])
                #pylab.show()
            except ValueError,e:
                print self._mean_data.shape, self._data[t:t+self._len,...].shape#, d_up.shape
                print "Shape missmatch occured. Seems like the last timepoint isn't completely recorded in the eeg and is therefore ignored for template-creation."
                if debug:
                     print e
        self._mean_data /= len(self._ts)
        #self._mean_data = downsample(self._mean_data)
        #if debug:
            #pylab.clf()
            #pylab.ioff()
            #for i in range(self._mean_data.shape[1]):
            #    pylab.subplot(8,self._mean_data.shape[1]/8,i+1)
            #    pylab.plot(self._mean_data[:,i])
            #pylab.show()
                
    def subtract_mean(self, n_ma=None):
        """Subtract the mean artifact from all artifact occurences
        n_ma: Number of nearby artifacts to include in average. 
          If None (default), do global subtraction
        """
        if n_ma==None:
            assert self._mean_data != None, "Before subtracting, first calulate the mean!"
            for i,t in enumerate(self._ts): #Second walk: Substracting mean
                try:
                    #self._data[t:t+self._len,...] -= downsample(self._mean_data,start=self._corrts[i])
                    self._data[t:t+self._len,...] -= self._mean_data
                except ValueError,e:
                    print "Shape missmatch occured during writing.", t
        else:
            n_ma = int(n_ma)
            assert n_ma>0 and n_ma<len(self._ts), "Wrong value for n_ma; need 0<n_ma<=len(timepoints)"
            for i,t in enumerate(self._ts): #Second walk: Substracting mean
                try:
                    i1 = i-n_ma/2 #Startindex for mean
                    i2 = i+n_ma/2 #Endindex for mean
                    #Correct these indices
                    if i1<0:
                        i1=0
                        i2 = n_ma
                    elif i2>len(self._ts):
                        i1 = len(self._ts)-n_ma
                        i2 = len(self._ts)
                    self._data[t:t+self._len,...] -= self._mean_data
                except ValueError,e:
                    print "Shape missmatch occured during writing.", t
        return self._data
    
    def frequency_filter(self):
        """Applies the filter that is given to the constructor using the keyword filter.
        If no filter is given, pass.
        """
        #TODO: Change behaviour. Filter each channel as a whole and use lowpass by default.
        if self._filter == None:
            pass
        else:
            for i in range(0,self._data.shape[0],10000):
                try:
                    self._data[i:i+10000,:] = self._filter(self._data[i:i+10000,:])
                except Exception, e:
                    if debug:
                        print "Error in MeanFilter.frequency_filter:", e 
        
    
    def filter(self):
        """Do the filtering"""
        self.check_before_filtering()
        self.calc_mean()
        self.subtract_mean()
        #self.frequency_filter()
        
    data = property(get_data,set_data)
    length = property(get_length,set_length)
            

#########################
# Convenience-functions #
#########################        
_meanfilter = None #Global Instance of MeanFilter for use with methods below

def filter(x,ts,len=None,filter=None,t_before=None):
    """Do the PCA, with array x."""
    global _meanfilter
    _meanfilter = MeanFilter(x,ts,len,filter,t_before)
    return _meanfilter.filter()


        
#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':

    #TODO: Make some example use
    print "Example code not implemented yet."
    