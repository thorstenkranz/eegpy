#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR, debug
from eegpy.analysis.phasespace import find_n_nearest_neighbors,find_n_nearest_neighbors_fast
from eegpy.analysis.wavelet import wavedec, wavedec_lin, waverec, waverec_lin
from freqfilt import filtfilt, filtfilt_high
import pcafilt
#from eegpy.filter.filt_misc import filterRecursively

#################
# Module-Import #
#################

try:
    import numpy as n
    np=n
    import scipy
    from scipy.signal import butter
    from scipy.optimize import leastsq, fmin
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    
########################
# Function definitions #
########################

def _highpass(x,Fs=1000.):
    return filtfilt_high(2.,x,Fs=Fs)

def bcg_filt_wavelet(x,ts,m=1024,wave_name="coif4",k=20, lambda_r=1,Fs=1000.0):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI"""
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]*window
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
        
        if debug:
            print "Find nearest neighbors"
        neighbors = find_n_nearest_neighbors(bcgs.T,k)
        
#        if debug:
#            print "Testing neighborhood: Do some plots..."
#        for i in range(0,neighbors.shape[0],50):
#            boolslice = n.zeros((bcgs.shape[1]),n.bool)
#            #neighborhood = n.zeros((m,k),"d")
#            for j in range(k):
#                boolslice[neighbors[i,j]] = True
#            neighborhood = bcgs[:,boolslice]
            
            #p.plot(neighborhood.mean(axis=1))
            #p.plot(bcgs[:,:].mean(axis=1))
            
        if debug:
            print "Wavelet transform"
        wavelet_tr_s_lin = wavedec_lin(bcgs,wave_name)
        #walk through all points
        #Calculate the center-of-mass and sd for all wavelets
        sqrt_kp1 = n.sqrt(k+1)
        
        wts_to_shrink = wavelet_tr_s_lin.copy()
        #print "Shrinking for lambda_r =", lambda_r 
        for i in range(wts_to_shrink.shape[1]):
            #numShrinked.append(0)
            if debug:
                if i%100==0:
                    if debug:
                        print "i:", i
            boolslice = n.zeros((wts_to_shrink.shape[1]),n.bool)
            for ni,neighbor in enumerate(neighbors[i]):
                boolslice[neighbor] = True
            NH = wavelet_tr_s_lin[:,boolslice]            
            COM = NH.mean(axis=1)
            STD = NH.std(axis=1)
            #apply the shrinking condition
            for idx in range(wts_to_shrink.shape[0]):
                if abs(COM[idx])<(2*lambda_r*STD[idx]/sqrt_kp1):
                    wts_to_shrink[idx,i]= COM[idx]
                    
        if debug:
            print "Wavelet-reconstruction"
        bcgs_f = waverec_lin(wts_to_shrink,wave_name)
        
        if debug:
            print "Writing back results to eeg"
        #I have to subtract the filtered bcgs from the eeg as these are now the 
        #"pure" artifact.
        for i in range(bcgs.shape[1]):
            try:
                x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] -= bcgs_f[:,i]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e

def bcg_filt_average_nh(x,ts,m=1024,k=20,Fs=1000.0):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI using AAS 
    (average artifact subtraction), but in local neighborhood"""
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
        
        if debug:
            print "Find nearest neighbors"
        neighbors = find_n_nearest_neighbors(bcgs.T,k)
        
#        if debug:
#            print "Testing neighborhood: Do some plots..."
#        for i in range(0,neighbors.shape[0],50):
#            boolslice = n.zeros((bcgs.shape[1]),n.bool)
#            #neighborhood = n.zeros((m,k),"d")
#            for j in range(k):
#                boolslice[neighbors[i,j]] = True
#            neighborhood = bcgs[:,boolslice]
            
            #p.plot(neighborhood.mean(axis=1))
            #p.plot(bcgs[:,:].mean(axis=1))


        for i in range(bcgs.shape[1]):
            boolslice = n.zeros((bcgs.shape[1]),n.bool)
            for ni,neighbor in enumerate(neighbors[i]):
                boolslice[neighbor] = True
            try:
                x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] -= bcgs[:,boolslice].mean(axis=1) * window
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
                
def bcg_filt_average_nh_fit(x,ts,m=1024,k=20,Fs=1000.0):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI using AAS 
    (average artifact subtraction), but in local neighborhood. Additionally perform a fit
    of the template to each artifact."""
    def residuals(p,y,x):
        """Residuals used for leastsq-fitting"""
        factor = p[0]
        #err = (y*factor-offset)-x
        err = (y*factor)-x
        return abs(err**2).sum()
    
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x_hp = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x_hp[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
        
        if debug:
            print "Find nearest neighbors"
        neighbors = find_n_nearest_neighbors(bcgs.T,k)
        
        for i in range(bcgs.shape[1]):
            boolslice = n.zeros((bcgs.shape[1]),n.bool)
            for ni,neighbor in enumerate(neighbors[i]):
                boolslice[neighbor] = True
            try:
                template = bcgs[:,boolslice].mean(axis=1) * window
                p_lsq = fmin(residuals,[1],args=(template,x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]),disp=False)
                factor = p_lsq[0]
                #print i, p_lsq
                if factor>0:
                    #template = template*factor-offset
                    template = template*factor
                    x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] -= template
            except ValueError,e:
                print "Value Error, probably shape missmatch", e

def bcg_filt_average(x,ts,m=1024,Fs=1000.0):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI using AAS 
    (average artifact subtraction), but in local neighborhood"""
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
        
        for i in range(bcgs.shape[1]):
            try:
                x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] -= bcgs[:,:].mean(axis=1) * window
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
                
def bcg_filt_pca_nh(x,ts,m=1024,k=20):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI using PCA
    (principal component analysis), but in local neighborhood"""
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
        
        if debug:
            print "Find nearest neighbors"
        neighbors = find_n_nearest_neighbors(bcgs.T,k)

        for i in range(bcgs.shape[1]):
            boolslice = n.zeros((bcgs.shape[1]),n.bool)
            for ni,neighbor in enumerate(neighbors[i]):
                boolslice[neighbor] = True
            NH = bcgs[:,boolslice]
            components = pcafilt.unmix(NH)
            for k in range(components.shape[1]):
                if components[:,k].std()>100.0:#acc > 0.03:
                    if debug:
                        print "Channel %i, Component %i, std %f, mean(abs) %f" % (ch_num,k,components[:,k].std(), abs(components[:,k]).mean()), "removed"
                    components[:,k] = n.zeros(components[:,k].shape,"d")  
                else:
                    if k==0:
                        if debug:
                            print "Channel %i, Component %i, std %f, mean(abs) %f" % (ch_num,k,components[:,k].std(), abs(components[:,k]).mean()), "not removed"
            NH = pcafilt.mix(components)
            try:
                x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] = NH[:,0]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e  
                
def bcg_filt_pca(x,ts,m=1024,Fs=1000.0):
    """Filter BCG-artifacts from simultaneous recordings of EEG and fMRI using PCA
    (principal component analysis)"""
    assert len(x.shape) in [1,2], "x must be a 1d- or 2d-array"
    if len(x.shape) == 1:
        x=x.reshape((x.shape[0],1))
    x = _highpass(x,Fs=Fs)
    window = n.hanning(m)
    #Other possibilities for window:
    #n.ones((m),"d")
    #n.hamming(m)
    #n.bartlett(m)
    ts_bcg = ts

    for ch_num in range(x.shape[1]):
        if debug:
            print "Get all BCG-artifacts of channel %i in one array" % ch_num
        bcgs = n.zeros((m,len(ts_bcg)),"d")
        for i in range(bcgs.shape[1]):
            #print eeg[ts_bcg[i]-m/2:ts_bcg[i]+m/2,14].shape
            try:
                bcgs[:,i] = x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e
                
        components = pcafilt.unmix(bcgs)
        for k in range(components.shape[1]):
            if components[:,k].std()>100.0:#acc > 0.03:
                if debug:
                    print "Channel %i, Component %i, std %f, mean(abs) %f" % (ch_num,k,components[:,k].std(), abs(components[:,k]).mean()), "removed"
                components[:,k] = n.zeros(components[:,k].shape,"d")  
            else:
                if k==0:
                    if debug:
                        print "Channel %i, Component %i, std %f, mean(abs) %f" % (ch_num,k,components[:,k].std(), abs(components[:,k]).mean()), "not removed"
        bcgs = pcafilt.mix(components)
        for i in range(bcgs.shape[1]):
            try:
                x[ts_bcg[i]-m/2:ts_bcg[i]+m/2,ch_num] = bcgs[:,i]
            except ValueError,e:
                print "Value Error, probably shape missmatch", e    
                

class BCGFinder(object):
    """Find the timepoints where the BCG-artifacts are. 
    """
    rv_ts = None
    
    def __init__(self,x,ts, Fs=1000.0, width=300, del_width=700):                              
        """x: array-like to search in. Could be F32-object.
        ts: sample-timepoints. Used to define a pattern for the search.
        """
        assert len(x.shape) == 2, "x must be a 2d-arraylike object"
        assert len(ts) > 0, "ts must be a list of timepoints"
        #self._x is set below
        self._ts = [int(t) for t in ts]
        self._ts.sort()
        self._pulse_init = int(np.median(np.diff(self._ts)))
        self._width = width
        self._del_width = del_width
        #Setup Bandpass-Filter
        self._Fs = Fs
        [b,a]=butter(3,[1./(Fs/2),13./(Fs/2)], btype="band")
        self.bp_filt = lambda x: filtfilt(b,a,x)
        self._x = self.bp_filt(x[:])
    
    def find_all_bcg(self,pulse_cutoff=1500):
        """Searches for all BCG-artifacts in array"""
        #create ProgressBar
        positions = []
        #Calculate averages
        self._mean = n.zeros((self._width,self._x.shape[1]),"f")
        for t in self._ts:
            self._mean += self._x[t-self._width/2:t+self._width/2,:]
        self._mean /= len(self._ts)
        corr = self.correlate2d(self._x[:], self._mean)
        #print "corr", corr.shape
        pulse = self._pulse_init
        cutoff = 0.4
        alocs=[self._ts[0]]
        #Go back from first t
        t_tmp = self._ts[0]-pulse
        print pulse
        while t_tmp>pulse:
            pos_maxc = corr[t_tmp-pulse/3:t_tmp+pulse/3].argmax()
            aloc = pos_maxc+(t_tmp-pulse/3) #new Position
            if corr[aloc]>cutoff:
                alocs.insert(0,aloc)
            t_tmp=aloc-pulse
            print aloc
        #Go forward from first t
        t_tmp = self._ts[0]+pulse
        while t_tmp<corr.shape[0]-pulse:
            aloc = corr[t_tmp-pulse/3:t_tmp+pulse/3].argmax()+(t_tmp-pulse/3) #new Position
            if corr[aloc]>cutoff:
                alocs.insert(0,aloc)
            t_tmp=aloc+pulse
            print aloc
        self.rv_ts = alocs
        self.rv_ts.sort()
        return self.rv_ts
        
    def find_all_bcg_old(self,pulse_cutoff=1500):
        """Searches for all BCG-artifacts in array"""
        #create ProgressBar
        positions = []
        #Calculate averages
        self._mean = n.zeros((self._width,self._x.shape[1]),"f")
        for t in self._ts:
            self._mean += self._x[t-self._width/2:t+self._width/2,:]
        self._mean /= len(self._ts)
        corr = self.correlate2d(self._x[:], self._mean)
        #print "corr", corr.shape
        pulse = 1e50
        cutoff = 0.4
        while pulse>pulse_cutoff:
            try: #Finde Artefakte
                cutoff*=0.75
                corrC = corr.copy()
                #self.a.clear()
                #self.a.plot(corrC[::10])
                alocs = []
                allMax = 0
                i=0
                #First round
                while True:#for i in range():
                    aloc = corrC.argmax()
                    if i == 0:
                        allMax = corrC[aloc]
                    #print "corrC[aloc]/allMax", corrC[aloc]/allMax
                    if corrC[aloc]/allMax>cutoff:
                        alocs.append(aloc)
                    else:
                        break
                    #Bereich um lokales Maximum auslöschen
                    dw = self._del_width/2
                    startDel = aloc-dw
                    endDel = aloc+dw
                    if startDel<0:
                        startDel=0
                    if endDel>corrC.shape[0]:
                        endDel=corrC.shape[0]
                    corrC[startDel:endDel] = n.zeros((endDel-startDel),"f")
                    i+=1
                alocs.sort()
            except Exception, e:
                print "Fehler beim finden der Artefakte:", e
            pulse = n.median(n.diff(n.array(alocs)))
            print "Pulse for cutoff=%.2f: %i" %(cutoff,int(pulse))
        self.rv_ts = [x for x in alocs]
        #Second round: go through the list of timepoints, find additional timepoints
        #2.1: Go from first tp to beginning
        alocs = []
        t = self.rv_ts[0]
        while True:
            t -= pulse
            if t<pulse/2:
                break #at the beginning of the 
            t_c = corr[t-pulse/2:t+pulse/2].argmax()
            if t_c==0 or t_c > pulse-2:
                continue #Maximum at beginning/end of segment -> No real max found
            alocs.append(int(t-pulse/2+t_c))
        print "Tps added before first tp:", alocs
        #2.1: Go from each tp to the next one and check if anotherone should be inserted
        for i,tp in enumerate(self.rv_ts):
            t = tp
            while True:
                t+= pulse
                if i<len(self.rv_ts)-1:
                    if t>self.rv_ts[i+1]-pulse/2:
                        break #too close to next tp
                else:
                    if t>self._x.shape[0]-pulse/2:
                        break #too close to end of _x
                t_c = corr[t-pulse/2:t+pulse/2].argmax()
                if t_c==0 or t_c > pulse-2:
                    continue #Maximum at beginning/end of segment -> No real max found
                alocs.append(int(t-pulse/2+t_c))        
        print "Tps added in correction:", alocs
        #Append correction-tps to rv_ts
        for t in alocs:
            self.rv_ts.append(t)
        self.rv_ts.sort()
        return self.rv_ts
        
    def correlate2d(self, x,y):
        assert len(x.shape)==2, "Array x must be 2d"
        assert len(y.shape)==2, "Array y must be 2d"
        #print x.shape, y.shape
        assert x.shape[1] == y.shape[1], "Number of channels must be equal in both arrays."
        
        #Bandpass on array y (x is self.data, already filtered)
        y = self.bp_filt(y)
        #weights = n.array([0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        #weights = n.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
        #weights = n.zeros((x.shape[1]),"d")
        #for i in self.plot.chList:
        #    weights[i] = 1
        #print weights
        #weights = weights / (weights.mean()*weights.shape[0])
        corr = n.correlate(x[:,0],y[:,0],"same")
        for i in range(1,x.shape[1]):
            corr += n.correlate(x[:,i],y[:,i],"same")
        corr/=x.shape[1]
        #corr /= x.shape[1]
        #print corr.shape
        return corr
    
#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':

    print "This is the bcgfilt-module. No examples implemented yet"
