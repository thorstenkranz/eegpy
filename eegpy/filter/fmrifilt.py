#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Mean filtering as used for simultaneous EEG/fMRI"""

from random import shuffle
import os
import time
from progressbar import ProgressBar

import eegpy
from eegpy.misc import FATALERROR
from eegpy.misc import debug
from eegpy.misc import show_progressbar
from eegpy.filter.meanfilt import MeanFilter, _t_before
from eegpy.filter import pcafilt
from eegpy.filter.smoothing import smooth,_smooth, smooth_windowed_eeg, smooth_windowed_eeg_power
from eegpy.filter.freqfilt import filtfilt_high, filtfilt_low, filtfilt_band
from eegpy.helper import find_maxs, find_all_maxs, find_max_overlap, findFirstMax, upsample, upsample_to_memmap, tmp_memmap, downsample, prepend_append_zeros#, ProgressBar
#from eegpy.filter.filt_misc import filterRecursively
if debug:
    import pylab

#################
# Module-Import #
#################

try:
    import numpy as n
    np = n
    from scipy.signal import hilbert
    from scipy.optimize import leastsq, fmin
    from scipy.stats import scoreatpercentile
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

########################
# Function definitions #
########################
test = None

class FMRIFilter(MeanFilter):
    """Implementation of a filter using combined mean/PCA stuff!"""
    _t_slices = None
    
    def __init__(self, x=None, ts=None, len=None, filter=None, t_before=None, num_slices=70, slice_width=86, Fs=1000.0, f_flp=90, ch_for_search=14):
        MeanFilter.__init__(self,x,ts,len,filter,t_before)
        self._num_slices = num_slices
        self._slice_width = slice_width
        self._Fs = Fs      
        self._f_flp = f_flp
        self._ch_for_search = ch_for_search
    
    def set_timepoints_auto(self,start,end=10e50,step=5500,ch_num=0,width=100):
        """If triggers are not good, one can try to automatically find the timepoints.
        For fMRI, search for first volume artifact."""
        #print "Setting timepoints automatically"
        assert ch_num<self._data.shape[1] and ch_num>=0, "ch_num is not valid"
        ts = []
        t = int(start)
        offset=0
        while t<self._data.shape[0]-step and t<end:
            if t==int(start):
                #template = self._data[t-width/2:t+width/2,ch_num]
                searchdata = abs(hilbert(self._data[t:t+500,ch_num]))
                searchdata = smooth(searchdata,self._slice_width/2)
                bp = np.where(searchdata>200)[0][0]
                t = (bp+searchdata[bp:bp+self._slice_width/2].argmax())-self._slice_width/2
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
    
    def find_slices_in_mean(self,num_slices = None,ch_num=None, width=None):
        """In the calculated average data, find all positions of slice-acquisition artifacts relative to trigger.
        Two rounds: first, use a smoothed version (lower risk of failure), 
        second, improve values on unsmoothed signal."""
        if num_slices==None:
            num_slices=self._num_slices
        if width==None:
            width=self._slice_width
        if ch_num==None:
            ch_num=self._ch_for_search
        #First round
        smooth_data = smooth(self._mean_data[:,ch_num]**2,width)
        #global test
        #test=smooth_data
        #smooth_data = (smooth_data**2).mean(axis=1)
        #if debug:
            #print smooth_data.shape, smooth_data
            #pylab.clf()
            #pylab.ioff()
            #pylab.plot(smooth_data[:])
            #pylab.show()
            #raw_input()
        maxs = find_maxs(smooth_data,num_maxs = num_slices, width=width)
        if debug:
            print "Coarse Maxs: ", maxs
        template = self._mean_data[maxs[3]-width/2:maxs[3]+width/2]
        for i in range(len(maxs)):
            try:
                maxs[i] += find_max_overlap(template,self._mean_data[maxs[i]-width/2:maxs[i]+width/2])
            except Exception, e:
                print "Error while refining max-position no. ", i, ", ", e
        if debug:
            if debug:
                print "Finer Maxs: ", maxs
        self._t_slices = maxs
    
    def remove_slice_components(self,width=None):    
        assert self._t_slices != None, "Before removing slice components, first find the slice-positions!"
        if width==None:
            width=self._slice_width
        
        print "Removing slice components"
        if show_progressbar:
            pbar = ProgressBar(self._data.shape[1])
        for cn in range(self._data.shape[1]):
            if show_progressbar:
                pbar.update(cn)
            if debug:
                print "Channel %i", cn
            try:
                for i in self._ts:
                    one_c_slices = n.zeros((width,len(self._t_slices)),"d")
                    for k,j in enumerate(self._t_slices):
                        if debug:
                            print k,i,j
                        one_c_slices[:,k] = self._data[i+j-width/2:i+j+width/2,cn]
                        components = pcafilt.unmix(one_c_slices)
                    for k in range(components.shape[1]):
                            if debug:
                                print "Channel %i, Component %i, std %f, mean(abs) %f" % (cn,k,components[:,k].std(), abs(components[:,k]).mean()),
                            if components[:,k].std()>100.0 or k<2:#acc > 0.03:
                                if debug:
                                    print " removed"
                                components[:,k] = n.zeros(components[:,k].shape,"d")  
                            else:
                                if debug:
                                    print " not removed"
                    #pylab.ioff()
                    #pylab.hist(accs,20)
                    #pylab.show()
                    one_c_slices = pcafilt.mix(components)
                    k=0
                    for k,j in enumerate(self._t_slices):
                        self._data[i+j-width/2:i+j+width/2,cn] = one_c_slices[:,k]
            except Exception,e:
                if debug:
                    print "Error in remove_slice_components,", e
        #pylab.ioff()
        #pylab.hist(accs)
        #pylab.show()
        if show_progressbar:
            pbar.update(self.data.shape[1])
        return self._data
                    
    def remove_common_components(self):
        assert self._mean_data != None, "Before removing components, first calulate the mean!"
        
        print "Removing common components"
        if show_progressbar:
            pbar = ProgressBar(self._data.shape[1])
        for cn in range(self._data.shape[1]):
            #pbar.progress(cn)
            if show_progressbar:
                pbar.update(cn)
            one_c_data = n.zeros((self._len,len(self._ts)),"d")
            for i,t in enumerate(self._ts):
                one_c_data[:,i] = self._data[t:t+self._len,cn]
            components = pcafilt.unmix(one_c_data)
            accs = []
            for i in range(components.shape[1]):
                    acc = abs(n.corrcoef(self._mean_data[:,cn],components[:,i])[0,1])
                    accs.append(acc)
                    if debug:
                        print "Channel %i, Component %i, corrcoef %f, std %f, mean(abs) %f" % (cn,i,acc, components[:,i].std(), abs(components[:,i]).mean()),
                    if (components[:,i].std()>100.0 and i<10) or i<2:#acc > 0.03:
                        if debug:
                            print " removed"
                        components[:,i] = n.zeros(components[:,i].shape,"d")  
                    else:
                        if debug:
                            print " not removed"
            #pylab.ioff()
            #pylab.hist(accs,20)
            #pylab.show()
            one_c_data = pcafilt.mix(components)
            for i,t in enumerate(self._ts):
                try:
                    self._data[t:t+self._len,cn] = one_c_data[:,i]
                except ValueError, e:
                    print "Error at t=%i, cannot write data back after PCA, shapes are %s and %s"%(t,self._data[t:t+self._len,cn].shape,one_c_data[:,i].shape), e
        if debug:
            pylab.ioff()
            pylab.hist(accs)
            pylab.show()
        if show_progressbar:
            pbar.update(self.data.shape[1])
        return self._data
    
    def final_lowpass(self, f=None):
        """Performs simple lowpass, removing sharp edges"""
        if f == None:
            f = self._f_flp
        eeg = self._data
        eeg[:] = filtfilt_low(f,eeg[:],Fs=self._Fs)
    
    def filter(self):
        """Do the filtering"""
        self.check_before_filtering()
        self.calc_mean()
        self.subtract_mean()
        self.find_slices_in_mean()
        self.remove_common_components()
        self.remove_slice_components()
        self.final_lowpass()
        #self.frequency_filter()
        #self.calc_mean()

            
class OBSFilter(FMRIFilter):
    
    upsf = 10 # Upscale-factor
    ch_for_slice_alignment = 14
    _ups_ts = None    
    Ya = None
    Yah = None
    Ya2 = None
    templates = None
    
    def realign_triggers(self):
        if self.Ya == None:
            raise ValueError("OBSFilter: No upsampled data are available - self.Ya not set yet")
        self._ups_ts = [self._ts[0]*self.upsf]
        template = self.Ya[self._ts[0]*self.upsf:(self._ts[0]+self._len)*self.upsf]
        for i,t in enumerate(self._ts[1:]):
            #pylab.clf()
            artifact = self.Ya[t*self.upsf:(t+self._len)*self.upsf]
            shift = find_max_overlap(template, artifact, self.upsf*10)
            #pylab.plot(template[1000:3000])
            #pylab.plot(artifact[1000:3000])
            #pylab.text(40000,500,str(shift))
            #raw_input("Press enter")
            self._ups_ts.append(t*self.upsf+shift)
        print "_ts", len(self._ts), self._ts
        print "_ups_ts", len(self._ups_ts), self._ups_ts
        
        
    def subtract_mean_from_Ya2(self, n_ma=20):
        """Subtract the mean from the upsampled data. 
        A local template for each artifact is created
        """
        def residuals(p,y,x):
            """Residuals used for leastsq-fitting"""
            offset, factor = p
            err = (y*factor-offset)-x
            return abs(err**2).sum()
        
        n_ma = int(n_ma)
        if self.Yah == None or self.Ya == None:
            raise ValueError("OBSFilter: No upsampled data are available - self.Ya not set yet")
        if self.Ya2 == None:
            self.Ya2 = self.Ya.copy()
            #self.Ya2 = n.memmap("tmp/Ya2.np",mode="w+",dtype=n.float64,shape=self.Ya.shape)
        assert n_ma>0 and n_ma<len(self._ts), "Wrong value for n_ma; need 0<n_ma<=len(self._ts)"
        self.templates = n.zeros((self._len*self.upsf,len(self._ups_ts)),n.float32)
        for i,t in enumerate(self._ups_ts): 
            if i%1000==0:
                print i,t,
            try:
                i1 = i-n_ma/2 #Startindex for mean
                i2 = i+n_ma/2 #Endindex for mean
                #Correct these indices
                if i1<0:
                    i1=0
                    i2 = n_ma
                elif i2>len(self._ups_ts):
                    i1 = len(self._ups_ts)-n_ma
                    i2 = len(self._ups_ts)
                for k in range(i1,i2):
                    self.templates[:,i] += self.Yah[self._ups_ts[k]:self._ups_ts[k]+self._len*self.upsf] 
                self.templates[:,i] /= i2-i1
                if i%1000==0:
                    print self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf].shape, abs(self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf]).mean(),
                    print self.templates[:,i].shape, abs(self.templates[:,i]).mean(),
                #TODO: Fit the template to the data
                p_lsq = fmin(residuals,[0,1],args=(self.templates[:,i],self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf]),disp=False)
                offset, factor = p_lsq#[0]
                if i%1000==0:
                    print "Offset:",offset,"Factor:",factor,
                self.templates[:,i]= self.templates[:,i]*factor-offset
                self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf] -= self.templates[:,i]
                if i%1000==0:
                    print self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf].shape, abs(self.Ya2[self._ups_ts[i]:self._ups_ts[i]+self._len*self.upsf]).mean()
            except ValueError,e:
                print "Shape missmatch occured during writing.", t
        return self._data
    
    def subtract_residual_obs(self):
        def obs_fit_error(p,y,x):
            err = n.zeros((x.shape[0]),"d")
            err += p[0]
            for i in range(5):
                err += y[:,i]*p[i+1]
            err -= x
            return abs(err**2).sum()
        assert self.Ya2 != None, "Not yet ready for subtracting residuals. Need to subtract average template first (self.subtract_mean_from_Ya2)"
        self.res_artifacts = n.zeros((self._len*self.upsf,len(self._ups_ts)),n.float64)
        #self.res_artifacts = n.memmap("tmp/res_art.np",mode="w+",dtype=n.float64,shape=(self._len*self.upsf,len(self._ups_ts)))
        #print "Memmapping residual-artifacts"
        for i,t in enumerate(self._ups_ts):
            if i%1000==0:
                print i,
            self.res_artifacts[:,i] = self.Ya2[t:t+self._len*self.upsf]
        #take only 200 artifacts in for obs-calculation 
        idx_ups_ts = range(len(self._ups_ts))
        shuffle(idx_ups_ts)
        bsl = n.zeros((self.res_artifacts.shape[1]),n.bool)
        for i in range(200):
            bsl[idx_ups_ts[i]] = True
        components = pcafilt.unmix(self.res_artifacts[:,bsl])
        for i in range(5):
            pylab.plot(components[:,i],label=str(i))
        pylab.legend()
        self.obs = components[:,:5].copy()
        #Iterate over all residual artifacts and fit OBS to it
        for i,t in enumerate(self._ups_ts):
            p = [0,0,0,0,0,0] 
            p_lsq = fmin(obs_fit_error,p,args=(self.obs,self.res_artifacts[:,i]),maxiter=1e5,maxfun=1e5)#[0]
            print i,t,"p_lsq", p_lsq
            #pylab.clf()
            #pylab.plot(self.res_artifacts[:,i])
            self.fit = n.zeros((self.obs.shape[0]),"d")
            self.fit +=p_lsq[0]
            for j in range(5):
                self.fit += self.obs[:,j]*p_lsq[j+1]
            #pylab.plot(self.fit)
            #pylab.xlim([10000,15000])
            #pylab.savefig("/home/thorsten/tmp/test.png")
            #raw_input()
            print self.res_artifacts[:,i].mean(), self.fit.mean(),
            self.res_artifacts[:,i] -= self.fit
            print self.res_artifacts[:,i].mean()
        for i,t in enumerate(self._ups_ts):
            print i,t,self.Ya2[t:t+self._len*self.upsf].mean(), self.res_artifacts[:,i].mean(),
            self.Ya2[t:t+self._len*self.upsf] = self.res_artifacts[:,i]
            print self.Ya2[t:t+self._len*self.upsf].mean()
    
    def filter(self):
        "Do the actual filtering"
        #Slice-alignment
        self.check_before_filtering()
        self.Ya = upsample(self._data[:,self.ch_for_slice_alignment],self.upsf)
        self.realign_triggers()
        for i in range(self._data.shape[1]):
            self.Ya = upsample(self._data[:,i],self.upsf)
            #ff.Yah = N.memmap("tmp/Yah.np",mode="w+",dtype=N.float64,shape=ff.Ya.shape)
            self.Yah = filtfilt_high(1.0,self.Ya,Fs=10000.0)
            self.subtract_mean_from_Ya2(n_ma=20)
            #ar1 = ff.Ya2.copy()
            self.subtract_residual_obs()
            

class FmriFilterContinuous(object):
    """New approach to filtering, for continuous data only.
    Takes the eeg, finds slices automatically and does a kind of OBS-Filtering."""
    def __init__(self, eeg=None, slice_width=78, k=20, obs_size = 5, f_lp = 70):
        self.eeg = eeg
        self._slice_width = slice_width
        self.ch_for_slice_alignment = 14
        #self.cfsa = self.ch_for_slice_alignment
        self.num_neighbors = k
        self.obs_size = obs_size
        self.f_lp = f_lp
        self.slice_shift=0 #Introduced for analysis, shifting of slice-timepoints
        self._pa_zeros = 5 #Number of zeros to pre-/append to slices
        
    def get_cfsa(self):
        return self.ch_for_slice_alignment
    
    def set_cfsa(self,cfsa):
        self.ch_for_slice_alignment = cfsa
    
    cfsa = property(get_cfsa,set_cfsa)
    
    def find_all_slice_artifacts(self):
        def update_pbar(num):
            """Callback for find_all_maxs"""
            if show_progressbar:
                pbar.update(num/2)
            
        eeg = self.eeg       
        #y = abs(smooth_windowed_eeg(eeg,[self.cfsa],self._slice_width))[:,0]
        #y = smooth_windowed_eeg_power(eeg,[self.cfsa],self._slice_width)[:,0]
        y = filtfilt_band(1,eeg.Fs/self._slice_width,eeg[:,self.cfsa])
        #pylab.plot(y[::10])
        #pylab.plot(eeg[:,14])
        print y.shape, self._slice_width
        #import pylab
        #pylab.ion()
        #pylab.plot(y[0:20000:1])
        #pylab.show()
        #raw_input()
        #pylab.plot(y[13000:20000:1])
        slcs_raw = find_all_maxs(y[:1000],ratio=0.6) # First segment
        slcs_raw.sort()
        #print "t=", t
        #slcs_raw.append(t)
        offset=0
        t=int(0.5*self._slice_width)
        while (t>0.4*self._slice_width or (y.shape[0]-offset)>10000) and (y.shape[0]-offset)>self._slice_width*2: 
            #print (y.shape[0]-offset)
            #print t, offset, "-",
            offset = slcs_raw[-1]+self._slice_width/2
            #print t, offset, "-",
            #pylab.plot(y[offset:offset+self._slice_width])
            #pylab.show()
            #raw_input()
            t=y[offset:offset+self._slice_width].argmax()
            slcs_raw.append(offset+t)   
            #print slcs_raw[-1], slcs_raw[-1]-slcs_raw[-2], " - ",
            #time.sleep(0.1)
            #print t, offset
        print ""
        #pylab.plot(y[::10])
        if show_progressbar:
            pbar = ProgressBar(maxval=eeg.shape[0]/self._slice_width).start()
        #slcs_raw = find_all_maxs(y[:,0],0.3,self._slice_width,20,callback=update_pbar)
        print "Raw slice-positions found", len(slcs_raw), np.mean(np.diff(slcs_raw)), np.min(slcs_raw), np.max(slcs_raw)
        slcs_raw_diff = np.diff(slcs_raw)
        print "slcs_raw_diff: ", scoreatpercentile(slcs_raw_diff,5), scoreatpercentile(slcs_raw_diff,50), scoreatpercentile(slcs_raw_diff,95)
        #raise Exception("Abbruch")
        y , fn = upsample_to_memmap(eeg[:,self.cfsa],10)
        slcs_raw_ups = [x*10 for x in slcs_raw]
        t = slcs_raw_ups[len(slcs_raw)/2]
        template = y[t-self._slice_width*10/2:t+self._slice_width*10/2]
        for i in range(5):
            t = slcs_raw_ups[len(slcs_raw)/2+i]
            template += y[t-self._slice_width*10/2:t+self._slice_width*10/2]
        template /= 6
        offsets = []
        for i,t in enumerate(slcs_raw_ups):
            #offset = find_max_overlap(template, eeg[t-self._slice_width/2:t+self._slice_width/2,self.cfsa], 100)
            offset = find_max_overlap(template, y[t-self._slice_width*10/2:t+self._slice_width*10/2], 100)
            offsets.append(offset)
        self.slcs_ups = [slcs_raw_ups[i]+offsets[i]+self.slice_shift for i in range(len(slcs_raw_ups))]
        if show_progressbar:
            pbar.finish()
        print "Refined slice-positions found. Finished.", len(offsets), np.mean(offsets), np.median(offsets), np.min(offsets), np.max(offsets)
        print "Percentile 0.5,5,95,99.5 of offsets: ", scoreatpercentile(offsets,0.5), scoreatpercentile(offsets,5), scoreatpercentile(offsets,95), scoreatpercentile(offsets,99.5)
        #Adjusting _slice_width...
        print "Old slice_width:", self._slice_width
        self._new_slice_width = int(n.ceil(n.mean(n.diff(self.slcs_ups))))/10
        self._new_slice_width += 3 # Make slice wider to have no zombie-timepoints
        self._new_slice_width = self._new_slice_width+self._new_slice_width%2
        #self._new_slice_width = (self._new_slice_width/2)*2 # make sw%2==0 (divisible by 2)
        print "New slice_width:", self._new_slice_width
        #raise Exception("Abbruch")
        return [x/10 for x in self.slcs_ups]
        
    def remove_artifacts(self):
        """Tries to remove the artifacts. 
        Uses mixture of mean-subtraction and OBS-PCA-Subtraction
        Needs self.slcs to be set.
        Saves example data in self.examples to make evaluation of filter quality possible"""
        def obs_fit_error(p,y,x):
            err = n.zeros((x.shape[0]),"d")
            err += p[0]
            for i in range(self.obs_size):
                err += y[:,i]*p[i+1]
            err -= x
            return abs(err**2).sum()
        
        def obs_fit_error_lsq(p,y,x):
            err = n.zeros((x.shape[0]),"d")
            err += p[0]
            for i in range(self.obs_size):
                err += y[:,i]*p[i+1]
            err -= x
            return err
        
        def las_fit_error_lsq(p,y,x):
            """Local artifact subtraction"""
            err = p[0]*y
            err -= x
            return err
        
        #Shortnames
        eeg = self.eeg
        slcs = self.slcs_ups
        sw = self._new_slice_width *10
        k = self.num_neighbors
        self.examples=[]
        
        num_examples=10
        
        #Loop over channels
        if show_progressbar:
            pbar = ProgressBar().start()
        for i_ch in range(eeg.num_channels):
            #Make Highpass-Version of channel
            #ch = eeg[:,i_ch].copy()
            y , fn = upsample_to_memmap(eeg[:,i_ch],10)
            y_hp,fn_hp = tmp_memmap(dtype=y.dtype,shape=y.shape,mode="w+")
            y_hp[:] = filtfilt_high(1, y, Fs=10000.0)
            y_out,fn_out = tmp_memmap(dtype=y.dtype,shape=y.shape,mode="w+")
            y_out[:] = y[:]
            #ch_hp = filtfilt_high(1.0,ch,Fs=eeg.Fs)
            neighbors = n.zeros((sw+2*self._pa_zeros,k))
            #Prefill
            for i in range(k):
                #print neighbors[:,i].shape, eeg[slcs[i]-sw/2:slcs[i]+sw/2,i_ch].shape
                neighbors[:,i] = prepend_append_zeros( y_hp[slcs[i]-sw/2:slcs[i]+sw/2] , self._pa_zeros)
            #Loop over slices and filter
            next_subst = 0
            for i,t in enumerate(slcs):
                try:
                    if i>k/2 and i<len(slcs)-k/2:
                        neighbors[:,next_subst] = prepend_append_zeros( y_hp[slcs[i+k/2]-sw/2:slcs[i+k/2]+sw/2] , self._pa_zeros)
                        next_subst+=1
                        next_subst=next_subst%k
                    tmp = prepend_append_zeros( y[t-sw/2:t+sw/2] , self._pa_zeros)
                    #Subtraction
                    #Shift/scale template
                    template = neighbors.mean(axis=1)
                    p = [1.0, 0]
                    p_las = leastsq(las_fit_error_lsq,p,args=(template,tmp))[0]
                    #print p_las[0],
                    tmp -= p_las[0]*template
                    
                    #Beispieldaten speichern: Teil 1
                    if i_ch == self.ch_for_slice_alignment:
                        if i%(len(slcs)/num_examples)==(len(slcs)/num_examples)/2:
                            print "examples, Teil 1"
                            example = {}
                            example["raw"] = prepend_append_zeros(y[t-sw/2:t+sw/2].copy() , self._pa_zeros )
                            example["mean"] = p_las[0]*template
                    #OBS-Fit
                    components = pcafilt.unmix(neighbors) #OBS will be first 5 components
                    components -= components.mean(axis=0).reshape(1,-1).repeat(components.shape[0],0) #demeaning column-wise
                    obs = components[:,:self.obs_size].copy()
                    #Fit OBS to artifact
                    p = [0]+[0]*self.obs_size
                    #p_lsq = fmin(obs_fit_error,p,args=(obs,eeg[t-sw/2:t+sw/2,i_ch]),maxiter=1e5,maxfun=1e5)#[0]
                    p_lsq = leastsq(obs_fit_error_lsq,p,args=(obs,tmp))[0]
                    #print i,t,"p_lsq", p_lsq
                    fit = n.zeros((obs.shape[0]),"d")
                    fit +=p_lsq[0]
                    for j in range(self.obs_size):
                        fit += obs[:,j]*p_lsq[j+1]
                    tmp -= fit
                    try:
                        #eeg[t/10-sw/10/2:t/10+sw/10/2,i_ch] = tmp[self._pa_zeros:-self._pa_zeros][::10]
                        y_out[t-sw/2:t+sw/2] = tmp[self._pa_zeros:-self._pa_zeros][:]
                    except ValueError, ve:
                        print i_ch, i,t, eeg[t/10-sw/10/2:t/10+sw/10/2,i_ch].shape, tmp[self._pa_zeros:-self._pa_zeros][::10].shape
                    
                    #Beispieldaten speichern: Teil 2
                    if i_ch == self.ch_for_slice_alignment:
                        if i%(len(slcs)/num_examples)==(len(slcs)/num_examples)/2:
                            print "examples, Teil 2"
                            #example["fit"] = n.zeros(example["raw"].shape) #fit.copy()
                            example["fit"] = fit.copy()
                            example["obs"] = obs.copy()
                            example["filt1"] = (tmp + fit).copy()
                            example["filt2"] = tmp.copy()
                            self.examples.append(example)
                
                    if show_progressbar:
                        pbar.update((i_ch+i/len(slcs))*100/eeg.num_channels)
                except Exception, e:
                    print "Error occurred at slice at t=",t,", ignored"
                    print e
            y_out[:] = filtfilt_low(self.f_lp,y_out[:],Fs=eeg.Fs*10)
            eeg[:,i_ch] = y_out[::10]
                    
        if show_progressbar:
            pbar.finish()
    
    def final_lowpass(self, f=90):
        """Performs simple lowpass, removing sharp edges"""
        eeg = self.eeg
        eeg[:] = filtfilt_low(f,eeg[:],Fs=eeg.Fs)
        
                        
                
                
            
           
       
#########################
# Convenience-functions #
#########################        
_fmrifilter = None #Global Instance of MeanFilter for use with methods below

def filter(x,ts,len=None,filter=None,t_before=None):
    """Do the fmri-filter, with array x."""
    global _filter
    _fmrifilter = FMRIFilter(x,ts,len,filter,t_before)
    return _fmrifilter.filter()


        
#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':
    
    #TODO: Make some example use
    print "Example code not implemented yet."
    