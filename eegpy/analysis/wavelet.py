#!/usr/bin/env python
# -*- coding: utf-8 -*-

from eegpy.helper import fconv, demean, taper, is_power_of_2

import sys
import numpy as n
import pylab as p
from time import sleep, time

try:
    from pywt import Wavelet, WaveletPacket, dwt, dwt2, dwt_coeff_len, dwt_max_level, families, idwt, idwt2, swt, swt_max_level, upcoef, wavedec, wavedec2, wavelist, waverec, waverec2
except ImportError, e:
    print "The pywavelet-package from Filip Wasiliewski is not available on your system. Please check http://www.pybytes.com/pywavelets/."

convolve = fconv #n.convolve #Choose which convolution to use: n.convolve oder fconv

#---- Analyse mit Morlet-Wavelet 
def make_wavelets(freqs,Fs=1000.0,wtSampleWidth=0.8):
    """Creates a set of Morlet-wavelets for analysis"""
    
    assert len(freqs)>0, "freqs must be an iterable object containing the frequencies"
    assert Fs > 0.0, "Fs, the sampling rate, must be larger than zero"
    
    sl=4./(2*n.pi*freqs)
    
    xw = n.arange(-wtSampleWidth,wtSampleWidth,1/Fs)  
    w = n.zeros((len(freqs),xw.shape[0]),"D")
    
    for fi in range(len(freqs)):
        #definition of the Morlet-Wavelets, from Mike
        w[fi,:]=n.exp(2*complex(0,1)*n.pi*freqs[fi]*xw)*n.exp(-xw**2./(2*sl[fi]**2))
    
    return w

def wt_analyze(x,freqs=None,Fs=1000.0,wtSampleWidth=0.8):
    """Analyzes a given dataset via wavelet-analyis. Returns an 2d-Array, with the timeserieses
at the frequencies specified by freqs"""

    assert len(x)>0, "x must be an array containing the data"
    if not freqs==None:
        assert len(freqs)>0, "freqs must be an iterable object with the frequencies"
    else:
        freqs = n.logspace(0.4,2.1,80)
    #tStart = time()
    wts = make_wavelets(freqs,Fs=Fs,wtSampleWidth=wtSampleWidth)
    #print "--", time()-tStart
    x = taper(demean(x)) #preprocessing...
    rv = n.zeros((wts.shape[1]+x.shape[0]-1,len(freqs)), "D")#wts.shape[1]),"D")
    for fi in range(len(freqs)):
        #print convolve(x,wts[fi,:]).shape, rv.shape
        #sys.stdout.flush()
        rv[:,fi] = convolve(x,wts[fi,:])#,"full")
    rv = rv[(wts.shape[1]-1)/2 : rv.shape[0]-1-(wts.shape[1]-1)/2,:]
    return rv
#---- Analyse mit Morlet-Wavelet - Ende

def wt_power_baseline(x,freqs=None,Fs=1000.0,baseline=slice(800,1000),style="dB",wtSampleWidth=0.8):
    """Convenience-function to automatically perform baseline-correction for 
    Morlet-wavelet power-analysis
    style: "dB", "fraction"
      Which output style to use"""
    wt_power = abs(wt_analyze(x,freqs,Fs,wtSampleWidth))**2
    wt_power /= wt_power[baseline,:].mean(axis=0).reshape(1,-1).repeat(wt_power.shape[0],axis=0)
    if style=="dB":
        wt_power = 10*n.log10(wt_power)
    return wt_power

def wavedec_lin(x,wavelet,mode="per",level=None,return_sl=False): 
    """Performs a wavelet-decomposition using pywt.wavedec, but returns the coefficients
    as one 1d-array.
    Uses different default-values for mode and level.
    Behaviour for the case "len(x) not a power of 2" is not complete!"""       
    assert len(x.shape) in [1,2], "Only 1d and 2d-arrays supported up to now!"
    if level == None:
        level = int(round(n.log2(x.shape[0])))
    if len(x.shape) == 1:
        wts = wavedec(x,wavelet,mode=mode,level=level)
        sc_lengths = [len(x) for x in wts]
        rv = n.zeros((sum(sc_lengths)),"d")
        for j in range(len(wts)):
            offset = sum(sc_lengths[:j])
            rv[offset:offset+sc_lengths[j]]=wts[j]
    else: #len(x.shape)==2
        wts = wavedec(x[:,0],wavelet,mode=mode,level=level)
        sc_lengths = [len(i) for i in wts]
        rv = n.zeros((sum(sc_lengths),x.shape[1]),"d")
        for i in range(x.shape[1]):
            if i>0:
                wts = wavedec(x[:,i],wavelet,mode=mode,level=level)
            for j in range(len(wts)):
                offset = sum(sc_lengths[:j])
                rv[offset:offset+sc_lengths[j],i]=wts[j]
    if return_sl:
        return rv, sc_lengths
    else:
        return rv
            
def waverec_lin(x,wavelet,mode="per",sc_lengths=None):
    """Reconstructs the timeseries from the wavelet-coefficients as returned by wavedec_lin.""" 
    assert len(x.shape) in [1,2], "Only 1d and 2d-arrays supported up to now!"
    if sc_lengths != None:
        assert sum(sc_lengths) == x.shape[0], "The values given for sc_lengths don't make sense!"
    else:
        sc_lengths = []
        l = x.shape[0]
        while l>1:
            l/=2
            sc_lengths.insert(0,l)
        sc_lengths.insert(0,l)
        assert sum(sc_lengths) == x.shape[0], "The values automatically calculated for sc_lengths don't make sense!"
            
    if len(x.shape) == 1:
        wts = []
        for j in range(len(sc_lengths)):
            offset = sum(sc_lengths[:j])
            wts.append(x[offset:offset+sc_lengths[j]])
        rv = waverec(wts,wavelet,mode=mode)
    else: #len(x.shape)==2
        for i in range(x.shape[1]):
            wts = []
            for j in range(len(sc_lengths)):
                offset = sum(sc_lengths[:j])
                wts.append(x[offset:offset+sc_lengths[j],i])
            if i==0:
                ar = waverec(wts,wavelet,mode=mode)
                rv = n.zeros((ar.shape[0],x.shape[1]),"d")
                rv[:,i] = ar
            rv[:,i] = waverec(wts,wavelet,mode=mode)
    return rv

def wavepower_lin(x,wavelet,mode="per",level=None,return_sl=False,normalise=True,norm_fraction=None,as_dB=True): 
    """Performs a wavelet-decomposition using pywt.wavedec, but returns the power
    as one 1d-array.
    Uses different default-values for mode and level.

    If normalise==False, len(x) must be a power of two.

    If normalise==True, len(x)*2/3  must be a power of two.
    In this case two Wavelet-transforms are performed, one from 0 to 2/3 
    and one from 1/3 to 1. The center-value of the first wt-power is then used
    for normalization of the second wt-power. The normalization is performed for each 
    scale separately.

    norm_fraction: Not used anymore!
    
    If as_dB, convert normalised power to dB. Ignored if normalise=False
    """       
    assert len(x.shape) in [1,2], "Only 1d and 2d-arrays supported up to now!"
    if normalise:
        assert is_power_of_2(len(x)*2/3), "len(x)*2/3 must be a power of 2"
        if level == None:
            level = int(round(n.log2(x.shape[0]*2/3)))
    else:
        assert is_power_of_2(len(x)), "len(x) must be a power of 2"
        if level == None:
            level = int(round(n.log2(x.shape[0])))
    lx = x.shape[0]
    if len(x.shape) == 1:
        if normalise:
            wts_bl = wavedec(x[:lx*2/3],wavelet,mode=mode,level=level)
            wts = wavedec(x[lx/3:],wavelet,mode=mode,level=level)
        else:
            wts = wavedec(x[:],wavelet,mode=mode,level=level)
        sc_lengths = [len(x) for x in wts]
        rv = n.zeros((sum(sc_lengths)),"d")
        for j in range(len(wts)):
            wts_power = wts[j]**2
            if normalise:
                wts_power /= (wts_bl[j][len(wts_bl[j])/2])**2
                if as_dB:
                    wts_power[:] = 10 * n.log(wts_power)
            offset = sum(sc_lengths[:j])
            rv[offset:offset+sc_lengths[j]]=wts_power[:]
    else: #len(x.shape)==2
        # create array of return-values
        if normalise:
            wts = wavedec(x[lx/3:,0],wavelet,mode=mode,level=level)
        else:
            wts = wavedec(x[:,0],wavelet,mode=mode,level=level)
        sc_lengths = [len(i) for i in wts]
        rv = n.zeros((sum(sc_lengths),x.shape[1]),"d")
        #iterate over dim1
        for i in range(x.shape[1]):
            if i>0:
                # for ch0, wts was already calculated
                if normalise:
                    wts = wavedec(x[lx/3:,i],wavelet,mode=mode,level=level)
                else:
                    wts = wavedec(x[:,i],wavelet,mode=mode,level=level)
            for j in range(len(wts)):
                wts_power = wts[j]**2
                if normalise:
                    wts_bl = wavedec(x[:lx*2/3,i],wavelet,mode=mode,level=level)
                    wts_power /= (wts_bl[j][len(wts_bl[j])/2])**2
                    if as_dB:
                        wts_power[:] = 10 * n.log(wts_power)
                offset = sum(sc_lengths[:j])
                rv[offset:offset+sc_lengths[j],i]=wts_power[:]
    if return_sl:
        return rv, sc_lengths
    else:
        return rv

def wavepower_db_secondhalf(x,wavelet,mode="per",level=None):
    """Unused? Maybe remove it."""
    if level == None:
        level = int(round(n.log2(x.shape[0])))
    if len(x.shape) == 1:
        wts = wavedec(x,wavelet,mode=mode,level=level)
        sc_lengths = [len(x) for x in wts]
        rv = n.zeros((sum(sc_lengths)/2),"d")
        for j in range(len(wts)):
            wts_power = wts[j]**2
            idx_norm = max( 0 , int(wts_power.shape[0]/2)-1 )
                #print "idx_norm:", idx_norm
            wts_power /= wts_power[idx_norm]
            wts_power = 10* np.log(wts_power)
            offset = sum(sc_lengths[:j])/2
            print rv[offset:offset+sc_lengths[j]/2].shape, wts_power[idx_norm+1:].shape
            rv[offset:offset+sc_lengths[j]/2]=wts_power[idx_norm+1:]
    else: #len(x.shape)==2
        return NotImplemented
        wts = wavedec(x[:,0],wavelet,mode=mode,level=level)
        sc_lengths = [len(i) for i in wts]
        rv = n.zeros((sum(sc_lengths),x.shape[1]),"d")
        for i in range(x.shape[1]):
            if i>0:
                wts = wavedec(x[:,i],wavelet,mode=mode,level=level)
            for j in range(len(wts)):
                wts_power = wts[j]**2
                if normalise:
                    idx_norm = max( 0 , int(wts_power.shape[0]*norm_fraction)-1 )
                    wts_power /= wts_power[idx_norm]
                offset = sum(sc_lengths[:j])
                rv[offset:offset+sc_lengths[j],i]=wts_power[:]
    return rv
        
    
#def wavedec_lin_normalise(x,base):
#    """Takes two wavelet-decompositions"""
    

if __name__=='__main__':
    #freqs = n.logspace(0.4,2.1,80)
    #wts = makeWavelets(freqs,wtSampleWidth=0.5)
    #d5 = p.load("/home/thorsten/workspace/Python/src/Spielwiese/data5.txt")
    #print d5.shape
    #res = analyze(d5,Fs=200.0)
    #p.imshow(n.abs(res)**2.,aspect="auto",origin="lower",interpolation="bilinear")
    #p.colorbar()
    #p.show()
    #print freqs
    #print res
    #print res.shape
    #ar = n.random.random((1024+512,10))
    #wtp1 = wavepower_lin(ar[512:],"db4",normalise=False)
    #wtp2 = wavepower_lin(ar,"db4",normalise=True,as_dB=True)
    #wtp3 = wavepower_lin(ar,"db4",normalise=True,as_dB=False)
    import pylab as p
    import eegpy
    from eegpy.plot.wavelet import plot_wavedec_lin
    eeg = eegpy.F32("/media/story/SchlafGed/iEEG/data/canseven_bp.f32")
    et = eegpy.EventTable("/media/story/SchlafGed/iEEG/data/evts/canseven_bp_manuell.evt")
    data = eeg.get_data_for_events(et.events_by_names(*["T %i"%i for i in range(239,256)]),start_end=(-512,1024),channels=[10])
    
    wts = [0,0,0]
    wts[0] = wavepower_lin(data[512:,0,:],"db4",normalise=False)
    wts[1] = wavepower_lin(data[:,0,:],"db4",normalise=True,as_dB=False)
    wts[2] = wavepower_lin(data[:,0,:],"db4",normalise=True,as_dB=True)
    for i in range(3):
        p.subplot(3,1,i+1)
        plot_wavedec_lin(wts[i].mean(axis=1))
        p.colorbar()




    
                      
    
    
