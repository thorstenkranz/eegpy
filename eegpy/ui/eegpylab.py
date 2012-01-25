#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Collection of high-level functions for interactive use"""

import os.path

import eegpy
from eegpy.filter.freqfilt import filtfilt_low, filtfilt_high,filtfilt_band,filtfilt_bandstop


def freqfilt_eeg(fn_in,fn_out,btype="lp",fl=None,fh=None,border=3,windowed=False):
    eeg = eegpy.open_eeg(fn_in)
    out = eegpy.F32(fn_out,"w+",shape=eeg.shape,cNames=eeg.channel_names,Fs=eeg.Fs)
    if btype == "lp":
        if not windowed:    
            out[:,:] = filtfilt_low(fl,eeg[:,:],Fs=eeg.Fs,border=border)
        else:
            for i in range(0,out.num_datapoints,100000):
                out[i:min(i+100000,out.num_datapoints),:] = filtfilt_low(fl,eeg[i:min(i+100000,out.num_datapoints),:],Fs=eeg.Fs,border=border)
    elif btype == "hp":
        #for i_c in range(eeg.num_channels):
        if not windowed:    
            out[:,:] = filtfilt_high(fh,eeg[:,:],Fs=eeg.Fs,border=border)
        else:
            for i in range(0,out.num_datapoints,100000):
                out[i:min(i+100000,out.num_datapoints),:] = filtfilt_high(fh,eeg[i:min(i+100000,out.num_datapoints),:],Fs=eeg.Fs,border=border)
    elif btype == "bp":
        if not windowed:    
            out[:,:] = filtfilt_band(fl,fh,eeg[:,:],Fs=eeg.Fs,border=border)
        else:
            for i in range(0,out.num_datapoints,100000):
                out[i:min(i+100000,out.num_datapoints),:] = filtfilt_band(fl,fh,eeg[i:min(i+100000,out.num_datapoints),:],Fs=eeg.Fs,border=border)
    elif btype == "bs":
        if not windowed:    
            out[:,:] = filtfilt_bandstop(fl,fh,eeg[:,:],Fs=eeg.Fs,border=border)
        else:
            for i in range(0,out.num_datapoints,100000):
                out[i:min(i+100000,out.num_datapoints),:] = filtfilt_bandstop(fl,fh,eeg[i:min(i+100000,out.num_datapoints),:],Fs=eeg.Fs,border=border)
        
def extract_part_from_eeg(fn_in,fn_out,start,end,et_in=None,et_out=None):
    """Extract a part of an EEG, write it to new file"""
    eeg = eegpy.open_eeg(fn_in)
    out = eegpy.F32(fn_out,"w+",shape=(end-start,eeg.shape[1]),cNames=eeg.channel_names,Fs=eeg.Fs)
    out[:] = eeg[start:end,:]
    if not et_in==None:
        if et_out==None:
            et_out=os.path.join(os.path.splitext(fn_out)[0],".evt")
        et = eegpy.EventTable(et_in)
        et-=start
        et.save(et_out)


