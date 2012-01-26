#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Model for EEG-BOLD-Coupling via EEG-Power"""

import numpy as np
from scipy.interpolate import splrep,splev
from scipy.signal import hilbert
from scipy.integrate import odeint
from eegpy.filter.freqfilt import filtfilt_band
from eegpy.filter.smoothing import smooth
from eegpy.models.eegbold import Hrf2d
import cPickle
import nose

import pylab as p
from eegpy.analysis.wavelet import wt_power_baseline
from scipy.stats import scoreatpercentile

int_range = ir =[0,2000]
nnodes=100
positions = np.zeros((nnodes,2))
for i in range(nnodes):
    x=i/10
    y=i%10
    positions[i,0] = x*1
    positions[i,1] = y*1
hrfm = Hrf2d(n_nodes=nnodes,positions=positions,input_bands=[(0.5,1.0,1.0)],C=1,L=1,tau=50,e_scale=0.05,ext_scale=0.08)
ar = hrfm.integrate(np.arange(ir[0],ir[1],0.1))#[:,nnodes/2]
cPickle.dump(ar,open("/home/thorsten/Dokumente/eegbold-sim/simple-bold.pck","w+"),2)
cPickle.dump(hrfm._signal_rstc,open("/home/thorsten/Dokumente/eegbold-sim/simple-eeg.pck","w+"),2)

    #p.subplot(3,1,1)
    #pd1 = hrfm._signal_rstc[:]+np.arange(0,nnodes*15,15).reshape(1,-1).repeat(hrfm._signal_rstc.shape[0],axis=0)
    #p.plot(np.arange(ir[0],ir[1],0.1),pd1)
    #p.xlim(ir[0],ir[1]-1)
    #p.ylim(pd1[500:].min(),pd1[500:].max())
    #p.subplot(3,1,2)
    #p.plot(np.arange(ir[0],ir[1],0.1),ar[:,::2])#nnodes/2*2])#+np.arange(0,nnodes*2*15,15).reshape(1,-1).repeat(ar.shape[0],axis=0))
    #p.xlim(ir[0],ir[1]-1)
    #p.ylim(ar[500:,::2].min(),ar[500:,::2].max())
    #p.subplot(3,1,3)
    #wt_power = wt_power_baseline(hrfm._signal_rstc[:,nnodes/2],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(1500,1800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    ##print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="bicubic",extent=[ir[0],ir[1],0.05,1.475],vmin=vmin,vmax=vmax)
    ##p.colorbar()
    #p.xlim(ir[0],ir[1])
    
