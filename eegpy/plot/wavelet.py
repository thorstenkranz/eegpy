#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for plotting wavelet-coefficients for complete 
wavelet-decomposition.
"""

__docformat__ = "restructuredtext en"

#################
# Module-Import #
#################

#eegpy-modules
import eegpy
from eegpy.misc import FATALERROR, debug

#Third-party
try:
    import numpy as N
    from scipy.stats import scoreatpercentile
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pylab as P
    from matplotlib import cm
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

def plot_wavedec_lin(coeff,Fs=1000.0,log_scale=False,**kwargs):
    """
    Plots the coefficients of a 1d wavelet-decomposition, as returned by 
    method eegpy.analysis.wavelet.wavedec_lin using imshow.
    The main work is to rewrite the coefficients as 2d-array. 
    """
    #print "coeff:", type(coeff)
    if not len(coeff.shape)==1:
        raise ValueError("coeff must be an 1d-array")
    length = coeff.shape[0]
    if not N.log2(length)==int(N.log2(length)): #length is not a power of two
        raise ValueError("Length of coeff must be a power of 2")
    if kwargs=={}:
        kwargs= {"origin": "upper",
                 "aspect": "auto",
                 "interpolation": "nearest",
                 }
    #Start building array, handle case of masked arrays
    if isinstance(coeff,N.ma.MaskedArray):
        data = make_plot_array(coeff.data,log_scale)
        mask = make_plot_array(coeff.mask,log_scale)
        plot_ar = N.ma.MaskedArray(data,mask)
    else:
        plot_ar = make_plot_array(coeff,log_scale)
    if not kwargs.has_key("vmin"):
        kwargs["vmin"]=scoreatpercentile(plot_ar.flatten(),5)
    if not kwargs.has_key("vmax"):
        kwargs["vmax"]=scoreatpercentile(plot_ar.flatten(),95)
    if not kwargs.has_key("extent"):
        kwargs["extent"]=[0,length*1000/Fs,0,Fs/2]
    #extent = kwargs["extent"]
    #P.subplot(9,1,9)
    #print "ar:", type(ar)
    P.imshow(plot_ar.T,**kwargs)
    if log_scale:
        #TODO: implement yticks
        print "yticks for log_scale=False not implemented yet"
        P.yticks([])
    #P.imshow(plot_ar.T, aspect="auto", vmin=-2, vmax=2)
    return plot_ar
                
def make_plot_array(coeff,log_scale=False):
    """Construct the array for plotting"""
    length = coeff.shape[0]
    i = 0
    pot = length/2
    while pot>1:
        #print "i,pot", i,pot
        ar = coeff[pot:2*pot].copy()
        #print "ar:", type(ar)
        if i==0:
            plot_ar = ar.reshape((-1,1))
            if not log_scale:
                plot_ar = plot_ar.repeat(pot/2,axis=1)
            #plot_ar.repeat(pot/2,axis=1)
        else:
            ar = ar.repeat(2**(i)).reshape(-1,1)
            if not log_scale:
                ar = ar.repeat(pot/2,axis=1)
            #print plot_ar.shape, ar.shape
            plot_ar = N.concatenate([plot_ar,ar],axis=1)
        #print "Shape:", plot_ar.shape
        #P.subplot(9,1,i+1)
        #P.imshow(plot_ar.T, aspect="auto", vmin=-2, vmax=2)
        i+=1
        pot/=2
    return plot_ar


if __name__=="__main__":
    for i in range(4,8):
        P.subplot(4,1,i-4+1)
        ar = N.arange(2**i)
        plot_ar = plot_wavedec_lin(ar,log_scale=True)
        #P.ylim((450,500))
    P.show()


