#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Frequency-filtering for timeseries-data.
The Butterworth-Algorithm from scipy.signal is used to determine filter-coefficients. 
This module supplies functions for filtering forwars and backwards in time in order
to maintain phase-relations.
"""
__docformat__ = "restructuredtext en"

import eegpy
from eegpy.misc import FATALERROR
#from eegpy.filter.filt_misc import filterRecursively
from eegpy.filter.remezord import remezord

#################
# Module-Import #
#################

try:
    import numpy as np
    from numpy import vstack, hstack, eye, ones, zeros, linalg, \
    newaxis, r_, flipud, convolve, matrix, array
    from scipy.signal import lfilter, butter, remez
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    
########################
# Function definitions #
########################

def filtfilt(b,a,x):
    """Double-wrapper around function _filfilt for filtering forwards and backwards in time.
    
This function makes a copy of array x and then calls _filfiltRecursively, which
is executed recursively on all dimensions of this array until it is a 1d-Array.
Then function _filtfilt is called which does the actual filtering.

Based on http://www.scipy.org/Cookbook/FiltFilt, extended to multiple dimensions.

:Parameters:
    - `b`, (array) Filter coefficients to use for filtering. Can be determined using `butter`.
    - `a`: (array) Filter coefficients to use for filtering. Can be determined using `butter`.
    - `x`: (array) data to filter. Can be array of arbitrary dimension. The first index must \
        be the time index.

:Returns:
    - y: (array) the filtered data.
"""

    #filterRecursively(x,_filtfilt,b,a)
    #print len(x.shape)                
    y = np.array(x[:]).copy()
    _filtfilt_recursively(b,a,y)
    return y

def _filtfilt_recursively(b,a,x):
    if len(x.shape)==1:
        x[:]=_filtfilt(b,a,x[:])
    else:
        for i in range(x.shape[-1]):
            _filtfilt_recursively(b,a,x[...,i])


def _filtfilt(b,a,x):
    #For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3

    if x.ndim != 1:
        raise ValueError, "_filtfilt is only accepting 1 dimension arrays."

    #x must be bigger than edge
    if x.size < edge:
        raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."

    if len(a) < ntaps:
        a=r_[a,zeros(len(b)-len(a))]

    if len(b) < ntaps:
        b=r_[b,zeros(len(a)-len(b))]

    zi=__lfilter_zi(b,a)

    #Grow the signal to have edges for stabilizing 
    #the filter with inverted replicas of the signal
    s=r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems 
    # both are needed for filtfilt

    (y,zf)=lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=lfilter(b,a,flipud(y),-1,zi*y[-1])

    return flipud(y[edge-1:-edge+1])

def filtfilt_band(fl,fh,x,Fs=1000.0,border=3):
    """Simple interface for bandpass-filter.

:Parameters:  
    fl: 
        float, lower edge-frequency for filter
    fh: 
        float, upper edge-frequency for filter
    x: 
        array, data to filter. Can be array of arbitrary dimension. The first index must 
        be the time index.
    Fs: 
        float,Sampling-frequency of the data
    border: 
        int, Order of the Butterworth-filter that is constructed. 
        Due to filtering forwards and backwards in time the order of the 
        resulting filter is 2*border
:Returns:
    y: array
        The filtered data
"""
    [b,a]=butter(border,[fl/(Fs/2),fh/(Fs/2)], btype="band")
    #z=lfilter(b,a,xn)
    y=filtfilt(b,a,x)
    return y

def filtfilt_bandstop(fl,fh,x,Fs=1000.0,border=3):
    """Simple interface for bandstop-filter.

:Parameters:  
    fl: 
        float, lower edge-frequency for filter
    fh: 
        float, upper edge-frequency for filter
    x: 
        array, data to filter. Can be array of arbitrary dimension. The first index must 
        be the time index.
    Fs: 
        float,Sampling-frequency of the data
    border: 
        int, Order of the Butterworth-filter that is constructed. 
        Due to filtering forwards and backwards in time the order of the 
        resulting filter is 2*border
:Returns:
    y: array
        The filtered data
"""
    [b,a]=butter(border,[fl/(Fs/2),fh/(Fs/2)], btype="bs")
    #z=lfilter(b,a,xn)
    y=filtfilt(b,a,x)
    return y

def filtfilt_low(fl,x,Fs=1000.0,border=3):
    """Simple interface for lowpass-filter.

:Parameters:  
    fl: 
        float, lower edge-frequency for filter
    x: 
        array, data to filter. Can be array of arbitrary dimension. The first index must 
        be the time index.
    Fs: 
        float, Sampling-frequency of the data
    border: 
        int, Order of the Butterworth-filter that is constructed. 
        Due to filtering forwards and backwards in time the order of the 
        resulting filter is 2*border
:Returns:
    y: array, The filtered data
"""
    [b,a]=butter(border,[fl/(Fs/2)], btype="low")
    #z=lfilter(b,a,xn)
    y=filtfilt(b,a,x)
    return y

def filtfilt_high(fh,x,Fs=1000.0,border=3):
    """Simple interface for highpass-filter.

:Parameters:  
    fh: 
        float, upper edge-frequency for filter
    x: 
        array, data to filter. Can be array of arbitrary dimension. The first index must 
        be the time index.
    Fs: 
        float, Sampling-frequency of the data
    border: 
        int, Order of the Butterworth-filter that is constructed. 
        Due to filtering forwards and backwards in time the order of the 
        resulting filter is 2*border
:Returns:
    y: array, The filtered data
"""
    [b,a]=butter(border,[fh/(Fs/2)], btype="high")
    #z=lfilter(b,a,xn)
    y=filtfilt(b,a,x)
    return y

#Some filters for special situations

def hp_filt_fmri(x,TR,c_t=100,order=3):
    """Used as highpass-filter for fMRI-data. Can be used to the whole 4d-Volume, as the 
    implementation of filtfilt is recursive. The Butterworth-algorithm is used.
    Paramters:
    x: Data to filter, may be 4d. Time must be first index.
    TR: time of repetition of experiment, in seconds
    c_t: time-constant in seconds. The cutoff-frequency is the inverse
    order: order of the filter"""
    assert TR>0.0
    assert c_t>0
    assert type(order) == type(3)
    assert order>0
    #Finished prerequisities, actual filtering
    b,a = butter(order,[1.0/c_t*TR*2],btype="high")
    return filtfilt(b,a,x)

def lowpass_remez(fl1,fl2,x,map=2,mas=40,Fs=1000.0):
    """Simple interface for lowpass-filter using Parks-McClellan algorithm.

:Parameters:  
    fl1: 
        float, lower edge-frequency for transition band
    fl2: 
        float, upper edge-frequency for transition band
    x: 
        array, data to filter. Can be array of arbitrary dimension. The first index must 
        be the time index.
    map: 
        float, maximum attenuation in pass band
    mas: 
        float, minimum attenuation in stop band
    Fs: 
        float, Sampling-frequency of the data
:Returns:
    y: array, The filtered data
"""
    map = float(-abs(map))
    mas = float(-abs(mas))
    n,f,a,w = remezord([fl1,fl2],[1,0],[1-(10**(map/10)),10**(mas/10)],Hz=Fs)
    b=remez(n,f,a,w)
    #z=lfilter(b,a,xn)
    y=filtfilt(b,[1],x)
    return y


#####################
# private functions #
#####################

def __lfilter_zi(b,a):
    #compute the zi state from the filter parameters. see [Gust96].

    #Based on:
    # [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
    # filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
    # Volume 44, Issue 4

    n=max(len(a),len(b))

    zin = (  eye(n-1) - hstack( (-a[1:n,newaxis],
                                 vstack((eye(n-2), zeros(n-2))))))

    zid=  b[1:n] - a[1:n]*b[0]

    zi_matrix=linalg.inv(zin)*(matrix(zid).transpose())
    zi_return=[]

    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
      zi_return.append(float(zi_matrix[i][0]))

    return array(zi_return)

#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':

    from scipy.signal import butter
    from scipy import sin, arange, pi, randn

    from pylab import plot, legend, show, hold

    t=arange(-1,1,.01)
    x=sin(2*pi*t*.5+2)
    #xn=x + sin(2*pi*t*10)*.1
    xn=x+randn(len(t))*0.05
    
    for fl in [0.02,0.1]:#,0.4,0.8]:
        [b,a]=butter(3,[fl,10*fl], btype="band")
    
        z=lfilter(b,a,xn)
        #y=filtfilt(b,a,xn)
        y=_filtfilt(b,a,xn)
    
    
        plot(x,'c')
        hold(True)
        plot(xn,'k')
        plot(z,'r')
        plot(y,'g')

    legend(('original','noisy signal','lfilter - butter 3 order','filtfilt - butter 3 order'))
    hold(False)
    show()
