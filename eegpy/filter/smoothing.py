# -*- coding: utf-8 -*-
#from numpy import *

import eegpy
from eegpy.misc import debug
from eegpy.misc import FATALERROR

try:
    import numpy as n
    convolve = n.convolve
except ImportError:
    raise FATALERROR('NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    from scipy.signal import fftconvolve
    convolve = fftconvolve
except ImportError:
    if debug:
        print "Scipy not found. Using numpy.convolve."



_smooth_docstring = """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

def smooth(x,window_len=10,window='hanning'):
    """Double-wrapper around function _smooth for smoothing data. 
    Smoothing is a kind of low-pass filtering.
    
    This function makes a copy of array x and then calls _smooth_recursively, which
    is executed recursively on all dimensions of this array until it is a 1d-Array.
    Then function _smooth is called which does the actual filtering.
    
    %s""" % _smooth_docstring
    y = x.copy()
    _smooth_recursively(y,window_len,window)
    return y

def smooth2d(A, sigma=3):
    #window_len = max(int(sigma), 3)*2+1
    window_len = sigma*2+1
    A1 = n.array([smooth(x, window_len) for x in n.asarray(A)])
    A2 = n.transpose(A1)
    A3 = n.array([smooth(x, window_len) for x in A2])
    A4 = n.transpose(A3)
    return A4

def smooth_windowed_eeg(eeg,channels=[],window_len=10,window='hanning',moving_window_length=10000):
    """For smoothing whole eeg-parts. Uses windowing.
    Will only work for eeg-files, output is array."""
    mwl=moving_window_length
    if channels==None:
        channels = range(eeg.num_channels)
    y = n.zeros((eeg.num_datapoints,len(channels)),"d")
    for win,i_dp in eeg.moving_windows(channels,mwl,window_len*2):
        if i_dp==0:
            y[:mwl-window_len,:] = smooth(win,window_len,window)[:-window_len]
        else:
            y[i_dp+window_len:i_dp+mwl-window_len] = smooth(win,window_len,window)[window_len:-window_len]
    #Now, take care of the last piece of data that doesn't fit into the  moving window
    start = i_dp+mwl-window_len*2
    last_piece = eeg[start:,channels]
    y[start+window_len:] = smooth(last_piece,window_len,window)[window_len:]
    return y

def smooth_windowed_eeg_power(eeg,channels=[],window_len=10,window='hanning',moving_window_length=10000):
    """For smoothing the power of whole eeg-parts. Uses windowing.
    Will only work for eeg-files, output is array."""
    mwl=moving_window_length
    if channels==None:
        channels = range(eeg.num_channels)
    y = n.zeros((eeg.num_datapoints,len(channels)),"d")
    for win,i_dp in eeg.moving_windows(channels,mwl,window_len*2):
        if i_dp==0:
            y[:mwl-window_len,:] = smooth(win**2,window_len,window)[:-window_len]
        else:
            y[i_dp+window_len:i_dp+mwl-window_len] = smooth(win**2,window_len,window)[window_len:-window_len]
    #Now, take care of the last piece of data that doesn't fit into the  moving window
    start = i_dp+mwl-window_len*2
    last_piece = eeg[start:,channels]
    y[start+window_len:] = smooth(last_piece,window_len,window)[window_len:]
    return y
        

def _smooth_recursively(x,window_len,window):
    if len(x.shape)==1:
        x[:]=_smooth(x,window_len,window)
    else:
        for i in range(x.shape[-1]):
            _smooth_recursively(x[...,i],window_len,window)
            
def _smooth(x,window_len=10,window='hanning'):
    """%s""" % _smooth_docstring

    if x.ndim != 1:
        raise ValueError, "_smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is none of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=n.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=n.ones(window_len,'d')
    else:
        w=eval('n.'+window+'(window_len)')

    y=convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


