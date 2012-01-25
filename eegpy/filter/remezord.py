#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Supplies remezord method according to Scipy Ticket #475
http://projects.scipy.org/scipy/ticket/475"""
import numpy
from numpy import atleast_1d, poly, polyval, roots, real, asarray, allclose, \
    resize, pi, absolute, logspace, r_, sqrt, tan, log10, arctan, arcsinh, \
    cos, exp, cosh, arccosh, ceil, conjugate, zeros, sinh, hstack, mod
from numpy import mintypecode
from scipy import special, optimize
from scipy.misc import comb
import string, types
 
abs = absolute
 
def findfreqs(num, den, N):
    m = numpy.arange(0,N)
    h = win*special.sinc(cutoff*(m-alpha))
    return h / numpy.sum(h,axis=0)

def oddround(x):
    """Return the nearest odd integer from x."""

    return x-mod(x,2)+1

def oddceil(x):
    """Return the smallest odd integer not less than x."""

    return oddround(x+1)
    
def remlplen_herrmann(fp,fs,dp,ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.

    Uses approximation algorithm described by Herrmann et al.:

    O. Herrmann, L.R. Raviner, and D.S.K. Chan, Practical Design Rules for
    Optimum Finite Impulse Response Low-Pass Digital Filters, Bell Syst. Tech.
    Jour., 52(6):769-799, Jul./Aug. 1973.
    """

    dF = fs-fp
    a = [5.309e-3,7.114e-2,-4.761e-1,-2.66e-3,-5.941e-1,-4.278e-1]
    b = [11.01217, 0.51244]
    Dinf = log10(ds)*(a[0]*log10(dp)**2+a[1]*log10(dp)+a[2])+ \
           a[3]*log10(dp)**2+a[4]*log10(dp)+a[5]
    f = b[0]+b[1]*(log10(dp)-log10(ds))
    N1 = Dinf/dF-f*dF+1

    return int(oddround(N1))

def remlplen_kaiser(fp,fs,dp,ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.

    Uses approximation algorithm described by Kaiser:

    J.F. Kaiser, Nonrecursive Digital Filter Design Using I_0-sinh Window
    function, Proc. IEEE Int. Symp. Circuits and Systems, 20-23, April 1974.
    """

    dF = fs-fp
    N2 = (-20*log10(sqrt(dp*ds))-13.0)/(14.6*dF)+1.0

    return int(oddceil(N2))

def remlplen_ichige(fp,fs,dp,ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.
    
    Uses approximation algorithm described by Ichige et al.:
    
    K. Ichige, M. Iwaki, and R. Ishii, Accurate Estimation of Minimum
    Filter Length for Optimum FIR Digital Filters, IEEE Transactions on
    Circuits and Systems, 47(10):1008-1017, October 2000.
    """
    
    dF = fs-fp
    v = lambda dF,dp:2.325*((-log10(dp))**-0.445)*dF**(-1.39)
    g = lambda fp,dF,d:(2.0/pi)*arctan(v(dF,dp)*(1.0/fp-1.0/(0.5-dF)))
    h = lambda fp,dF,c:(2.0/pi)*arctan((c/dF)*(1.0/fp-1.0/(0.5-dF)))
    Nc = ceil(1.0+(1.101/dF)*(-log10(2.0*dp))**1.1)
    Nm = (0.52/dF)*log10(dp/ds)*(-log10(dp))**0.17 
    N3 = ceil(Nc*(g(fp,dF,dp)+g(0.5-dF-fp,dF,dp)+1.0)/3.0) 
    DN = ceil(Nm*(h(fp,dF,1.1)-(h(0.5-dF-fp,dF,0.29)-1.0)/2.0))
    N4 = N3+DN
    
    return int(N4)

def remezord(freqs,amps,rips,Hz=1,alg='ichige'):
    """Filter parameter selection for the Remez exchange algorithm.

    Description:

      Calculate the parameters required by the Remez exchange algorithm to
      construct a finite impulse response (FIR) filter that approximately
      meets the specified design. 
      
    Inputs:

      freqs --- A monotonic sequence of band edges specified in Hertz. All
                elements must be non-negative and less than 1/2 the
                sampling frequency as given by the Hz parameter.
      amps  --- A sequence containing the amplitudes of the signal to be
                filtered over the various bands.
      rips  --- A sequence specifying the maximum ripples of each band.
      alg   --- Filter length approximation algorithm. May be
                'herrmann', 'kaiser', or 'ichige'.

    Outputs:

      numtaps,bands,desired,weight -- See help for the remez function.   
    """

    # Make sure the parameters are floating point numpy arrays:
    freqs = asarray(freqs,'d')
    amps = asarray(amps,'d')
    rips = asarray(rips,'d')

    # Scale ripples with respect to band amplitudes:
    rips /= (amps+(amps==0.0))

    # Normalize input frequencies with respect to sampling frequency:
    freqs /= Hz

    # Select filter length approximation algorithm:
    if alg == 'herrmann':
        remlplen = remlplen_herrmann
    elif alg == 'kaiser':
        remlplen = remlplen_kaiser
    elif alg == 'ichige':
        remlplen = remlplen_ichige
    else:
        raise ValueError('Unknown filter length approximation algorithm.')
    
    # Validate inputs:
    if any(freqs > 0.5):
        raise ValueError('Frequency band edges must not exceed the Nyquist frequency.')
    if any(freqs < 0.0):
        raise ValueError('Frequency band edges must be nonnegative.')
    if any(rips < 0.0):
        raise ValueError('Ripples must be nonnegative.')
    if len(amps) != len(rips):
        raise ValueError('Number of amplitudes must equal number of ripples.')
    if len(freqs) != 2*(len(amps)-1):
        raise ValueError('Number of band edges must equal 2*((number of amplitudes)-1)')

    # Find the longest filter length needed to implement any of the
    # low-pass or high-pass filters with the specified edges:
    f1 = freqs[0:-1:2]
    f2 = freqs[1::2]
    L = 0
    for i in range(len(amps)-1):
        L = max((L,
                 remlplen(f1[i],f2[i],rips[i],rips[i+1]),
                 remlplen(0.5-f2[i],0.5-f1[i],rips[i+1],rips[i])))

    # Cap the sequence of band edges with the limits of the digital frequency
    # range:
    bands = hstack((0.0,freqs,0.5))

    # The filter design weights correspond to the ratios between the maximum
    # ripple and all of the other ripples:
    weight = max(rips)/rips
    
    return [L,bands,amps,weight]

#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':
    pass