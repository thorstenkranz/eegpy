# -*- coding: utf-8 -*-
"""Try out a new measure for synchronization / directionality.
Calculate the stability of the time delay in the cross correlation function"""

import eegpy
from eegpy.misc import FATALERROR, debug
from eegpy.helper import demean

try:
    import numpy as np
    from scipy.stats import ttest_1samp
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


def xcorr_stability(x1,x2,window_len=200,overlap=190):
    """Go through the arrays x1,x2 in a moving-window-fashion.
    For each window, calculate the crosscorrelation function,
    determine tau (time delay).
    Then, for the tau of all windows, calculate t-stats to test
    H0 that the taus are zero. If not, we have sync (and 
    directionality?)
    """
    if not len(x1.shape)==len(x2.shape)==1:
        raise ValueError("Input arrays must be 1d")
    if not x1.shape[0]==x2.shape[0]:
        raise ValueError("Input arrays must have same length.")
    #TODO: check sanity of window_len and overlap
    ndp = x1.shape[0] #number of datapoints
    taus = np.zeros((ndp-window_len)/(window_len-overlap))
    start = 0
    for i in range(len(taus)): #(0,ndp-window_len,window_len-overlap)):
        ar1 = x1[start:start+window_len]
        ar2 = x2[start:start+window_len]
        xcorr = np.correlate(ar1,ar2,mode="same")
        tau = xcorr.argmax()-xcorr.shape[0]/2
        taus[i]=tau
        start += window_len-overlap
    #print taus, taus.mean(),taus.std()
    t,p = ttest_1samp(taus,0)
    return (t,p), taus, np.tanh(t)


if __name__ == "__main__":
    import pylab as p
    from eegpy.models.roessler import TwoRoessler

    for e in np.arange(-0.2,0.2,0.01):
        if e<0:
            tworoe = TwoRoessler(omega1=1.0,omega2=1.0,e1=-e,e2=0.0)
        else:
            tworoe = TwoRoessler(omega1=1.0,omega2=1.0,e1=0.0,e2=e)
        signal6d = tworoe.integrate(np.arange(100,509.6,0.1))
        signal1 = signal6d[:,0]
        signal2 = signal6d[:,3]
        (t,p),taus,mymeasure = xcorr_stability(signal1,signal2,500,490)
        print e, ">>>", t,p,mymeasure


