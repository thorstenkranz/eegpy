# -*- coding: utf-8 -*-

"""Plotting parts of eegs"""

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
    from scipy.signal import detrend
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pylab as P
except ImportError:
    raise FATALERROR('Pylab not found.\nPlease visit matplotlib.sf.net for further information')

def plot_eeg_segment(data,channel_names,xs=None,offset=100,do_detrend=False,ax=P,*args,**kwargs):
    """Plots single segment of eeg data
    """
    if do_detrend:
        data = detrend(data,axis=0,type="constant")
    plotdata = N.zeros(data.shape,"d")
    for d in range(data.shape[1]):
        plotdata[:,d] = data[:,d]-d*offset 
    if xs==None:
        xs = N.arange(data.shape[0])
    ax.plot(xs,plotdata,*args,**kwargs)
    ax.yticks(-offset*N.arange(len(channel_names)),channel_names)
    ax.ylim((-offset*data.shape[1]+1,offset))
    ax.xlim((xs[0],xs[-1]))

