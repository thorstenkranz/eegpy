#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR

#Third-party
try:
    import numpy
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

try:
    from matplotlib.axes import Subplot
    # uncomment to select /GTK/GTKAgg/GTKCairo
    from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
    from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
    import matplotlib
    #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, NavigationToolbar
    from matplotlib.figure import Figure, SubplotParams
    from matplotlib.axis import Axis
    import matplotlib.cm
except ImportError:
    raise FATALERROR('Error while importing matplotib. Please visit http://matplotlib.sf.net for more information.')

import sys
import time

class MultichannelSubplot(Subplot):
    _chNames = []
    data=None
    def __init__(self, fig, *args, **kwargs):
        Subplot.__init__(self,fig,*args,**kwargs)
        print "Hello, MultichannelSubplot here."
        
              
    def plotChannels(self, x, fmt="k", chNames = None):
        """Takes a 2d-array and plots it."""
        assert len(x.shape)==2, "plotChannels: x is not a 2d-Array"
        if not chNames == None:
            if len(chNames)==x.shape[1]:
                self._chNames = chNames
        if not len(self._chNames) == x.shape[1]:
            self._chNames = range(x.shape[1])
        
        #print x.mean(), x.std(), max(abs(x.max()),abs(x.min()))
        offset = 1#abs(x.mean())+2*x.std()#abs(x).max()
        #print offset
        if offset > 0:
            ytics = [-(y*offset) for y in range(len(self._chNames))]
            self.set_yticks(ytics)
            self.set_yticklabels(self._chNames)
        else:
            self.set_yticks([])
            
        for i in range(x.shape[1]):
            self.data = (x[:,i]-x[:,i].mean())
            self.data /= self.data.std()*2
            self.plot(self.data-i*offset,fmt)
        
        if offset > 0: 
            self.set_ylim(0-(len(self._chNames)-1)*offset-offset,offset)
        else:
            self.set_ylim(-offset,offset)
            
        self.grid(color='#AAAAAA', linestyle='--', linewidth=0.5)
            
        
        
        
              