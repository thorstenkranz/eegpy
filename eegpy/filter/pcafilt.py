#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR
from eegpy.misc import debug
#from eegpy.filter.filt_misc import filterRecursively

#################
# Module-Import #
#################

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import mdp
except ImportError:
    raise FATALERROR('Modular Toolkit for Dataprocessing in Python not found!\nPlease visit mdp.sf.net for more information.')
    
################
# Gerneral PCA #
################

class PCAnalyzer():
    """Class as a wrapper around mdp.nodes.PCANode for easy use. For general use."""
    
    node=None
    data=None
    outdata=None
    
    def __init__(self, x=None):
        if x!=None:
            assert len(x.shape) == 2, "x is not an array"
            self.set_data(x)
                
    def set_data(self, x):
        self.node = mdp.nodes.PCANode(svd=True)
        assert len(x.shape) == 2, "x is not an array"
        self.data = x
        self.node.train(x)
        
    def unmix(self,x=None):
        """Do the PCA, either with the data set or custom data."""
        if x==None and self.data != None:
            try:
                self.outdata = self.node.execute(self.data)
            except mdp.NodeException, ne:
                if debug:
                    print "Error while unmixing,", ne
        elif type(x) == type(n.zeros((1,1))):
            try:
                self.outdata = self.node.execute(x)
            except NodeException, ne:
                if debug:
                    print "Error while unmixing,", ne
        return self.outdata
        
    def mix(self,y=None,mask=None):
        """Remix the components"""
        if y == None:
            y=self.outdata
        y = y.copy()
        if mask != None:
            assert len(mask) == y.shape[1]
            for i,m in enumerate(mask):
                y[:,i] *= float(m)
        return self.node.inverse(y)
            

########################
# Convenience-functions #
########################        
_pcanalyzer = None #Global Instance of PCAnalyzer for use with methods below

def unmix(x):
    """Do the PCA, with array x."""
    global _pcanalyzer
    _pcanalyzer = PCAnalyzer(x)
    return _pcanalyzer.unmix()
    
        
def mix(y,mask=None):
    """Remix the components"""
    if type(_pcanalyzer)==type(PCAnalyzer()):
        return _pcanalyzer.mix(y,mask)

###############################
# Method of optimal basis set #
###############################


        
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
        y=filtfilt_band(0.01,0.05,xn, border=1)
    
    
        plot(x,'c')
        hold(True)
        plot(xn,'k')
        plot(z,'r')
        plot(y,'g')

    legend(('original','noisy signal','lfilter - butter 3 order','filtfilt - butter 3 order'))
    hold(False)
    show()