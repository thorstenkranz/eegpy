#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR
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
    
########################
# Function definitions #
########################

class ICAnalyzer():
    """Class as a wrapper around mdp.nodes.ICANode for easy use."""
    
    node=None
    data=None
    outdata=None
    
    def __init__(self, x=None):
        if x!=None:
            self.set_data(x)
                
    def set_data(self, x):
        #self.node = mdp.nodes.FastICANode()
        self.node = mdp.nodes.CuBICANode()
        assert x.ndim == 2, "x must be array with ndim=2"
        self.data = x
        self.node.train(x)
        
    def unmix(self,x=None):
        """Do the PCA, either with the data set or custom data."""
        if x==None and self.data != None:
            self.outdata = self.node.execute(self.data)
        elif type(x) == type(n.zeros((1,1))):
            self.outdata = self.node.execute(x)
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
# Conenience-functions #
########################        
_icanalyzer = None #Global Instance of PCAnalyzer for use with methods below

def unmix(x):
    """Do the PCA, with array x."""
    global _icanalyzer
    _icanalyzer = ICAnalyzer(x)
    return _icanalyzer.unmix()
    
        
def mix(y,mask=None):
    """Remix the components"""
    if type(_icanalyzer)==type(ICAnalyzer()):
        return _icanalyzer.mix(y,mask)

        
#######################################
# If called directly, do some example #
#######################################    
if __name__=='__main__':
    pass
