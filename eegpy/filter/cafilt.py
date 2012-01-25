#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module should implement the recursive filtering for arbitrary dimensional arrays"""

try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.filter.pcafilt import PCAnalyzer
    from eegpy.filter.icafilt import ICAnalyzer
#from eegpy.filter.filt_misc import filterRecursively
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')


#################
# Module-Import #
#################

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


class RecursiveCAnalyzer():
    """CLASS"""
    
    analyzer_ar = None
    analyzer = None
    data=None
    outdata=None
    
    def __init__(self, type_="PCA", x=None):
        assert type_ in ["PCA","ICA"], "Unknown type of component-analysis chosen. Choose one of PCA or ICA."
        if x!=None:
            try:
                test = len(x.shape)
                self.data = x
            except AttributeError, e:
                raise AttributeError(e)
            
        if type_ == "PCA":
            self.analyzer = PCAnalyzer
        elif type_ == "ICA":
            self.analyzer = ICAnalyzer
        
            
    def unmix(self,x=None):
        if x!=None:
            try:
                test = len(x.shape)
                self.data = x
            except AttributeError, e:
                raise AttributeError(e)
        
        #Todo: Das rekursive Filtern muss irgendwie geregelt werden.