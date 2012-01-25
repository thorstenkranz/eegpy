#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Classes and funtions for making symbolic 
representations of time-series"""

from copy import copy
import itertools

import eegpy
from eegpy.misc import FATALERROR, EegpyBase, debug
from eegpy.analysis.phasespace import make_time_delay_embedding

try:
    import numpy as N
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


###########
# Classes #
###########

#############
# Functions #
#############

def make_permutation_symbols(x, w=4, t=1):        
    y = make_time_delay_embedding(x,w,t)
    pot_of_10 = N.array([10**x10 for x10 in range(w)])
    if len(x.shape)==1:
        for i in range(y.shape[0]):
            idx = y[i,:].argsort()
            y[i,0] = (idx*pot_of_10).sum()
            #y[i,0] = float("".join([str(x) for x in idx]))
            #y[i,i1,0] = hash(str(ixd))
    elif len(x.shape)==2:
        idx = y[:,:,:].argsort(axis=2)
        for i1 in range(y.shape[1]):
            for i in range(y.shape[0]):
                #idx = y[i,i1,:].argsort()
                y[i,i1,0] = (idx[i,i1,:]*pot_of_10).sum() #float("".join([str(x) for x in idx[i,i1,:]]))
                #y[i,i1,0] = hash(str(idx))
    y = y[...,0]
    return y#, symbols
        
def make_all_permutations(li):
    if len(li)>1:
        perms = []
        lower_perms = make_all_permutations(li[1:])
        for lp in lower_perms:
            for j in range(len(lp)+1):
                perm = copy(lp)
                perm.insert(j,li[0])
                perms.append(perm)
    else:
        perms = [li]
    #print "rv:", perms
    return perms

def make_dict_of_permutations(w=4):
    li = range(10,10+w)
    all_p = make_all_permutations(li)
    rv = {}
    for i,p in enumerate(all_p):
        print i, p
        rv[float("".join(map(str,p)))] = i
    return rv

def make_all_permutation_symbols(w):
    if type(w)==str:
        if w.startswith("w"): #make symbols as in make_permutation_symbols
            w=int(w[1:])
    if type(w)==int:
        symbols = []
        pot_of_10 = N.array([10**x10 for x10 in range(w)])
        for ii in itertools.permutations(range(w)):
            symbols.append((N.array(ii)*pot_of_10).sum())
        return symbols

    else:
        raise ValueError("Definition of permutation symbols not understood")
        
        
        
        
        
        
        
        
