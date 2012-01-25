#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Classes and funtions for making symbolic 
representations of time-series"""

import itertools

import eegpy
from eegpy.helper import factorial
from eegpy.misc import FATALERROR, EegpyBase, debug
#from eegpy.analysis.phasespace import make_time_delay_embedding
from eegpy.analysis.entropy.symbols import make_permutation_symbols, make_all_permutation_symbols

try:
    import numpy as N
    np = N
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


###########
# Classes #
###########

#############
# Functions #
#############

def entropy(x,symbols=None):        
    """For a given set of symbols, returns the entropy.
    H = - \sum_i=0^N p_i log(p_i)"""
    ps = prob(x)
    tmp = -(ps*N.log2(ps))
    return np.ma.masked_invalid(tmp).sum()
    

def perm_entropy(x,w=4,t=1,normalize=True):
    """Calculate the permutation-entropy of a time-series x."""
    pe = entropy(make_permutation_symbols(x,w,t))
    if normalize:
        pe=pe/N.log(factorial(w))
    return pe

def joint_entropy(x,y,symbols=None):
    jps = j_prob(x,y,symbols=symbols)
    tmp = -(jps*np.log2(jps))
    return np.ma.masked_invalid(tmp).sum()

def perm_M(x,y,w=4,t=1):
    """Permutation mutual information
    Creates symbolic timeseries from x and y and calculates the mutual information
    """
    s1 = make_permutation_symbols(x,w=w,t=t)
    s2 = make_permutation_symbols(y,w=w,t=t)
    return M(s1,s2,symbols="w%i"%w)

def perm_d(x,y,w=4,t=1):
    """Metric based on mutual information. Normalized to [0,1]
    See http://en.wikipedia.org/wiki/Mutual_information#Metric
    """
    s1 = make_permutation_symbols(x,w=w,t=t)
    s2 = make_permutation_symbols(y,w=w,t=t)
    je = joint_entropy(s1,s2,symbols="w%i"%w) 
    mi = M(s1,s2,symbols="w%i"%w)
    return (je-mi)

def perm_D(x,y,w=4,t=1):
    """Metric based on mutual information. Normalized to [0,1]
    See http://en.wikipedia.org/wiki/Mutual_information#Metric
    """
    s1 = make_permutation_symbols(x,w=w,t=t)
    s2 = make_permutation_symbols(y,w=w,t=t)
    je = joint_entropy(s1,s2,symbols="w%i"%w) 
    mi = M(s1,s2,symbols="w%i"%w)
    return (je-mi)/je

def perm_TE_ensemble(y,x,w=4,tau=1,test_surr=False):
    """Permutation transfer entropy of two sets of timeseries x and y
    Measures information transfer from y to x.
    
    Parameters:
     y: 2d-array, first index timepoints, second index ensemble
     x: equal-size array
     w: word-size for creation of permutation-symbols
     tau: time-delay for creation of permutation-symbols
     test_surr: bool. If a surrogate-test should be perfomed or not
    
    Returns: 1d-array of length y.shape[0]-1"""
    s1 = make_permutation_symbols(y,w=w,t=tau)
    s2 = make_permutation_symbols(x,w=w,t=tau)
    return TE_ensemble(s1, s2)

def M(x,y,symbols=None):
    """Mutual information of two symbolic timeseries x and y"""
    jp = j_prob(x,y,symbols=symbols)
    px = prob(x,symbols=symbols)
    py = prob(y,symbols=symbols)
    rv = 0
    px = px.reshape(-1,1).repeat(jp.shape[0],1)
    py = py.reshape(1,-1).repeat(jp.shape[1],0)
    tmp = jp*N.log2(jp/(px*py))
    rv = N.ma.masked_invalid((tmp*(jp>0).astype("i"))).sum()
    return rv

def TE(y,x,symbols=None):
    """Transfer entropy of two symbolic timeseries x and y
    Measures information transfer from x to y"""
    #rv = j_prob(x[1:],x[:-1],y[:-1])\
    #    *N.log2(c_prob(x[1:],x[:-1],y[:-1])/c_prob(x[1:],x[:-1]))
    jp_xxy = j_prob(x[1:],x[:-1],y[:-1],symbols=symbols)
    jp_xx = j_prob(x[1:],x[:-1],symbols=symbols)
    jp_xy = j_prob(x[:-1],y[:-1],symbols=symbols)
    p_x = prob(x,symbols=symbols)
    #make all the size of jp_xxy to do vectorised calculations
    #print jp_xx.shape, jp_xy.shape, jp_xxy.shape
    jp_xx = N.tile(jp_xx,[jp_xxy.shape[2],1,1])
    jp_xx = jp_xx.swapaxes(0,1)
    jp_xx = jp_xx.swapaxes(1,2)
    jp_xy = N.tile(jp_xy,[jp_xxy.shape[0],1,1])
    p_x = N.tile(p_x,[jp_xxy.shape[0],jp_xxy.shape[2],1])
    p_x = p_x.swapaxes(1,2)
    rv = N.ma.masked_invalid(jp_xxy * N.log10(jp_xxy*p_x/(jp_xy*jp_xx))).sum()
#    for ix1 in range(jp_xxy.shape[0]):
#        for iy in range(jp_xxy.shape[1]):
#            for ix in range(jp_xxy.shape[2]):
#                #try:
#                tmp = jp_xxy[ix1,ix,iy] * N.log10(jp_xxy[ix1,ix,iy]*p_x[ix]/(jp_xy[ix,iy]*jp_xx[ix1,ix]))
#                if N.isfinite(tmp):
#                    rv+=tmp
                #print ix, iy, ix1, rv
                #except Exception, e:
                    #print "error"
                    #pass
    return rv#N.ma.masked_invalid(rv).sum()

def perm_TE(y,x,w=4,t=1):
    """Creates symbolic timeseries from real data and calculates Transfer Entropy"""
    s1 = make_permutation_symbols(x,w=w,t=t)
    s2 = make_permutation_symbols(y,w=w,t=t)
    #return TE(s2,s1,symbols="w%i"%w)
    return TE(s2,s1,symbols=make_all_permutation_symbols(w))

def TE_ensemble(y,x,test_surr=False):
    """Transfer entropy of two sets of symbolic timeseries x and y
    Measures information transfer from x to y.
    
    :Parameters:
    y: 2d-array, first index timepoints, second index ensemble
    x: equal-size array
    test_surr: bool. If a surrogate-test should be perfomed or not
    
    :Returns: 1d-array of length y.shape[0]-1
    """
    rv = N.zeros((y.shape[0]-1),"d")
    for idx in range(y.shape[0]-1):
        jp_xxy = j_prob(x[idx+1,:],x[idx,:],y[idx,:])
        jp_xx = j_prob(x[idx+1,:],x[idx,:])
        jp_xy = j_prob(x[idx,:],y[idx,:])
        p_x = prob(x[idx])
        #make all the size of jp_xxy to do vectorised calculations
        jp_xx = N.tile(jp_xx,[jp_xxy.shape[2],1,1])
        jp_xx = jp_xx.swapaxes(0,1)
        jp_xx = jp_xx.swapaxes(1,2)
        jp_xy = N.tile(jp_xy,[jp_xxy.shape[0],1,1])
        p_x = N.tile(p_x,[jp_xxy.shape[0],jp_xxy.shape[2],1])
        p_x = p_x.swapaxes(1,2)
        #print idx, jp_xxy.shape, jp_xx.shape, jp_xy.shape, p_x.shape
        rv[idx] = N.ma.masked_invalid(jp_xxy * N.log10(jp_xxy*p_x/(jp_xy*jp_xx))).sum() 
#        jp_xx = j_prob(x[idx+1,:],x[idx,:])
#        jp_xy = j_prob(x[idx,:],y[idx,:])
#        p_x = prob(x[idx])

#        for ix1 in range(jp_xxy.shape[0]):
#            for iy in range(jp_xxy.shape[2]):
#                for ix in range(jp_xxy.shape[1]):
#                    tmp = jp_xxy[ix1,ix,iy] * N.log10(jp_xxy[ix1,ix,iy]*p_x[ix]/(jp_xy[ix,iy]*jp_xx[ix1,ix]))
#                    if N.isfinite(tmp):
#                        rv[idx]+=tmp
    return rv


#Probability-functions
def prob(x,symbols=None):
    """Calculate the probabilities for each symbol"""
    #print x
    if symbols==None:
        symbols = np.unique(x)
    elif type(symbols) == str:
        symbols = make_all_permutation_symbols(symbols)
    symbdict = dict([(k,0) for k in symbols])
    for xi in x:
        symbdict[xi] += 1
    ks = symbols
    ks.sort()
    ps = np.array([symbdict[k] for k in ks]).astype(np.float)
    ps /= x.shape[0]
    return ps

def j_prob(x,y,z=None,symbols=None):
    """The joint probability of the symbols in x and y"""
    #print "symbls", symbols
    #extract variables from args
    #print args
    if z==None:
        assert len(x) == len(y), "Sequences must have same length"
    else:
        assert len(x) == len(y) and len(x) == len(z), "Sequences must have same length"
    #Which symbols exist
    if symbols==None:
        symbols = np.unique(np.r_[x,y])
    elif type(symbols) == str:
        symbols = make_all_permutation_symbols(symbols)
    if z==None:
        ps = N.zeros((len(symbols),len(symbols)),"d")
        ks = list(symbols)
        ks.sort()
        idx_dict = dict([(k,ks.index(k)) for k in ks])
        for i in range(len(x)):
            ps[idx_dict[x[i]],idx_dict[y[i]]] +=1
        ps /= x.shape[0]
    else:
        ps = N.zeros((len(symbols),len(symbols),len(symbols)),"d")
        ks = list(symbols)
        ks.sort()
        idx_dict = dict([(k,ks.index(k)) for k in ks])
        for i in range(len(x)):
            ps[idx_dict[x[i]],idx_dict[y[i]],idx_dict[z[i]]] +=1
        ps /= x.shape[0]
    return ps
    
def c_prob(x,y,z=None,symbols=None):
    """Conditional probability. Defined as P(X|Y) = P(X,Y)/P(Y)"""
    jp = j_prob(x,y,z,symbols)
    if z==None:
        ps = jp/jp.mean(axis=1)
    else:
        ps = jp/j_prob(y,z)
        #ps /= ps.mean(axis=1)
    return ps
    
#Helper-functions

def join_symb_lists(li1,li2):    
    """Not used..."""
    rv = []
    for li in [li1,li2]:
        for v in li:
            if not v in rv:
                rv.append(v)
    return rv
        
    
if __name__ == "__main__":
    #from eegpy.analysis.entropy.base import j_prob, prob, c_prob
    #from eegpy.analysis.entropy.symbols import make_permutation_symbols
    from eegpy.models.roessler import TwoRoessler
    import pylab as p
    p.figure(1)
    for i_r,rausch in enumerate(np.arange(0,3,0.3)):
        ar1 = np.sin(np.arange(0,100,0.1)) + np.random.random((1000))*rausch
        ar2 = np.sin(np.arange(0,100,0.1)) + np.random.random((1000))*rausch
        p.subplot(10,1,i_r+1)
        p.plot(ar1)
        p.plot(ar2)
        d = perm_d(ar1,ar2,w=4,t=50)
        D = perm_D(ar1,ar2,w=4,t=50)
        p.text(10,-0.5,"d=%.2f, D=%.2f"%(d,D),bbox={"fc":"w","alpha":0.7})
    p.figure(2)
    for i_e,e in enumerate(np.arange(0,0.4,0.04)):
        tr = TwoRoessler(e1=e,e2=e)
        ar = tr.integrate(np.arange(0,1000,0.1))[-4096:]
        ar1 = ar[:,0]
        ar2 = ar[:,3]
        p.subplot(10,1,i_e+1)
        p.plot(ar1)
        p.plot(ar2)
        d = perm_d(ar1,ar2,w=4,t=50)
        D = perm_D(ar1,ar2,w=4,t=50)
        p.text(10,-0.5,"d=%.2f, D=%.2f"%(d,D),bbox={"fc":"w","alpha":0.7})


    #s1 = make_permutation_symbols(ar,4,1)

    #probs = prob(s1,symbols="w4")
    #print probs, probs.sum()

    #ar3 = N.sin(N.arange(0,5,0.005))
#    eeg = eegpy.F32("/home/thorsten/443_as_fsl.f32","r")
#    t = N.zeros((25,25),"d")
#    for i in range(25):#range(eeg.num_channels):
#        for j in range(25):#range(eeg.num_channels):
#            s1 = make_permutation_symbols(eeg[10000:20000,i],t=10)
#            s2 = make_permutation_symbols(eeg[10000:20000,j],t=10)
#            t[i,j] = TE(s1,s2)
#            print "TE(%i->%i)"%(i,j), t[i,j]
    #eeg = eegpy.open_eeg("/media/Extern/public/Experimente/AudioStroop/2008-07-29_Rieke/eeg/rieke_as.f32")
    #evt = eegpy.EventTable("/media/Extern/public/Experimente/AudioStroop/2008-07-29_Rieke/eeg/rieke_as.vmrk")
    #print evt.keys()
    #ts = evt["S65282"]
    #data = eeg.get_data_for_events(ts, (-200,1000))
    #s1 = make_permutation_symbols(data[:,13,:],t=10)
    #print data.shape, data.mean()
    #s1 = make_permutation_symbols(data[:,13,:],w=3,t=10)
    #s2 = make_permutation_symbols(data[:,14,:],w=3,t=10)
    #w=3
    #tau=10
    #te = N.zeros((data.shape[0]-((w-1)*tau+1),data.shape[1],data.shape[1]),"d")
    #import time
    #t_start = time.time()
    #for ch1 in range(data.shape[1]):
    #    for ch2 in range(data.shape[1]):
    #        te[:,ch1,ch2] = perm_TE_ensemble(data[:,ch1,:], data[:,ch2,:],w=w,tau=tau)
    #        print "<PTE(%i,%i)> ="%(ch1,ch2), te[:,ch1,ch2].mean()
    #print "Duration of calculation:", time.time()-t_start
    
    
            
#    ar2 = ar3[::-1]
#    ar1 = -ar3
#    s1 = make_permutation_symbols(ar1)
#    s2 = make_permutation_symbols(ar2)
#    s3 = make_permutation_symbols(ar3)
#    jp2 = j_prob(s1,s2)
#    jp3 = j_prob(s1,s2,s3)
#    cp2 = c_prob(s1,s2)
#    cp3 = c_prob(s1,s2,s3)
#    TExy = TE(s1,s2)
#    TEyx = TE(s2,s1)
#    print "TExy", TExy,
#    print "TEyx", TEyx
        
        
        
