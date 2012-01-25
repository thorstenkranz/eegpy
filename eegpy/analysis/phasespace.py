#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module containing routines for working in phase space. 
Classes:
KdNode: for kd-Trees

methods:

"""
#################
# Module-Import #
#################

#eegpy-modules
try:
    import eegpy
    from eegpy.misc import debug
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')

#Third-party
try:
    import numpy
    #from scipy.signal import lfilter, butter
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


def make_time_delay_embedding(x,m=5,tau=1):
    """Performs time-delay embedding.
    Returns a new array with embedding-vectors. First index is tp, second index is component.
    Up to now only 1d- and 2d-arrays supported."""
    assert len(x.shape) in [1,2], "Only 1d- and 2d-arrays are supported at the moment."
    m = int(round(m))
    tau = int(round(tau))
    assert m>1, "m must be larger than 1"
    assert tau>0, "tau must be larger than 0"
    
    if len(x.shape)==1:
        rv_ar = numpy.zeros((x.shape[0]-(m-1)*tau,m)) #Make the new array
        #print rv_ar.shape
        for mi in range(m): #Clever embedding via slices
            tmp = x[mi*tau:(rv_ar.shape[0])+mi*tau]
            #print rv_ar[:,mi].shape, tmp.shape
            rv_ar[:,mi] = tmp
        return rv_ar
    elif len(x.shape)==2:
        rv_ar = numpy.zeros((x.shape[0]-(m-1)*tau,x.shape[1],m)) #Make the new array
        #print rv_ar.shape
        for i1 in range(x.shape[1]):
            for mi in range(m): #Clever embedding via slices
                tmp = x[mi*tau:(rv_ar.shape[0])+mi*tau,i1]
                #print rv_ar[:,mi].shape, tmp.shape
                rv_ar[:,i1,mi] = tmp
        return rv_ar

def unmake_time_delay_embedding(x,tau=1,mode="mean"):
    """Performs time-delay embedding.
    Returns a new array with embedding-vectors. First index is tp, second index is component.
    Up to now only 1d-arrays supported.
    mode - one of "mean","median","flat" """
    assert len(x.shape) in [2,3], "Only TDE of 1d- and 2d-arrays are supported at the moment."
    if len(x.shape)==2:
        m = x.shape[1]
    else:
        m = x.shape[2]
    tau = int(round(tau))
    assert m>1, "m must be larger than 1"
    assert tau>0, "tau must be larger than 0"
    if len(x.shape)==2:
        if mode in ["mean","median"]:
            tmp = numpy.ones((x.shape[0]+(m-1)*tau,m))*numpy.NaN
            #print rv_ar.shape
            for mi in range(m): #Clever embedding via slices
                tmp[mi*tau:mi*tau+x.shape[0],mi] = x[:,mi]
            #rv_ar = numpy.zeros((tmp.shape[0]),"d")
            #for i in range(rv_ar.shape[0]):
            #    aslice = tmp[i,:]
            if mode=="mean":
                rv_ar = numpy.ma.masked_invalid(tmp).mean(axis=1)
            else:
                rv_ar = np.median(numpy.ma.masked_invalid(tmp),axis=1)
                #print rv_ar[:,mi].shape, tmp.shape
        else: #Old method - 
            rv_ar = numpy.r_[x[:,0],x[-1,1:]]
        return rv_ar
    elif len(x.shape)==3:
        if mode in ["mean","median"]:
            tmp = numpy.ones((x.shape[0]+(m-1)*tau,x.shape[1],m))*numpy.NaN
            #print rv_ar.shape
            rv_ar = numpy.zeros((tmp.shape[0],tmp.shape[1]),"d")
            for i1 in range(x.shape[1]):
                for mi in range(m): #Clever embedding via slices
                    tmp[mi*tau:mi*tau+x.shape[0],i1,mi] = x[:,i1,mi]
                if mode=="mean":
                    rv_ar[:,i1] = numpy.ma.masked_invalid(tmp[:,i1,:]).mean(axis=1)
                else:
                    rv_ar[:,i1] = np.median(numpy.ma.masked_invalid(tmp[:,i1,:]),axis=1)
                    #rv_ar[i,i1] = aslice[True-numpy.isnan(aslice)].mean()
                    #print rv_ar[:,mi].shape, tmp.shape
        else: #Old method - 
            rv_ar = numpy.r_[x[:,:,0],x[-1,:,1:]]
        return rv_ar

def find_n_nearest_neighbors(embed,n=10,fast=True):
    """Searches the n nearest neighbors in state-space.
    Takes an array of TD-Vectors as returned by make_time_delay_embedding
    where first index is tp, second index is component.
    Returns an array of shape (embed.shape[0],n) containing the indices 
    of the nearest neighbors for each tp."""
    if fast:
        try:
            return find_n_nearest_neighbors_fast(embed,n)
        except Exception,e:#ImportError,e:
            print "find_n_nearest_neighbors_fast didn't work. Probably the required library (scikits.ann) is not available"
    #Fallback if scikits.ann is not installed: brute force method
    assert len(embed.shape) == 2, "Only 2d-arrays (one array of TD-Vectors) are supported at the moment."
    n = int(round(n))
    #First, calculate all pairwise distances
    distances = numpy.zeros((embed.shape[0],embed.shape[0]),"d")
    for i in range(embed.shape[0]):
        if debug:
            print i
        for j in range(i):
            #print i, j
            #dist = 0
            #for k in range(embed.shape[1]):
            #    dist += (embed[i,k]-embed[j,k])**2
            #dist = numpy.sqrt(dist)
            #dist = numpy.sqrt(((embed[i,:]-embed[j,:])**2).sum())
            d = embed[i,:]-embed[j,:]
            dist = numpy.dot(d,d)
            distances[i,j] = dist
            distances[j,i] = dist
    #print distances
    #return distances
    #Now, walk through every line and find the n smallest values
    rv_ar = numpy.zeros(((embed.shape[0],n)),numpy.int32)
    BIGNUM = 9e99
    for i in range(distances.shape[0]):
        distances[i,i] = BIGNUM
        for ni in range(n):
            if debug:
                print i,ni
            idx = distances[i].argmin()
            rv_ar[i,ni] = idx
            distances[i,idx] = BIGNUM
    return rv_ar

def find_n_nearest_neighbors_fast(embed,n=10):
    """Searches the n nearest neighbors in state-space.
    Takes an array of TD-Vectors as returned by make_time_delay_embedding
    where first index is tp, second index is component.
    Returns an array of shape (embed.shape[0],n) containing the indices 
    of the nearest neighbors for each tp."""
    assert len(embed.shape) == 2, "Only 2d-arrays (one array of TD-Vectors) are supported at the moment."
    n = int(round(n))
    try:
        import scikits.ann as ann
        k = ann.kdtree(embed)
        rv_ar = k.knn(embed,n+1)[0][:,1:]
    except Exception, e:
        from scipy.spatial import KDTree
        kdt = KDTree(embed)
        rv_ar = kdt.query(kdt.data,k=n)[1][:,1:]
    return rv_ar         
    
    
if __name__ == "__main__":
    pass
