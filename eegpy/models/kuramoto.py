#!/usr/bin/env python
# -*- coding: utf-8 -*-

from eegpy.misc import FATALERROR
from integrators import rk4int, eulerint
import sys
import time
import pylab as p

try:
    import numpy as np
    from numpy.fft import *
    from scipy.integrate import odeint
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pyosd
except ImportError:
    pyosd=None

class Kuramoto(object):
    """A single coupled Kuramoto-oscillator"""
    def __init__(self,y=None,omega=(1.0,.05)):
        if type(omega)==float:
            self.omega = omega
        elif type(omega[0])==float and type(omega[1])==float:
            self.omega = np.random.normal(loc=omega[0],scale=omega[1])
        else:
            raise ValueError("omega must be float or 2-tuple of floats (loc,scale)")
        if y==None:
            self.y = np.random.random((1))*2*np.pi
        else:
            self.y = y
    
    def ode(self,y,t):
        dy = np.zeros((1))
        dy[0] = self.omega
        return dy

    def integrate(self,ts):
        rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv

class NStochasticKuramoto(object):
    """A single coupled Kuramoto-oscillator"""
    def __init__(self,N_=2,y=None,omegas=(1.0,.05),sigmas=None,ks=None):
        #Number of Kuramotos
        self.N = N_
        #frequencies
        if type(omegas)==np.ndarray:
            if not omegas.shape[0]==N_:
                raise ValueError("Number of omegas is not equal to N")
            self.omegas = omegas
        elif type(omegas)==tuple:
            self.omegas = np.random.normal(loc=omegas[0],scale=omegas[1],size=(N_))
        else:
            raise ValueError("omega must be float or 2-tuple of floats (loc,scale)")
        #Initial values
        if y==None:
            self.y = np.random.random((N_))*2*np.pi
        else:
            self.y = y
        #Noise
        if sigmas==None:
            self.sigmas=np.zeros((N_))
        else:
            if not sigmas.shape[0] == N_:
                raise ValueError("sigmas must have shape (N)")
            self.sigmas=sigmas
        #Couplings
        if type(ks)==np.ndarray:
            if not ks.shape == (N_,N_):
                raise ValueError("If ks are given, must be a NxN-array")
            self.ks = ks
        else: # Make bidirectional chain
            ks = np.zeros((N_,N_))
            for i in range(N_):
                try:
                    ks[N_-1,N_] = np.random.random()*0.05+0.05 # 0.05-0.1
                except Exception,e:
                    pass
                try:
                    ks[N_+1,N_] = np.random.random()*0.05+0.05 # 0.05-0.1
                except Exception,e:
                    pass
            self.ks = ks
    
    def ode(self,y,t):
        dy = self.omegas.copy()
        dy += self.sigmas*np.random.normal(0,1,size=(self.N))
        for i in range(self.N):
            dy += 1./self.N * self.ks[i,:]*np.sin(y-y[i])
        return dy

    def integrate(self,ts):
        rv = rk4int(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return np.sin(rv)

class Kuramoto2ndOrder(object):
    """A single coupled 2nd order Kuramoto-oscillator"""
    def __init__(self,y=None,omega=(1.0,.05)):
        if type(omega)==float:
            self.omega = omega
        elif type(omega[0])==float and type(omega[1])==float:
            self.omega = np.random.normal(loc=omega[0],scale=omega[1])
        else:
            raise ValueError("omega must be float or 2-tuple of floats (loc,scale)")
        if y==None:
            self.y = np.random.random((1))*2*np.pi
        else:
            self.y = y
    
    def ode(self,y,t):
        dy = np.zeros((1))
        dy[0] = self.omega
        return dy

    def integrate(self,ts):
        rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv

if __name__ == "__main__":
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    mpl.rcParams['legend.fontsize'] = 10

    ##Test 1: Ein Kuramoto
    #fig = plt.figure(1)
    #ko = Kuramoto()
    #xs = np.arange(0,10,0.01)
    #ys = ko.integrate(xs)
    #ko2 = Kuramoto()
    #ys2 = ko2.integrate(xs)
    #plt.plot(xs,np.sin(ys),label="omega = %.2f"%ko.omega)
    #plt.plot(xs,np.sin(ys2),label="omega = %.2f"%ko2.omega)
    #plt.legend()

    ##Test 2: Zwei symmetrisch gekoppelte Kuramoto
    from eegpy.analysis.phases import phase_coherence
    ks = np.arange(0.1,0.15,0.05)
    sigmas = np.arange(1.,12,5)
    n_real = 1
    pcs = np.zeros((n_real,sigmas.shape[0],ks.shape[0]))
    for i_k,k in enumerate(ks):
        for i_s,sigma in enumerate(sigmas):
            print "k=%.2f, s=%.2f"%(k,sigma)
            for i_r in range(n_real):
                tr = NStochasticKuramoto(omegas=np.array([0.9,1.1]),sigmas=np.ones((2))*sigma,ks=np.array([[0,k],[k,0]]))
                v = tr.integrate(np.arange(0,70,0.01))[-5000:]
                pcs[i_r,i_s,i_k] = phase_coherence(v[:,0],v[:,1])
                print " R=%.3f" % pcs[i_r,i_s,i_k]
                p.subplot(3,1,i_s+1)
                p.plot(v[:,0])
                p.plot(v[:,1])
                
    #p.imshow(pcs.mean(axis=0),vmin=0,vmax=1)
    #p.bar(np.arange(0.1,ks.shape[0],1),pcs.mean(axis=0),width=0.8)
    #p.errorbar(np.arange(0.5,ks.shape[0],1),pcs.mean(axis=0),pcs.std(axis=0),fmt="none",ecolor="r",elinewidth=2.,capsize=6.)
    #p.xticks(np.arange(0.5,ks.shape[0],1),ks)


