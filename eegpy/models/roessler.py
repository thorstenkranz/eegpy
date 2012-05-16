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

def roessler_ode(y,t,omega=1,a=0.165,b=0.2,c=10):
    dy = np.zeros((3))
    dy[0] = -1.0*(omega*y[1] + y[2]) #+ e1*(y[3]-y[0])
    dy[1] = omega * y[0] + a * y[1]
    dy[2] =  b + y[2] * (y[0] - c)
    return dy

class Roessler(object):
    """A single coupled Roessler oscillators"""
    def __init__(self, y=None, omega=1.0, a=0.165,b=0.2,c=10):
        self.omega = omega
        self.a = a
        self.b = b
        self.c = c
        if y==None:
            self.y = np.random.random((3))+0.5
        else:
            self.y = y

    def ode(self,y,t):
        dy = roessler_ode(y[:],t,self.omega,self.a,self.b,self.c)
        return dy

    def integrate(self,ts):
        rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv

class TwoRoessler(object):
    """Two coupled Rössler oscillators"""
    def __init__(self, y=None, omega1=1.03, omega2=0.97, e1=0, e2=0, a1=0.165, a2=0.165, b1=0.2, b2=0.2, c1=10, c2=10):
        self.omega1 = omega1
        self.omega2 = omega2
        self.e1 = e1
        self.e2 = e2
        self.a1 = a1
        self.b1 = b1
        self.c1 = c1
        self.a2 = a2
        self.b2 = b2
        self.c2 = c2
        if y==None:
            self.y = np.random.random((6))+0.5
        else:
            self.y = y

    def ode(self,y,t):
        dy = np.zeros((6))
        dy[:3] = roessler_ode(y[:3],t,self.omega1,self.a1,self.b1,self.c1)
        dy[3:] = roessler_ode(y[3:],t,self.omega2,self.a2,self.b2,self.c2)
        dy[0] += self.e1 * (y[3]-y[0])
        dy[3] += self.e2 * (y[0]-y[3])
        return dy

    def integrate(self,ts):
        rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv
        
class NRoessler(object):
    """Implements a network of N coupled Rössler oscillators"""
    def __init__(self, omegasOrN = None, es=None, as_=None, bs=None, cs=None):
        if omegasOrN == None and es == None:
            raise ValueError("I need either omegas, epsilons or number of oscillators to go on.")
        if type(omegasOrN) == type(np.ones((2))):
            self.omegas = omegasOrN
            self.N = self.omegas.shape[0]
            if es != None:
                assert self.omegas.shape[0] == es.shape[0]
                self.es = es
            else:
                self.es = np.zeros((self.N,self.N))
            self.y = np.random.random((3*self.N))+0.5
        elif omegasOrN == None:
            self.es = es
            self.N = self.es.shape[0]
            self.omegas = np.random.normal(1.0,0.03,(self.N))
            self.y = np.random.random((3*self.N))+0.5
        elif type(omegasOrN) == type(1):
            self.N = omegasOrN
            self.omegas = np.random.normal(1.0,0.03,(self.N))
            self.es = np.zeros((self.N,self.N))
            self.y = np.random.random((3*self.N))+0.5
        else:
            raise ValueError("Illegal combination of paramters")
        # Store attractor-parameters
        if as_ == None:
            self._as = np.ones((self.N))*0.165
        else:
            self._as = as_
        if bs == None:
            self._bs = np.ones((self.N))*0.2
        else:
            self._bs = bs
        if cs == None:
            self._cs = np.ones((self.N))*10.
        else:
            self._cs = cs

    def ode(self,y,t):
        dy = np.zeros((self.N*3))
        for i in range(self.N):
            # Dynamics
            dy[(i*3):((i+1)*3)] = roessler_ode(y[(i*3):((i+1)*3)],t,self.omegas[i],self._as[i],self._bs[i],self._cs[i]) 
            # Couplings
            for j in range(self.N):
                if i!=j:
                    dy[i*3] += self.es[i,j] * (y[j*3]-y[i*3])
        return dy

    def integrate(self,ts):
        rv = odeint(self.ode,self.y,ts)#,rtol=1e-4,atol=1e-8) #,hmin=0.0001)
        self.y = rv[-1,:]
        return rv

class NRoesslerTransientCouplings(NRoessler):
    """Network of N coupled Roessler oscillators. The first one is regarded as external oscillator.
        During the transient disturbations, the coupling from this external oscillator is enabled.
        The external oscillator is reset at the beginning of each transient coupling.
    """
    def __init__(self, omegasOrN = None, es=None, as_=None, bs=None, cs=None, transient_couplings = [], reset_offset = 100):
        """Like NRoessler, but with transient forces. 
            trans_couplings: list of tuples (t_start, t_end, factor)
        """
        NRoessler.__init__(self,omegasOrN,es,as_,bs,cs)
        self.original_es = self.es.copy()
        print "original_es", self.original_es
        #for tf in transient_couplings:
        #    assert tf[0]>=0 and tf[0]<self.N
        self.transient_couplings = transient_couplings
        self.__last_t = 0
        self._reset_offset = reset_offset
        self._external_reset = None
        self._external_reset_power = 5
    
    def ode(self,y,t):
        # set reset values?
        if self._external_reset == None:
            if t>self._reset_offset:
                print "Setting reset_offset:", t, y[:3]
                self._external_reset = y[:3].copy()
        #adapt epsilons, with or without external oscillator?
        self.es[:,0] = np.zeros((self.N))
        do_pull = False
        for tf in self.transient_couplings:
            if t>=tf[0] and t<tf[1]:
                self.es[:,0] = tf[2] * self.original_es[:,0]
                if abs(t-tf[0])<1:
                    do_pull = True
                break
        #then normal dynamics
        dy = NRoessler.ode(self,y,t)
        # pulling external oscillator to reset value
        if self._external_reset != None:
            if do_pull:
                #print "pulling ext. osc. to", self._external_reset, "now:", y[:3] 
                dy[:3] = (self._external_reset[:3]-y[:3])*self._external_reset_power
        if time.time()-self.__last_t > 10:
            print t
            self.__last_t = time.time()
        return dy


class RoesslerBrainModel(NRoesslerTransientCouplings):
    """This system is designed to be a model for (visual) stimulus processing in the brain.
    It consists of (i) an external oscillator, during stimulus coupled to (ii) a collection
    of coupled oscillators repr. the visual cortex. From there, connections exist to (iii + iv) two other
    clusters of oscillators, repr. higher order brain structures. Each of these clusters has stronger
    internal couplings and some connections to the other cluster. No couplings go back to VC.
    """
    def __init__(self,n_vc=20,n_c1=20,n_c2=20,e_v1=0.3,e_v2=0.3,e_12=0,e_21=0,transient_couplings=[],*args,**kwargs):
        n_all = 1+n_vc+n_c1+n_c2 # Gesamtzahl von Oszillatoren
        #Set connectivities
        epsilons = np.zeros((n_all,n_all),"d")
        epsilons[1+n_vc:1+n_vc+n_c1,1:1+n_vc] = np.random.random((n_c1,n_vc))*e_v1 # VC -> C1
        epsilons[1+n_vc+n_c1:1+n_vc+n_c1+n_c2,1:1+n_vc] = np.random.random((n_c2,n_vc))*e_v2 # VC -> C2
        epsilons[1+n_vc:1+n_vc+n_c1,1+n_vc+n_c1:1+n_vc+n_c1+n_c2] = np.random.random((n_c1,n_c2))*e_21 # C2 -> C1
        epsilons[1+n_vc+n_c1:1+n_vc+n_c1+n_c2,1+n_vc:1+n_vc+n_c1] = np.random.random((n_c2,n_c1)) # C1 -> C2
        NRoesslerTransientCouplings.__init__(self, n_all, es=epsilons, transient_couplings = transient_couplings, *args, **kwargs) 
        self.original_es[1+n_vc:] = 0 # Keine Kopplung von externem Oszillator auf C1 / C2

        
       
class TwoRoesslerPulseCoupling(TwoRoessler):
    """Same as normal coupled Roessler Oscillators, but other coupling.
    A pulse-coupling, like in neurons, is implemented. The z-component 
    of the driver makes the responder move outwards, away from the 
    center of rotation.
    After each 'flip' in z-direction, each responder has a 'dead time' 
    where it doesn't react to external driving."""
    def __init__(self, y=None, omega1=1.03, omega2=0.97, e1=0, e2=0, deadtime1=5, deadtime2=5):
        TwoRoessler.__init__(self,y,omega1,omega2,e1,e2)
        self.dtime1 = deadtime1
        self.dtime2 = deadtime2
        self.fire_ts = [0,0]
        # Parameters which are changeable but not in contructor
        self.z_thres = 3 # when is a flip in z-direction regarded as firing?

    def ode(self,y,t):
        """Implements the dynamics. First, normal Roessler. 
        Then coupling: if not in deadtime, the additional dy is in 
        direction of (dy[0],dy[1]) and is scaled by dy[5], the 
        z-component of the other oscillator"""
        dy = np.zeros((6))
        dy[:3] = roessler_ode(y[:3],t,self.omega1)
        dy[3:] = roessler_ode(y[3:],t,self.omega2)
        #Coupling; Movement in xy-plane
        dxy = y[0:2]/np.sqrt(np.sum(y[0:2]**2)) * self.e1 * y[5] * int((t-self.fire_ts[0])>self.dtime1)
        #print dxy, y[0:2]
        #print dxy, dy[0:2], np.sqrt(np.sum(dy[0:2]**2)), self.e1, y[5], int((t-self.fire_ts[0])>self.dtime1)
        #print "dxy:", dxy,
        dy[0] += dxy[0]
        dy[1] += dxy[1]
        dxy = y[3:5]/np.sqrt(np.sum(y[3:5]**2)) * self.e2 * y[2] * int((t-self.fire_ts[1])>self.dtime2)
        #print " ", dxy, y[3:5]
        #time.sleep(0.1)
        #print " ", dxy, dy[3:5], np.sqrt(np.sum(dy[3:5]**2)), self.e1, y[2], int((t-self.fire_ts[1])>self.dtime2)
        #print "dxy:", dxy
        dy[3] += dxy[0]
        dy[4] += dxy[1]
        if abs(y[2]-self.z_thres)<0.1 and dy[2]>0: #when z goes up and near thres: start dead-time
            self.fire_ts[0] = t
            #print t
        if abs(y[5]-self.z_thres)<0.1 and dy[5]>0: #when z goes up and near thres: start dead-time
            self.fire_ts[1] = t
            #print " ", t
        return dy
    #TODO: implement!

class StochasticRoessler(Roessler):
    """Roessler-Oszillator with intrinsic noise.
    Uses rk4 instead of odeint.
    """
    def __init__(self, y=None, omega=1.0, a=0.165,b=0.2,c=10,sigma=1.5):
        Roessler.__init__(self,y,omega,a,b,c)
        if not 0<=sigma:
            raise ValueError("StochasticRoessler: sigma must be non-negative float")
        self.sigma = sigma

    def ode(self,y,t):
        dy = roessler_ode(y[:],t,self.omega,self.a,self.b,self.c)
        dy[0] += self.sigma*np.random.normal(0,1) 
        return dy

    def integrate(self,ts):
        #rv = eulerint(self.ode,self.y,ts,h=10**-2)
        rv = rk4int(self.ode,self.y,ts,h=1e-3)
        #rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv

class TwoStochasticRoessler(TwoRoessler):
    """Two coupled stochastic Roessler oscillators"""
    def __init__(self, y=None, omega1=1.03, omega2=0.97, e1=0, e2=0, a1=0.165, a2=0.165, b1=0.2, b2=0.2, c1=10, c2=10, sigma1=1.5, sigma2=1.5):
        TwoRoessler.__init__(self,y,omega1,omega2,e1,e2,a1,a2,b1,b2,c1,c2)
        if not 0<=sigma1:
            raise ValueError("StochasticRoessler: sigma1 must be non-negative float")
        if not 0<=sigma2:
            raise ValueError("StochasticRoessler: sigma2 must be non-negative float")
        self.sigma1 = sigma1
        self.sigma2 = sigma2

    def ode(self,y,t):
        dy = np.zeros((6))
        dy[:3] = roessler_ode(y[:3],t,self.omega1,self.a1,self.b1,self.c1)
        dy[3:] = roessler_ode(y[3:],t,self.omega2,self.a2,self.b2,self.c2)
        dy[0] += self.e1 * (y[3]-y[0])
        dy[3] += self.e2 * (y[0]-y[3])
        dy[0] += self.sigma1*np.random.normal(0,1) 
        dy[3] += self.sigma2*np.random.normal(0,1) 
        return dy

    def integrate(self,ts):
        #rv = eulerint(self.ode,self.y,ts,h=10**-2)
        rv = rk4int(self.ode,self.y,ts,h=1e-3)
        #rv = odeint(self.ode,self.y,ts)
        self.y = rv[-1,:]
        return rv

        

#if __name__ == "__main__":
    #import matplotlib as mpl
    #from mpl_toolkits.mplot3d import Axes3D
    #import matplotlib.pyplot as plt

    #mpl.rcParams['legend.fontsize'] = 10
    ############################


    ##Test 1: Zwei unabhängige Rössler
    #fig = plt.figure(1)
    ##p.subplot(121)
    #ax = fig.gca(projection='3d')
    ##theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    ##z = np.linspace(-2, 2, 100)
    ##r = z**2 + 1
    ##x = r * np.sin(theta)
    ##y = r * np.cos(theta)
    #v = odeint(roessler_ode, np.random.random((3))*2-1, np.arange(1000,1200,0.01),args=(1.03,0.165,0.2,10))
    #x1 = v[:,0]
    #y1 = v[:,1]
    #z1 = v[:,2]
    #ax.plot(x1, y1, z1, label='omega = 1.03')
    #v = odeint(roessler_ode, np.random.random((3))*2-1, np.arange(1000,1200,0.01),args=(0.97,0.165,0.2,10))
    #x2 = v[:,0]
    #y2 = v[:,1]
    #z2 = v[:,2]
    #ax.plot(x2-30, y2, z2, label='omega = 0.97')
    #ax.legend() 
    #plt.figure(2)
    ##p.subplot(122)
    #p.plot(x1[::10])
    #p.plot(x2[::10])
    #plt.show()

    #Test 2: Zwei symmetrisch gekoppelte Rössler
    #fig = plt.figure(1)
    #plt.subplot(121)
    ##ax = fig.gca(projection='3d')
    #tr = TwoRoessler(e1=0.01,e2=0.01)
    #v = tr.integrate(np.arange(0,2200,0.01))[-20000:]
    #x1 = v[:,0]
    #y1 = v[:,1]
    #z1 = v[:,2]
    #x2 = v[:,3]
    #y2 = v[:,4]
    #z2 = v[:,5]
    #ax.plot(x1, y1, z1, label='omega = 1.03')
    #ax.plot(x2-30, y2, z2, label='omega = 0.97')
    #ax.legend() 
    #plt.figure(2)
    #p.subplot(122)
    #p.plot(x1[::10])
    #p.plot(x2[::10])
    #plt.show()
   # 
    #from eegpy.analysis.phases import phase_coherence
    #for e in np.arange(0.0,0.16,0.02):
    #    print "e=%f"%e
    #    tr = TwoRoessler(e1=e,e2=e)
    #    v = tr.integrate(np.arange(0,2200,0.01))[-20000:]
    #    pc = phase_coherence(v[:,0],v[:,3])
    #    print "R=%.3f"%pc

    ##Test 3: Zehn Rössler mit 2 Treibern
    #fig = plt.figure(1)
    #epsilons = np.zeros((10,10))
    #epsilons[:,3] = np.ones((10))*0.1
    #epsilons[:,-3] = np.ones((10))*0.2
    #epsilons[3,-3] = epsilons[-3,3] = 0
    #print "epsilons", epsilons
    #nr = NRoessler(es=epsilons)
    #v = nr.integrate(np.arange(0,2200,0.01))[-20000:]
    #from eegpy.analysis.phases import phase_coherence
    #pcs = np.ones((10,10))
    #for i in range(10):
    #    for j in range(i):
    #        pcs[i,j] = pcs[j,i] = phase_coherence(v[:,i*3],v[:,j*3])
    #p.hot()
    #p.imshow(pcs,interpolation="nearest",vmin=0,vmax=1)

    ##Test 4: Zehn Rössler mit 2 Treibern und transienter Kraft
    #n_roe = 30#30
    #omegas = np.random.normal(0.89,0.05,(n_roe))
    #epsilons = np.zeros((n_roe,n_roe))
    ##epsilons = np.random.random((n_roe,n_roe))*0.003
    #epsilons[2:10,2:10] = np.random.random((8,8))*0.03
    #epsilons[20:25,20:25] = np.random.random((5,5))*0.03
    #epsilons[0,:] = np.zeros((n_roe))
    #epsilons[:,0] = np.random.random((n_roe))*0.03+0.03
    ##epsilons[:,3] = np.ones((10))*0.1
    ##epsilons[:,-3] = np.ones((10))*0.06
    ##epsilons[3,-3] = epsilons[-3,3] = 0
    #print "omegas", omegas
    #print "epsilons", epsilons
    #transient_couplings = [(300,340,5.0)]
    #for i in range(7):
    #    offset = np.random.randint(80,120)
    #    transient_couplings.append((transient_couplings[-1][0]+offset, transient_couplings[-1][1]+offset,5.0))
    #nr = NRoesslerTransientCouplings(omegas,es=epsilons,transient_couplings=transient_couplings, reset_offset=100)
    #v = nr.integrate(np.arange(0,1100,0.1))
    #print v.shape
    #p.figure(1)
    #for tc in transient_couplings:
    #    p.axvspan(tc[0],tc[1],color="b",alpha=0.1)
    #for i in range(n_roe):
    #    p.plot(np.arange(0,1100,0.1),v[:,i*3]-i*30)
    #p.fill_between(np.arange(0,1100,0.1),v[:,::3].mean(axis=1)-(i+1)*30-v[:,::3].std(axis=1),v[:,::3].mean(axis=1)-(i+1)*30+v[:,::3].std(axis=1),alpha=0.2)
    #p.plot(np.arange(0,1100,0.1),v[:,::3].mean(axis=1)-(i+1)*30,"k-",lw=2)
    #p.figure(2)
    #from eegpy.analysis.phases import phase_coherence
    #coherence = np.ones((n_roe,n_roe))
    #for i in range(n_roe):
    #    for j in range(i):
    #        coherence[i,j] = coherence[j,i] = phase_coherence(v[:,i*3],v[:,j*3])
    #p.hot()
    #p.imshow(coherence,vmin=0,vmax=1,interpolation="nearest")
    
    ##Test 5: Stochastic Roessler
    #fig = plt.figure(1)
    #tr = StochasticRoessler(sigma=1.5)
    ##tr = Roessler()
    ##assert 1==0
    #v = tr.integrate(np.arange(0,200,0.01))[-20000:]
    ##assert 1==0
    #x1 = v[:,0]
    #y1 = v[:,1]
    #z1 = v[:,2]
    #plt.plot(x1, label='x')
    #plt.plot(y1, label='y')
    #plt.plot(z1, label='z')
    #plt.legend() 
    ##plt.figure(2)
    ##plt.show()

    ##Test 6: Two Stochastic Roessler
    #fig = plt.figure(1)
    #tr = TwoStochasticRoessler(e1=0,e2=0.05,sigma1=1.5,sigma2=1.5)
    ##tr = Roessler()
    ##assert 1==0
    #v = tr.integrate(np.arange(0,500,0.01))[-20000:]
    ##assert 1==0
    #x1 = v[:,0]
    #y1 = v[:,1]
    #z1 = v[:,2]
    #x2 = v[:,3]
    #plt.plot(x1, label='x1')
    ##plt.plot(y1, label='y')
    ##plt.plot(z1, label='z')
    #plt.plot(x2, label='x2')
    #plt.legend() 
    ##plt.figure(2)
    ##plt.show()

    #Test 7: Two Roessler, one with damping, transient coupling
    #roe = NRoesslerTransientCouplings(omegasOrN=np.array([1.0,3.0]),es=np.array([[0,0.0],[3.5,0]]),as_=[0.165,-0.05],transient_couplings=[(300,350,1.)])
    #v = roe.integrate(np.arange(10,1000,0.5))
    #p.plot(np.arange(10,1000,0.5),v[:,0],"b")
    #p.plot(np.arange(10,1000,0.5),v[:,3]-30,"r")
    #p.plot(np.arange(10,1000,0.5),v[:,4]-50,"r")
    #p.plot(np.arange(10,1000,0.5),v[:,5]-70,"r")


