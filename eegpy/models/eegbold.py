#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Model for EEG-BOLD-Coupling via EEG-Power"""

import numpy as np
from scipy.interpolate import splrep,splev
from scipy.signal import hilbert
from scipy.integrate import odeint
from eegpy.filter.freqfilt import filtfilt_band
from eegpy.filter.smoothing import smooth
from roessler import Roessler, NRoesslerTransientCouplings

class RoesslerSuperposTransientCouplings(NRoesslerTransientCouplings):
    """Superposition of signals from two or more Roesslers of diferent frequencies."""
    def __init__(self, omegasOrN = None, es=None, as_=None, bs=None, cs=None, transient_couplings = [], reset_offset = 100, amplitudes=None):
        """Create independent Roessler-systems from the list of arguments
        amplitudes: list of amplitude-factors, eg. [ [0.0,1.0,0.2] , [0.0,1.0,1.0] ] would create two output-signals
        """
        NRoesslerTransientCouplings.__init__(self,omegasOrN,es,as_,bs,cs,transient_couplings,reset_offset)
        
        if amplitudes==None:
            amplitudes = [np.ones((self.N))]
        for i_a in range(len(amplitudes)):
            assert self.N == len(amplitudes[i_a]), "Number of Roesslers and amplitudes must be equal for each amplitude-list"
        self._amplitudes = amplitudes
        print self._amplitudes
        self.v = None #integration-output of parent, not set yet

    def integrate(self,ts):
        self.v = NRoesslerTransientCouplings.integrate(self,ts)
        rv = np.zeros((self.v.shape[0],len(self._amplitudes)))
        for i_a1 in range(len(self._amplitudes)):
            for i_a2 in range(len(self._amplitudes[i_a1])):
                print rv[:,i_a1].shape, (self._amplitudes[i_a1][i_a2]*self.v[:,3*i_a2]).shape
                rv[:,i_a1] += self._amplitudes[i_a1][i_a2]*self.v[:,3*i_a2] #take x-components, sum up
        return rv

class HrfModel_RSTC_LC(object):
    """Model for the neurovascular coupling in the brain.
    Creates EEG-like signal from superposed Roesslers, 
    then models the HRF as response of a damped LC-resonator
    Oscillates if R^2/4L^2 < 1/LC
    if '=', aperiodical limit
    if '>', just decays.
    Must choose something like
      R^2/4L^2 + eps = 1/LC
    with eps small (i.e. strongly damped oscillation).
    """

    #TODO:
    # - smooth in integrate generalisieren

    def __init__(self, input_bands, C=1, L=1, R=None, tau=50, smooth_width=20, RSTC_kwargs=None):
        """Init.
         input_bands: list of 3-tuples
           which frequency bands (fl,fh,factor) add to the input.
         C,L,R: floats
           parameters for the LC-resonator
         tau: int
        """
        self._input_bands = input_bands
        self._C = C
        self._L = L
        if R == None:
            self._R = 2*np.sqrt(L/C) * 0.9 #See class docstring why. Maybe not *0.9, choose different value
        else:
            self._R = R
        self._tau = tau
        self._smooth_width = smooth_width
        if RSTC_kwargs==None: #provide defaults
            print "Make RSTC defaults"
            RSTC_kwargs = dict(
                omegasOrN=np.array([2.7,1.0,2.7]),
                es=np.array([[0.,0.,0.],[0.,0.,0.],[0.1,0.,0.]]), # 1 koppelt in 3 ein
                as_=[0.165,0.165,-0.05],
                transient_couplings=[(300,350,1.)],
                amplitudes=[[0.,0.8,0.2]]
            )
        self._RSTC = RoesslerSuperposTransientCouplings(**RSTC_kwargs)
        self.y = np.random.random((2))*0.2-0.1 # State-vector [Q,I], randomly initalized, close to zero

        #Empty initializations
        self._signal_rstc = None
        self._ts = None
        self._Fs = None
        self._spls = None # Spline-representations of the power-timeseries of the selected bands
    
    def ode(self,y,t):
        """ODE is 
          LQ'' + RQ' + Q/C = 0
        Can be rewritten as two first-order ODEs:
          Q' = I
          I' = -RI - Q/C
        """
        dy = np.zeros((2))
        dy[0] = y[1] + self.get_input(t) # Q' = I + external input, from rstc
        dy[1] = (-self._R*y[1] - y[0]/self._C) / self._L
        return dy

    def integrate(self,ts):
        #First, calculate EEG-timecourse (independently from BOLD)
        self._ts = ts
        self._Fs = 1./(ts[1]-ts[0])
        self._signal_rstc = self._RSTC.integrate(ts)[:,::3]
        self._spls = []
        for in_bd in self._input_bands:
            tmp = filtfilt_band(in_bd[0],in_bd[1],self._signal_rstc[:,0],Fs=self._Fs,border=2)
            power = abs(hilbert(tmp))**2
            smooth_power = smooth(power,int(round(self._smooth_width*self._Fs)))
            self._spls.append(splrep(ts,smooth_power))
        #print "Anzahl spls:", len(self._spls)
        rv = odeint(self.ode,self.y,self._ts)
        return rv

    def get_input(self,t):
        """Calculates the external input from the EEG, i.e. a current"""
        rv = 0
        try:
            for i_b in range(len(self._input_bands)):
                rv += splev(t-self._tau,self._spls[i_b]) * self._input_bands[i_b][2]
        except Exception, e: #time-delay macht probleme
            pass
        return rv
            


##########################
class Hrf1d(HrfModel_RSTC_LC):
    """Creates a chain of HfrModels, connection between the nodes is made be rewiring."""
    def __init__(self, n_nodes, input_bands, rewiring_prob=0.1, C=1, L=1, R=None, tau=50, smooth_width=20, e_scale=0.05, ext_connect=None, ext_scale=0.05,RSTC_kwargs=None):
        self._n_nodes = n_nodes
        self._rewiring_prob = rewiring_prob
        if RSTC_kwargs==None:
            #determine connectivity of base-oscillators according to Watts-Strogartz as small-world network
            connectivity = self.set_connectivity()
            self.connectivity = connectivity
            #Define omegas
            omegas = [2.7] #External oscillator
            for i in range(n_nodes):
                omegas.append(np.random.normal(1.0,0.05)) #Base-osc
                omegas.append(2.7) #damped osc.
            omegas = np.array(omegas)
            #Define es (Couplings)
            self.ext_connect = ext_connect
            self.e_scale = e_scale
            self.ext_scale = ext_scale
            es = self.set_es()
            #Define as
            as_ = [0.165] #External oscillator
            for i in range(n_nodes):
                as_.append(0.165) #Base-osc
                as_.append(-0.05) #damped osc.
            as_ = np.array(as_)
            amplitudes = []
            for i_n in range(n_nodes):
                #amplitudes.append([0]+i_n*[0.,0.]+[0.6,0.4]+(n_nodes-i_n-1)*[0.,0.])
                amplitudes.append([0]+i_n*[0.,0.]+[1.0,1.4]+(n_nodes-i_n-1)*[0.,0.])
            #Make RSTC_kwargs
            RSTC_kwargs = dict(
                omegasOrN=omegas,
                es=es,
                as_=as_,
                transient_couplings=[(200,250,1.),(400,450,0.7),(600,650,1.3),(800,850,1.6),(1000,1050,1.0),(1200,1250,1.),(1400,1450,0.7),(1600,1650,1.3),(1800,1850,1.6),(2000,2050,1.0),(2200,2250,1.),(2400,2450,0.7),(2600,2650,1.3),(2800,2850,1.6),(3000,3050,1.0),(3200,3250,1.),(3400,3450,0.7),(3600,3650,1.3),(3800,3850,1.6),(4000,4050,1.0),(4200,4250,1.),(4400,4450,0.7),(4600,4650,1.3),(4800,4850,1.6),(5000,5050,1.0),],
                amplitudes=amplitudes
            )

            
        super(Hrf1d,self).__init__( input_bands, C=C, L=L, R=R, tau=tau, smooth_width=smooth_width, RSTC_kwargs=RSTC_kwargs)
        self.y = np.random.random((2*n_nodes))*0.2-0.1 # State-vector [Q,I], randomly initalized, close to zero

    def set_connectivity(self):
        connectivity = np.zeros((self.n_nodes,self.n_nodes))
        for i in range(self.n_nodes-1):
            if np.random.random()<self._rewiring_prob:
                idx = i
                while idx==i:
                    idx = np.random.randint(0,self.n_nodes)
                    connectivity[i,idx] = connectivity[idx,i] = 1
            else:
                connectivity[i,i+1] = connectivity[i+1,i] = 1
        return connectivity

    def set_es(self):
        if self.ext_connect == None:
            self.ext_connect = np.zeros((self.n_nodes))
            self.ext_connect[self.n_nodes/2]=1
        print "ext_connect:\n", self.ext_connect
        #Define es (Couplings)
        es = np.zeros((2*self.n_nodes+1,2*self.n_nodes+1))
        for i in range(self.n_nodes):
            if self.ext_connect[i]:
                es[2*i+2,0] = (self.ext_scale+np.random.random()*self.ext_scale)
            for j in range(i):
                if self.connectivity[i,j]:
                    es[2*i+1,2*j+1] = es[2*j+1,2*i+1] = np.random.random()*self.e_scale+self.e_scale
                    es[2*i+2,2*j+2] = es[2*j+2,2*i+2] = np.random.random()*self.e_scale/10+self.e_scale/10
        self.es = es
        return es

    def ode(self,y,t):
        """ODE is 
          LQ'' + RQ' + Q/C = 0
        Can be rewritten as two first-order ODEs:
          Q' = I
          I' = -RI - Q/C
        """
        n_nodes = self._n_nodes
        dy = np.zeros((2*n_nodes))
        for i_n in range(n_nodes):
            dy[2*i_n+0] = y[2*i_n+1] + self.get_input(t,i_n) # Q' = I + external input, from rstc
            dy[2*i_n+1] = (-self._R*y[2*i_n+1] - y[2*i_n+0]/self._C) / self._L
        return dy

    def integrate(self,ts):
        #First, calculate EEG-timecourse (independently from BOLD)
        self._ts = ts
        self._Fs = 1./(ts[1]-ts[0])
        self._signal_rstc = self._RSTC.integrate(ts)[:,:]
        print self._signal_rstc.shape
        self._spls = np.zeros((len(self._input_bands),self._n_nodes),"O")
        for i_b, in_bd in enumerate(self._input_bands):
            for i_n in range(self._n_nodes):
                tmp = filtfilt_band(in_bd[0],in_bd[1],self._signal_rstc[:,i_n],Fs=self._Fs,border=2)
                power = abs(hilbert(tmp))**2
                smooth_power = smooth(power,int(round(self._smooth_width*self._Fs)))
                self._spls[i_b,i_n] = splrep(ts,smooth_power)
        #print "Anzahl spls:", len(self._spls)
        rv = odeint(self.ode,self.y,self._ts)
        return rv

    def get_input(self,t,i_n):
        """Calculates the external input from the EEG, i.e. a current
        For time t and node i_n."""
        rv = 0
        try:
            for i_b in range(len(self._input_bands)):
                rv += splev(t-self._tau,self._spls[i_b,i_n]) * self._input_bands[i_b][2]
        except Exception, e: #time-delay macht probleme
            pass
        return rv

    @property
    def n_nodes(self):
        return self._n_nodes


class Hrf2d(Hrf1d):
    """Organize oscillators as a 2d grid, making (not-rewired) connectivity depend on neighborhood"""
    def __init__(self, n_nodes, positions, input_bands, conn_dist=1.5, rewiring_prob=0.1, C=1, L=1, R=None, tau=50, smooth_width=20, e_scale=0.05, ext_connect=None, ext_scale=0.05,RSTC_kwargs=None):
        """
        conn_dist: max distance for connectivity if not rewired
        """
        self.positions=positions
        self.conn_dist = conn_dist
        Hrf1d.__init__(self,n_nodes,input_bands,rewiring_prob, C, L, R, tau, smooth_width, e_scale, ext_connect, ext_scale, RSTC_kwargs)

    def set_connectivity(self):
        connectivity = np.zeros((self.n_nodes,self.n_nodes))
        for i in range(self.n_nodes):
            for j in range(self.n_nodes):
                if np.sqrt(np.sum((self.positions[i]-self.positions[j])**2)) <self.conn_dist:
                    if np.random.random()<self._rewiring_prob:
                        #rewiring
                        idx = i
                        while idx==i:
                            idx = np.random.randint(0,self.n_nodes)
                            connectivity[i,idx] = connectivity[idx,i] = 1
                    else:
                        connectivity[i,j] = connectivity[j,i] = 1
        return connectivity




if __name__=="__main__":
    import pylab as p
    from eegpy.analysis.wavelet import wt_power_baseline
    from scipy.stats import scoreatpercentile

    ##Test 1: Überlagerte Rössler mit transienter Kopplung
    #roe = RoesslerSuperposTransientCouplings(
    #    omegasOrN=np.array([2.7,1.0,2.7]),
    #    es=np.array([[0.,0.,0.],[0.,0.,0.],[0.1,0.,0.]]), # 1 koppelt in 3 ein
    #    as_=[0.165,0.165,-0.05],
    #    transient_couplings=[(300,350,1.)],
    #    amplitudes=[[0.,0.8,0.2],[0.,1.0,1.0]])
    #ar = roe.integrate(np.arange(10,410,0.1))

    #p.subplot(2,2,1)
    #p.plot(np.arange(10,410,0.1),ar[:,0])
    #p.xlim(10,409)
    #p.subplot(2,2,3)
    #wt_power = wt_power_baseline(ar[:,0],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(2500,2800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    ##print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="nearest",extent=[9.975,409.5,0.05,1.475],vmin=vmin,vmax=vmax)
    #p.colorbar()
    #p.xlim(10,409)
    #p.subplot(2,2,2)
    #p.plot(np.arange(10,410,0.1),ar[:,1])
    #p.xlim(10,409)
    #p.subplot(2,2,4)
    #wt_power = wt_power_baseline(ar[:,1],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(2500,2800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    ##print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="nearest",extent=[9.975,409.5,0.05,1.475],vmin=vmin,vmax=vmax)
    #p.colorbar()
    #p.xlim(10,409)

    #Test 2: HRF-Modell
    #hrfm = HrfModel_RSTC_LC(input_bands=[(0.5,1.0,1.0)],C=1,L=1,tau=50,
    #        RSTC_kwargs = dict(
    #            omegasOrN=np.array([2.7,1.0,2.7]),
    #            es=np.array([[0.,0.,0.],[0.,0.,0.],[0.1,0.,0.]]), # 1 koppelt in 3 ein
    #            as_=[0.165,0.165,-0.05],
    #            transient_couplings=[(100,150,1.)],
    #            amplitudes=[[0.,0.4,0.6]])
    #        )
    #ar = hrfm.integrate(np.arange(10,410,0.1))[:,0]

    #p.subplot(2,1,1)
    #p.plot(np.arange(10,1210,0.1),ar)
    #p.xlim(10,1209)
    #p.subplot(3,1,1)
    #p.plot(np.arange(10,410,0.1),hrfm._signal_rstc[:])
    #p.xlim(10,409)
    #p.subplot(3,1,2)
    #p.plot(np.arange(10,410,0.1),ar)
    #p.xlim(10,409)
    #p.subplot(3,1,3)
    #wt_power = wt_power_baseline(hrfm._signal_rstc[:,0],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(1500,1800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    ##print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="nearest",extent=[9.975,409.5,0.05,1.475],vmin=vmin,vmax=vmax)
    ##p.colorbar()
    #p.xlim(10,409)

    #Test 3: HRF-Modell-Chain
    #nnodes=6#20
    #hrfm = Hrf1d(n_nodes=nnodes,input_bands=[(0.5,1.0,3.0)],C=1,L=1,tau=50)
    #ar = hrfm.integrate(np.arange(10,510,0.1))#[:,nnodes/2]

    #p.subplot(3,1,1)
    #pd1 = hrfm._signal_rstc[:]+np.arange(0,nnodes*15,15).reshape(1,-1).repeat(hrfm._signal_rstc.shape[0],axis=0)
    #p.plot(np.arange(10,510,0.1),pd1)
    #p.xlim(10,509)
    #p.ylim(pd1[500:].min(),pd1[500:].max())
    #p.subplot(3,1,2)
    #p.plot(np.arange(10,510,0.1),ar[:,::2])#nnodes/2*2])#+np.arange(0,nnodes*2*15,15).reshape(1,-1).repeat(ar.shape[0],axis=0))
    #p.xlim(10,509)
    #p.ylim(ar[500:,::2].min(),ar[500:,::2].max())
    #p.subplot(3,1,3)
    #wt_power = wt_power_baseline(hrfm._signal_rstc[:,nnodes/2],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(1500,1800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    #print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="bicubic",extent=[9.95,509.5,0.05,1.475],vmin=vmin,vmax=vmax)
    #p.colorbar()
    #p.xlim(10,509)
    
    #Test 4: HRF-Modell-Grid
    #nnodes=49#20
    #positions = np.zeros((49,2))
    #for i in range(49):
    #    x=i/7
    #    y=i%7
    #    positions[i,0] = x*1
    #    positions[i,1] = y*1
    #hrfm = Hrf2d(n_nodes=nnodes,positions=positions,input_bands=[(0.5,1.0,1.0)],C=1,L=1,tau=50)
    #ar = hrfm.integrate(np.arange(10,1010,0.1))#[:,nnodes/2]

    #p.subplot(3,1,1)
    #pd1 = hrfm._signal_rstc[:]+np.arange(0,nnodes*15,15).reshape(1,-1).repeat(hrfm._signal_rstc.shape[0],axis=0)
    #p.plot(np.arange(10,1010,0.1),pd1)
    #p.xlim(10,1009)
    #p.ylim(pd1[500:].min(),pd1[500:].max())
    #p.subplot(3,1,2)
    #p.plot(np.arange(10,1010,0.1),ar[:,::2])#nnodes/2*2])#+np.arange(0,nnodes*2*15,15).reshape(1,-1).repeat(ar.shape[0],axis=0))
    #p.xlim(10,1009)
    #p.ylim(ar[500:,::2].min(),ar[500:,::2].max())
    #p.subplot(3,1,3)
    #wt_power = wt_power_baseline(hrfm._signal_rstc[:,nnodes/2],freqs=np.arange(0.1,1.5,0.05),Fs=10.,baseline=slice(1500,1800),wtSampleWidth=20)
    #vmin = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),5)
    #vmax = scoreatpercentile(wt_power[ar.shape[0]/4*1:ar.shape[0]/4*3].flatten(),95)
    #print vmin,vmax,ar.shape[0]/4*1,ar.shape[0]/4*3
    #p.imshow(wt_power.T,aspect="auto",origin="lower",interpolation="bicubic",extent=[9.95,1009.5,0.05,1.475],vmin=vmin,vmax=vmax)
    #p.colorbar()
    #p.xlim(10,1009)
    
