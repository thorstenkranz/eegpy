#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Transfer function systems"""

import numpy as np
from roessler import TwoStochasticRoessler

class Winterhalder(object):
    """Implements a autoregressive (thus lowpass-filtered) and time-delayed 
    propagated signal. It uses a system of two diffusively coupled stochastic
    Roessler oscillators.
    cp. Winterhalder et al., PLA A 356 (2006), 26-34
    Equation:
    u(t) = a_1 u(t) + a_2 x(t-tau) + sigma n(t)
    """
    def __init__(self,a1=0.3,a2=0.7,tau=64,sigma=1.5,roessler_kwargs={},return_input=False):
        self._a1 = a1
        self._a2 = a2
        self._tau = tau
        self._sigma = sigma
        self._roessler_kwargs = roessler_kwargs
        self._return_input = return_input #wether also the Roessler-timeseries should be returned

        self._tsr = TwoStochasticRoessler(**roessler_kwargs)
    
    def integrate(self,ts):
        """Integrate the underlying system and propagate and dilate it."""
        ts = np.array(ts)
        timestep = ts[1]-ts[0] # Assume equally spaced times, needed for TwoStochasticRoessler anyways
        ts += 100
        ts = np.concatenate([np.arange(ts[0]-1000*timestep,ts[0],timestep),ts]) #iteriere 500 Schritte vor
        print ts, np.diff(ts).min(),np.diff(ts).max(),np.median(np.diff(ts)),ts[998:1002]
        roe_ts = self._tsr.integrate(ts)

        #Define Shorthands
        a1 = self._a1
        a2 = self._a2
        tau = self._tau
        sigma = self._sigma

        rv = roe_ts[:,0].copy()
        for i in range(tau,roe_ts.shape[0]):
            if i%100==0:
                print i
            rv[i] = a1*rv[i-1] + a2*roe_ts[i-tau,0] + sigma*np.random.normal(0,1)

        if self._return_input:
            return rv[1000:], roe_ts[1000:,:]
        else:
            return rv[1000:]


            

if __name__ == "__main__":
    winter = Winterhalder(a1=0.9,a2=0.1,tau=100,roessler_kwargs=dict(omega1=1+0.05,omega2=1-0.05),return_input=True)
    wh_ts, roe_ts = winter.integrate(np.arange(100,509.6,0.01))

