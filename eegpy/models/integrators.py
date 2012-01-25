#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def eulerint(ode,y,ts,h):
    """Integrate an ODE using the Euler-Method""" 
    assert abs(np.diff(ts).std())<=10**-10, "ts must be equally spaced"
    assert abs(round((ts[1]-ts[0])/h)-(ts[1]-ts[0])/h)<=10**-10, "step between ts must be a multiple of h"
    assert abs(round(ts[0]/h)-ts[0]/h)<=10**-10, "all ts must be a multiple of h"
    rv = np.zeros((ts.shape[0],y.shape[0]))
    t = 0
    y_tmp = y.copy()
    for i_t, next_t in enumerate(ts):
        print i_t, next_t
        while t<next_t:
            if t%10==0:
                print t
            dydt = ode(y_tmp,t)
            y_tmp += dydt*h
            t+=h
        rv[i_t,:] = y_tmp[:]
    return rv


def rk4int(ode,y,ts,h=None):
    """Integrate an ODE using 4th-order Runge-Kutta-Method""" 
    def rk4(y, dydt, t, h):
        i=0
        th=0.0
        hh=0.0
        h6=0.0
    
        n=len(y)
        
        hh=h*0.5
        h6=h/6.0
        th=t+hh
        yt=y+hh*dydt
        dyt = ode(yt,th)
        yt=y+hh*dyt
        dym = ode(yt,th)
        yt=y+h*dym
        dym += dyt
        dyt = ode(yt,t+h)
        yout=y+h6*(dydt+dyt+2.0*dym)
        return yout

    if h==None:
        h=ts[1]-ts[0]
    assert abs(np.diff(ts).std())<=10**-7, "ts must be equally spaced"
    assert abs(round((ts[1]-ts[0])/h)-(ts[1]-ts[0])/h)<=10**-10, "step between ts must be a multiple of h"
    assert abs(round(ts[0]/h)-ts[0]/h)<=10**7, "all ts must be a multiple of h"


    rv = np.zeros((ts.shape[0],y.shape[0]))
    t = 0
    y_tmp = y.copy()
    for i_t, next_t in enumerate(ts):
        #print i_t, next_t
        while t<next_t:
            #if abs(t%10)<10**-10:
            #    print t
            dydt = ode(y_tmp,t)
            y_tmp = rk4(y_tmp, dydt, t, h)
            t+=h
        rv[i_t,:] = y_tmp[:]
    return rv






