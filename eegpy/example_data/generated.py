# -*- coding: utf-8 -*-
"""Generating example-data from nonlinear model-systems"""

import numpy as np
from eegpy.models.roessler import NRoessler
from eegpy.events import EventTable

def simple_erp_like_data(t_max = 5000, Fs=2., peek_width=300):
    model = NRoessler(3)
    data = model.integrate(np.arange(0,t_max,1./Fs))
    events = EventTable({"before peeks" : np.arange(500,t_max*Fs-1000,1000),
                         "on peeks" : np.arange(500 + peek_width/2,t_max*Fs-1000,1000)})

    peek_data = np.hanning(peek_width)    
    
    for e in events["before peeks"]:
        data[e:e+peek_width,0] += peek_data*20
        data[e:e+peek_width,3] += peek_data*6
        data[e:e+peek_width,6] += peek_data*-14
    return data[:,0::3], events    
    