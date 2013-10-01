# -*- coding: utf-8 -*-
u"""
Create some artificial data using coupled Rößler-oszillators.
"""
import numpy as np
import matplotlib.pyplot as plt
from eegpy.models.roessler import TwoRoessler

model_system = TwoRoessler(e1 = 0.03, e2=0.03)
data1 = model_system.integrate(np.arange(0,500,0.1))


from eegpy.analysis.phases import phase_coherence
print phase_coherence(data1[:, 0], data1[:,3])

model_system.e1 = 0.07
model_system.e2 = 0.07
data2 = model_system.integrate(np.arange(0,500,0.1))
print phase_coherence(data2[:, 0], data2[:,3])   
    
data = np.concatenate([data1, data2],0)
for i in range(3):
    plt.subplot(3,1,i+1)
    plt.plot(data[:,i], "k-")
    plt.plot(data[:,i+3], "b-")
    
plt.show()
