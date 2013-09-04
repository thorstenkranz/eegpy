# -*- coding: utf-8 -*-

import os

import numpy as np
from eegpy import F32
from eegpy.temp import UnlinkingTempfile

with UnlinkingTempfile("f32") as filename:
    channel_names = ["F3", "Fz", "F3", "C3" ,"Cz", "C4", "P3", "Pz", "P4"]
    
    eeg = F32(filename, "w+", 
              shape=(10000,len(channel_names)), 
              cNames=channel_names, 
              Fs=250.)
    eeg[:,:] = np.random.random((10000,len(channel_names)))
    
    eeg[:,0] *= 5
    
    print eeg.channel_names
    print eeg.shape
    print eeg.num_datapoints, eeg.num_channels
    print eeg.Fs
    
    eeg.close()
    
    