# -*- coding: utf-8 -*-
"""Calculate Morlet-wavelet analysis on some generated data.

..hint::
    Generation of the data might take some time..."""
import numpy as np
import matplotlib.pyplot as plt
from eegpy.example_data.generated import simple_erp_like_data
from eegpy import F32
from eegpy.analysis.timefreq import WaveletPowerAnalyzer
from eegpy.temp import UnlinkingTempfile

with UnlinkingTempfile(".f32") as temp_fn:   
    # create some example data
    data, events = simple_erp_like_data()
    eeg = F32(temp_fn, "w+", cNames = ["A", "B", "C"], shape=data.shape, Fs=250.)
    eeg[:,:] = data
    
    condition = "before peeks"
    # plot whole 'eeg' with events
    plt.subplot2grid((2,1), (0,0))
    plt.plot(eeg[:,:], "k-")
    for event in events[condition]:
        plt.axvline(event, color="r", lw=2., alpha=0.5)
    
    # analyze ERPs
    analyzer = WaveletPowerAnalyzer(eeg, start_end=(-300,+600), baseline=(-100,0),
                               frequencies = [0.5,1.0,1.5,2.5,5.0,10])
                               
    analyzer.add_condition(condition, events[condition])

    # plot ERPs
    relative_times = analyzer.times_to_trigger
    rts = relative_times
    for channel_idx, channel in enumerate(analyzer.channel_names):
        ax = plt.subplot2grid((4,2), (2+channel_idx%2,channel_idx/2))
        plt.text(0, 1, "channel {}".format(channel), transform=ax.transAxes, va="top",
                 fontsize=12, bbox={"fc":"w", "ec":"none", "alpha":"0.8"})
        plt.imshow(analyzer[condition][:,channel_idx,:].T, 
                   extent = (rts[0], rts[-1], 0, len(analyzer.frequencies)),
                   origin="lower",
                   interpolation="nearest",
                   aspect="auto",
                   vmin=-10,vmax=10)
        plt.yticks(np.arange(0.5, len(analyzer.frequencies)), analyzer.frequencies)
        plt.colorbar()    
        
    eeg.close()
    
    plt.show()
    



