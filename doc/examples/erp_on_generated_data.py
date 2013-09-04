# -*- coding: utf-8 -*-
"""Calculate ERP-data on some generated data.

..hint::
    Generation of the data might take some time..."""

import matplotlib.pyplot as plt
from eegpy.example_data.generated import simple_erp_like_data
from eegpy import F32
from eegpy.analysis.timefreq import ERPAnalyzer
from eegpy.temp import UnlinkingTempfile

with UnlinkingTempfile(".f32") as temp_fn:   
    # create some example data
    data, events = simple_erp_like_data()
    eeg = F32(temp_fn, "w+", cNames = ["A", "B", "C"], shape=data.shape, Fs=250.)
    eeg[:,:] = data
    
    # plot whole 'eeg' with events
    plt.subplot(211)
    plt.plot(eeg[:,:], "k-")
    for event in events.all_events:
        plt.axvline(event, color="r", lw=2., alpha=0.5)
    
    # analyze ERPs
    analyzer = ERPAnalyzer(eeg, (-300,+600), baseline=(-200,0))
    for condition in events.keys():
        analyzer.add_condition(condition, events[condition])

    # plot ERPs
    relative_times = analyzer.times_to_trigger
    n_conditions = len(events.keys())
    for i_condition, condition in enumerate(events.keys()):
        ax = plt.subplot(2, n_conditions, i_condition+1+n_conditions)
        plt.text(0, 1, condition, transform=ax.transAxes, va="top")
        for channel_idx, channel_name in enumerate(analyzer.channel_names):
            plt.plot(relative_times, analyzer[condition][:,channel_idx], "k-", 
                     label=channel_name)
        plt.axvline(0, color="r", lw=2., alpha=0.5)
        plt.legend()
        
    eeg.close()
    
    plt.show()



