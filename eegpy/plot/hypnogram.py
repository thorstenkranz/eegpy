#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import pylab as p
import numpy as np

import eegpy

def plot_hypnogram(et,window=None,stage_map=None,mask_kwargs=None,*args,**kwargs):
    """Plots a hypnogram for EventTable"""
    if stage_map==None:
        stage_map = {"awake":0, "REM":-1, "Stage 1":-2, "Stage 2":-3, "Stage 3":-4, "Stage 4":-5, "mask":0}
    if mask_kwargs==None:
        mask_kwargs = dict(fc="k",alpha=1.0)
    #et = eegpy.EventTable(fn)
    evts = et.get_all_events_with_keys()
    
    if window==None:
        window = np.median(np.diff([x[1] for x in evts]))
    
    xs = np.zeros((len(evts)*2))
    ys = np.zeros((len(evts)*2))
    
    mask_segments = []
    for i,e in enumerate(evts):
        xs[2*i] = e[1]
        if len(evts)-i>1:
            xs[2*i+1] = evts[i+1][1]-1
        else:
            xs[2*i+1] = xs[2*i]+window
        if e[0] in stage_map.keys():
            ys[2*i] = stage_map[e[0]]
            ys[2*i+1] = stage_map[e[0]]                  
        else:
            ys[2*i] = stage_map["mask"]
            ys[2*i+1] = stage_map["mask"]                  
            mask_segments.append((xs[2*i],xs[2*i+1]))
    
    p.plot(xs,ys,*args,**kwargs)
    for ms in mask_segments:
        p.axvspan(ms[0],ms[1],**mask_kwargs)
    p.ylim((np.array(stage_map.values()).min()-1,np.array(stage_map.values()).max()+1))
    #print [x for x in stage_map.keys()],[stage_map[x] for x in stage_map.keys()]
    del stage_map["mask"]
    p.yticks([stage_map[x] for x in stage_map.keys()],[x for x in stage_map.keys()])
    
    return xs,ys

def __main__(fn):
    et = eegpy.EventTable(fn)
    plot_hypnogram(et)
    p.show()

if __name__ == "__main__":
    #__main__(sys.argv[1])
    events = []
    for i in range(20):
        events.append("awake")
        events.append("REM")
        events.append("Stage 1")
        events.append("Stage 2")
        events.append("Stage 3")
        events.append("foo")
    np.random.shuffle(events)
    et = eegpy.EventTable()
    for i,e in enumerate(events):
        et.add_trigger(e,i*20000)
    plot_hypnogram(et)

