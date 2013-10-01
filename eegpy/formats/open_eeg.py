#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR, debug
from eegpy.formats.edf import EDFReader
from eegpy.formats.f32 import F32Reader, F32
from eegpy.formats.dbpa import DBPA
#from eegpy.events import EventTable

from eegpy.ui import eegpylab

def open_eeg(fn):
    """Opens the file fn which is expected to contain eeg data
    Determines the type automatically and returns a wrapper"""
    #print "Test open_eeg"
    #TODO: some tests...
    f = None
    if fn.endswith(".f32"):
        try:
            f = F32(fn)
            return f
        except Exception, e:
            print "Couldn't open file as F32-file", e
    elif fn.endswith(".edf") or fn.endswith(".bdf"):
        try:
            f = EDFReader(fn)
            return f
        except IndexError,e:#Exception, e:
            print "Couldn't open file as EDF/BDF-file", e
    elif fn.endswith(".dat"):
        f = DBPA(fn) # JWP: prefers not to use try...except: I want the full error trace
        return f
    else:
        if debug:
            print "Format not recognized"
        raise ValueError("Unknown data format. Couldn'T open file %s"%fn)

