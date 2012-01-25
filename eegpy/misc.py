#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle

debug = True #Contols wether extensive debugging-information is printed to stdout or plotted
show_progressbar = True #Contols wether extensive debugging-information is printed to stdout or plotted

#########################
# FATAL ERROR EXCEPTION #
#########################

class FATALERROR(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
debug = False
#debug = True

class EegpyBase(object):
    def save(self, fn):
        """Saves the object to disk, using the pickle-module"""
        try:
            pickle.dump(self,open(fn,"wb"),2)
        except Exception, e:
            print "Cannot save event-object to disk, ", e  
            
first_pulse_offset = 25
