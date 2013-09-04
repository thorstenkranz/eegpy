#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: eegpy.misc
   :platform: Unix, Windows, OSX
   :synopsis: A useful module indeed.

.. moduleauthor:: Thorsten Kranz <thorstenkranz@gmail.com>


"""

import pickle

debug = False #Contols wether extensive debugging-information is printed to stdout or plotted
show_progressbar = True #Show progess-bar for long running tasks.

class EegpyBase(object):
    """Base class of all eegpy-objects.

    Adds some convenience features to our objects, e.g. for persistence.

    .. note::

      Uses the pickle module
    """

    def save(self, fn):
        """Saves the object to disk, using the pickle-module"""
        try:
            pickle.dump(self,open(fn,"wb"),2)
        except Exception, e:
            print "Cannot save eegpy-object to disk, ", e  
            

#########################
# FATAL ERROR EXCEPTION #
#########################

#Not so clever...
class FATALERROR(Exception):
    """Fatal error, indicating severe problems and suggesting abortion"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

first_pulse_offset = 25
