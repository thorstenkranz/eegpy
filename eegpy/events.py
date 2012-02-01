#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Classes and funtions dealing with various kinds of event-informations."""

from __future__ import with_statement
import os.path
import pickle
import copy

import eegpy
from eegpy.misc import FATALERROR, EegpyBase, debug

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

class EventTable(EegpyBase):
    """Contains a list of events. 
    
    Can be constructed from various file-types or manually."""
    
    #The next three variables define the events. In the first one, the original 
    #eventmarks are saved. The other two are used for an affine transformation, 
    #such that    event = (ur_event+_offset) * _factor
    _ur_eventdict = {}
    _lambdas = [] #list of functions that are applied to the timepoints in transform
    _inverse_lambdas = [] #list of functions that are applied to the timepoints in inverse-transform
    #_offset = 0
    #_factor = 1
    
    def __init__(self,x=None):
        self._ur_eventdict = {}
        self._lambdas = []
        self._inverse_lambdas = []
        #self._offset = 0
        #self._factor = 1 
        success = False
        if debug:
            print "EventTable.__init__"
        if x==None:
            success=True
        if not success:
            try:
                if type(x) == type(""):
                    #String
                    tmp = eegpy.load(x)
                    self._ur_eventdict = copy.deepcopy(tmp._ur_eventdict)
                    #try:
                    self._lambdas = tmp._lambdas
                    self._inverse_lambdas = tmp._inverse_lambdas
                    #except Exception, e:
                    #    if debug:
                    #        print "Error in EventTable.__init__ while setting lambdas"
                    success=True
            except Exception,e:
                pass
                #if debug:
                #    print "First try: Loading file as saved EventTable didn't work.", e
        if not success:
            try:
                if len(x.keys()) > 0:
                    #TODO: Check if structure of dict is o.k.
                    for k in x.keys():
                        it = iter(x[k]) #If x[k] is not iterable, it raises an error
                    self._ur_eventdict = x
                    success = True
            except Exception, e:
                pass
                #if debug:
                #    print "Second try: Open x as dictionary of iterables didn't work.", e
        if not success:
            try:
                if type(x) == type(""):
                    #String
                    self.read_from_bv(x)
                    if len(self.keys()) == 0:
                        success=False
                        raise ValueError("No BVA-file given")
                    else:
                        success=True
            except Exception,e:
                if debug:
                    print "Third try: Reading file as BrainVision Analyzer didn't work.", e
        if not success:
            raise ValueError("Construction of EventTable from given variable failed")            
        
        
        #TODO: importing from arrays or different filetypes
        
        
        all_events = property(self.get_all_events)
    
    
### Transformations-stuff    
    #numerical operators will set _offset and _factor
    def __add__(self,offset):
        """This method adjusts the value of _offset."""
        try:
            if abs(int(round(offset))-offset)<1:
                offset = int(round(offset))
                self._lambdas.append("lambda x: x+%i"%offset) #_offset += int(round(offset))
                self._inverse_lambdas.insert(0,"lambda x: x-%i"%offset)
            return self
        except Exception, e:
            if debug:
                print "Error while setting the offset", e
            raise NotImplemented()
    
    def __radd__(self,offset):
        return self.__add__(offset)
    
    def __sub__(self,offset):
        return self.__add__(-offset)
    
    def __rsub__(self,offset):
        if debug:
            print "Operations of the kind 1-EventTable don't make sense and are therefore not implemented"
        return NotImplemented#return self.__add__(-offset)
    
    def __mul__(self,factor):
        #self._factor *= factor
        assert factor!=0, "Multiplication with 0 doesn't make sense and is therefore not implemented"
        self._lambdas.append("lambda x: x*%i"%factor)
        self._inverse_lambdas.insert(0,"lambda x: float(x)/%i"%factor)
        return self
    
    def __rmul(self,factor):
        return self.__mul__(factor)
    
    def __div__(self,factor):
        self.__mul__(1.0/factor)
        return self
    
    def __rdiv__(self,offset):
        if debug:
            print "Operations of the kind 1/EventTable don't make sense and are therefore not implemented"
        return NotImplemented#return self.__add__(-offset)
    
    def _transform(self,x):
        x = copy.copy(x)
        for g in [eval(s) for s in self._lambdas]:
            x = g(x)
        return x
    
    def _invert_transform(self, x):
        x = copy.copy(x)
        for g in [eval(s) for s in self._inverse_lambdas]:
            x = g(x)
        return x
    
    def reset_transformation(self):
        #self._offset = 0
        #self._factor=1
        self._lambdas = []
        self._inverse_lambdas = []

### Adding new triggers    
    def add_trigger_type(self, key, ts):
        """Manually add list of triggers of type key at times ts to EventTable.
        The timepoints must be provided in the present time-transformation and 
        are back-transformed included in the _ur_eventdict"""
        if self._ur_eventdict.has_key(key):
            #if debug:
            print "Trigger-type %s is being replaced." % str(key)
        self._ur_eventdict[key] = [self._invert_transform(x) for x in ts]
        
    def add_trigger(self,key,t):
        """Manually add a triggers of type key at time t to EventTable.
        The timepoint must be provided in the present time-transformation and 
        are back-transformed included in the _ur_eventdict"""
        if not type(t) == int:
            raise ValueError("t must be integer")
        if not self._ur_eventdict.has_key(key):
            self._ur_eventdict[key] = []
        #if self._ur_eventdict.has_key(key):
            #if debug:
            #print "Trigger-type %s is being replaced." % str(key)
        self._ur_eventdict[key].append(self._invert_transform(t))
        #else:
            #raise ValueError("%s is not a valid trigger-type" %key)
            
    #Removal 
    def remove(self,k,tr=None):
        if tr==None:
            #remove a key
            if self._ur_eventdict.has_key(k):
                del self._ur_eventdict[k]
            else:
                raise KeyError("Key %s doesn't exist." % k)
        else:
            #remove one trigger:
            if self._ur_eventdict.has_key(k):
                for i in range(len(self._ur_eventdict[k])):
                    if self._transform(self._ur_eventdict[k][i])==tr:
                        del self._ur_eventdict[k][i]
                        return True
                raise ValueError("Trigger %s doesn't exist in key %s." % (tr,k))
            else:
                raise KeyError("Key %s doesn't exist." % k)    

    #Automatic removal of tps < 0
    def remove_negative(self):
        for k in self._ur_eventdict.keys():
            for tr in self._ur_eventdict[k]:
                if self._transform(tr)<0:
                    self.remove(k,self._transform(tr))
                

    def __setitem__(self,key,ts):
        self.add_trigger_type(key, ts)
        
### Access to elements    
    #Item-wise access provided via __get_item__-Method. It returns the set of 
    #triggers refered to by key in ur_eventdict correctly transformed
    def __getitem__(self,key):
        rv = [int(round(self._transform(y))) for y in self._ur_eventdict[key]]
        rv.sort()
        return rv
    
    def events_by_names(self,*args):
        rv = []
        try:
            if len(args)<1:
                raise ValueError()
            for key in args:
                #print key
                try:
                    rv += [self._transform(y) for y in self._ur_eventdict[key]]
                except KeyError, e:
                    if debug:
                        print "Key %s not in EventTable. Ignored." %key
            rv.sort()
            return rv
        except Exception, e:
            if debug:
                print "Error while getting items", e
       
    def events_by_names_with_keys(self,*args):
        rv = []
        try:
            if len(args)<1:
                raise ValueError()
            for key in args:
                #print key
                try:
                    rv += [(key,self._transform(y)) for y in self._ur_eventdict[key]]
                except KeyError, e:
                    if debug:
                        print "Key %s not in EventTable. Ignored." %key
            rv.sort(cmp=lambda x,y:int(x[1]-y[1]))
            return rv
        except Exception, e:
            if debug:
                print "Error while getting items", e
       
    def get_all_events(self):
        """Gets all event timepoints, sorted by time"""
        all_events = []
        for k in self._ur_eventdict.keys():
            all_events += [self._transform(y) for y in self._ur_eventdict[k]]
        all_events.sort()
        return all_events
    
    def get_all_events_with_keys(self):
        """Gets all event timepoints together with keys, sorted by time"""
        all_events = []
        for k in self._ur_eventdict.keys():
            all_events += [(k,self._transform(y)) for y in self._ur_eventdict[k]]
        all_events.sort(cmp=lambda x,y:int(x[1]-y[1]))
        return all_events
    
    def keys(self):
        return self._ur_eventdict.keys()
    
    def has_key(self,k):
        return self._ur_eventdict.has_key(k)
        
### Import-methods   
    def read_from_bv(self, x):
        """Reads a file as exported by BrainVision Analyzer
        x is checked to be a string"""
        #Check if file exists
        inMarkerInfo = False
        try:
            with open(x) as f:
                for line in f:
                    #if debug:
                    #    print line
                    if line.startswith("["):
                        inMarkerInfo = False
                    if line.startswith("[Marker Infos]"):
                        inMarkerInfo = True
                    if inMarkerInfo:
                        if line.startswith("Mk"):
                            #print line
                            ldAsList = line.split("=")[1].split(",")
                            #if ldAsList[0]=="Trigger": #take only events of type Trigger
                            description = ldAsList[1]
                            pos = int(ldAsList[2])
                            if not self._ur_eventdict.has_key(description):
                                self._ur_eventdict[description] = []
                            self._ur_eventdict[description].append(pos)            
        except IOError, ioe:
            print "File %s could not be found for importing events." % x  
    
### Helper-methods
    def insert_missing_triggers(self,key):
        """For equaly spaced triggers, detect missing ones and insert them
        Most frequently used for the correction of scanpulses"""
        trigs = list(self[key])
        med = n.median(n.diff(n.array(trigs)))
        start_len = len(trigs)
        for i in range(start_len-1):
            step = float(trigs[i+1]-trigs[i])/med
            error = abs(step-round(step))
            num_dist = int(round(step))
            if error>0.05:
                raise ValueError("The triggers cannot be corrected, more than 5 percent error encountered. num_dist, error = %i, %f" %(num_dist,error))
            else:
                for j in range(1,num_dist):
                    if debug:
                        print "Inserted trigger at ", int(round(trigs[i]+j*med))
                    trigs.append(int(round(trigs[i]+j*med)))
        trigs.sort()
        self[key] = trigs
    
    def find_closest_event(self, t, keys=None):
        """Find the event that is temporaly closest to the given time t.
        Restrict search to given keys."""
        if keys==None:
            keys = self.keys()
        mi = 10e50
        best = None
        best_key = None
        for k in keys:
            if self.has_key(k):
                for i,v in enumerate(self[k]):
                    if abs(t-v)<mi:
                        mi=abs(t-v)
                        best = v
                        best_key = k
                    #else: 
                    #    break
        return best, best_key
    
    def find_next_event(self,t,keys=None):
        """Like find_closest_event (see there for details), just that only 
        in positive direction of t is searched." 
        """
        if keys==None:
            keys = self.keys()
        mi = 10e50
        best = None
        best_key = None
        for k in keys:
            if self.has_key(k):
                for i,v in enumerate(self[k]):
                    if v-t<mi and v-t>0:
                        mi=v-t
                        best = v
                        best_key = k
                    #else: 
                    #    break
        return best, best_key
            
    #Two Methods for pickling
    def __getstate__(self):
        odict = self.__dict__
        #odict["_ur_eventdict"] = pickle.dumps(self._ur_eventdict)
        return odict
    
    def __setstate__(self,dict):
        #self.__dict__ = dict
        self._ur_eventdict = copy.deepcopy(dict["_ur_eventdict"])
        if dict.has_key("_lambdas"):
            self._lambdas = dict["_lambdas"]
            self._inverse_lambdas = dict["_inverse_lambdas"]
        else:
            self._lambdas = []
            self._inverse_lambdas = []
            #print dict.keys()
            if debug:
                print "EventTable: Converting from old format"
            self += dict["_offset"]
            self *= dict["_factor"]
            
    
def judge_rejection(et0,et_ref,et1, keys=None ):
    """To judge quality of artifact rejection algorithms.

    :Parameters:
      et0: EventTable
        Original EventTable
      et_ref: EventTable
        EventTable to compare against
      et1: EventTable
        EventTable to be judged
      keys: List
        Which keys should be taken into account. Standard: all in et0

    Returns:
      (accuracy, precision, specificity, sensitivity)
    """
    if keys==None:
        keys=et0.keys()
    events0 = et0.events_by_names_with_keys(*keys)
    events_ref = et_ref.events_by_names_with_keys(*keys)
    events1 = et1.events_by_names_with_keys(*keys)

    rejections_ref = set(events0).difference(set(events_ref))
    remains_ref = set(events0).intersection(set(events_ref))
    
    rejections1 = set(events0).difference(set(events1))
    remains1 = set(events0).intersection(set(events1))
    
    correct_alarms = rejections1.intersection(rejections_ref) #triggers that were rejected correctly
    false_alarms = rejections1.intersection(remains_ref) #triggers that were rejected by error
    correct_misses = remains1.intersection(remains_ref) #triggers that were left in EventTable correctly
    false_misses = remains1.intersection(rejections_ref) #triggers that were left in EventTable by error

    #print len(correct_alarms), len(false_alarms), len(correct_misses), len(false_misses)
    
    # Definitions from Wikipedia
    try:
        accuracy = float(len(correct_alarms)+len(correct_misses))/(len(correct_alarms)+len(false_alarms)+len(correct_misses)+len(false_misses))
    except ZeroDivisionError, zde:
        accuracy=n.NaN
    try:
        precision= float(len(correct_alarms))/(len(correct_alarms)+len(false_alarms))
    except ZeroDivisionError, zde:
        precision=n.NaN
    try:
        specificity = float(len(correct_misses))/(len(correct_misses)+len(false_alarms))
    except ZeroDivisionError, zde:
        specificity=n.NaN
    try:
        sensitivity = float(len(correct_alarms))/(len(correct_alarms)+len(false_misses))
    except ZeroDivisionError, zde:
        sensitivity=n.NaN
    
    return (accuracy, precision, specificity, sensitivity)

