#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for working with sensor-locations.
Will support reading of .ced and .loc files as eeglab uses.
"""
#################
# Module-Import #
#################
import os.path

#eegpy-modules
import eegpy
from eegpy.misc import FATALERROR, debug, EegpyBase

#Third-party
try:
    import numpy as N
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

standard_ced = os.path.join(os.path.split(__file__)[0],"sample_locs/Standard-10-20-Cap-Extended.ced")

class _SensorCed(EegpyBase):
    def __init__(self,*args):
        EegpyBase.__init__(self)
        if len(args) == 0:
            #No value given, return Cz
            self.set_values(0, 'Cz', 90.0, 0.0, 3.75e-33, -6.12e-17, 1.0, -90.0, 90.0, 1.0)
        elif len(args)==1:
            #String passed, as one line of .ced-file
            splt = args[0].split()
            splt[0] = int(splt[0])
            splt[2:] = [float(s) for s in splt[2:]]
            self.set_values(*splt)
        elif len(args)==10:
            self.set_values(*args)
        else:
            raise ValueError("Illegal input for creating a sensor.")
    
    def __repr__(self):
        return """_SensorCed object.
    idx: %i
    name: %s 
    polar2d: %s 
    cartesian3d: %s 
    polar3d: %s""" % (self.idx,self.name,str(self.polar2d),str(self.cartesian3d),str(self.polar3d))
        
    def set_values(self,idx,name,pol2d1,pol2d2,cart3d1,cart3d2,cart3d3,pol3d1,pol3d2,pol3d3):
            self.idx = idx
            self.name=name
            self.polar2d=N.array([N.deg2rad(pol2d1), pol2d2])
            self.cartesian3d=N.array([cart3d1, cart3d2, cart3d3])
            self.polar3d=N.array([N.deg2rad(pol3d1), N.deg2rad(pol3d2), pol3d3])
            
class SensorLocations(EegpyBase):
    def __init__(self,fn=None):
        EegpyBase.__init__(self)
        self._sensors = []
        self._fn = fn
        if self._fn != None:
            try:        
                self.read_ced(self._fn)
            except Exception, e:
                print "Error reading .ced-file:", e
        else:
            try:        
                self.read_ced(standard_ced)
            except Exception, e:
                print "Error reading standard .ced-file:", e
    
    def __repr__(self):
        return "SensorLocations object containing %i sensors" % len(self._sensors)
                
    def read_ced(self,fn):
        f = open(fn,"r")
        for i,line in enumerate(f.readlines()):
            if i>0:
                self._sensors.append(_SensorCed(line))
    
    def get_all_polar2d(self):
        rv = N.zeros((len(self._sensors),2),"d")
        for i,sens in enumerate(self._sensors):
            rv[i,:] = sens.polar2d
        return rv
    
    def get_all_cartesian3d(self):
        rv = N.zeros((len(self._sensors),3),"d")
        for i,sens in enumerate(self._sensors):
            rv[i,:] = sens.cartesian3d
        return rv
    
    def get_all_polar3d(self):
        rv = N.zeros((len(self._sensors),3),"d")
        for i,sens in enumerate(self._sensors):
            rv[i,:] = sens.polar3d
        return rv
    
    def __getitem__(self,key):
        if type(key) == str:
            for i in range(len(self._sensors)):
                if self._sensors[i].name.lower() == key.lower():
                    return self._sensors[i]
            raise ValueError("Sensor %s not found"%key)#return None
        elif type(key) == int:
            for i in range(len(self._sensors)):
                if self._sensors[i].idx == key:
                    return self._sensors[i]
            raise ValueError("Sensor No. %i not found"%key)#return None
        else:
            raise ValueError("Key must be of string or integer.")#return None
        
    def get_for_names_cartesian3d(self,names,catch_errors=True):
        rv = N.zeros((len(names),3),"d")
        for i in range(rv.shape[0]):
            try:
                rv[i,:] = self[names[i]].cartesian3d
            except Exception,e:
                if catch_errors:
                    rv[i,:] = [0.0,0.0,0.0]
                else:
                    raise e
        return rv
    
    def get_channel_boolslice(self,names):
        #for given names, check if positions can be found and return a boolslice for further usage
        rv = N.zeros((len(names)),N.bool)
        for i in range(rv.shape[0]):
            try:
                self[names[i]].cartesian3d
                rv[i] = True
            except Exception, e:
                pass
        return rv
        
                
        
if __name__=="__main__":
    print "module sensorlocations executing"
    sloc = SensorLocations("sample_locs/Standard-10-20-Cap-Extended.ced")
    
