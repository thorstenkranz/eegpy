#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for accessing F32 files remotely. Consist of a server (F32Server) and a proxy/client (F32Proxy). 
The server reads F32-files saved in a certain directory and exposes some rudimentary methods for reading.
The proxy communicates to this server and features a reduced interface similar to eegpy.F32-objects.
"""

import eegpy
from eegpy.misc import FATALERROR, debug
import eegpy.analysis.phases
from eegpy.formats.open_eeg import open_eeg

try:
    import numpy as n
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

from SimpleXMLRPCServer import SimpleXMLRPCServer
import xmlrpclib
import base64
import pickle
import socket
socket.setdefaulttimeout(200)
from sys import stdout


#import BaseHTTPServer
#import httpserver
#from pyftpdlib import ftpserver
from threading import Thread, Lock 
import os, os.path, time
from sys import stdout

port = 51235



class F32Server(SimpleXMLRPCServer):
    """The F32 server"""
    def __init__(self,dirname,*args,**kwargs):
        """Initialisation. Needs two arguments:
        Parameters:
         * dirname
            name of the directory where the F32 files are stored
         * ip-configuration
            tuple of (host, port) where to listen. Use ("",port) for listening on all interfaces.
        """
        SimpleXMLRPCServer.__init__(self,*args,**kwargs)
        self.timeout=60
        if not os.path.exists(dirname):
            raise ValueError("Illegal value for dirname: path doesn't exist!")
        #self._locks
        self._dirname = dirname
        self.register_function(self.getData)
        self.register_function(self.getChannelData)
        self.register_function(self.getChannelNames)
        self.register_function(self.getNumDatapoints)
        self.register_function(self.getShape)
        self.register_function(self.getSamplRate)

    def getData(self,fn,start,length=None,stride=1,channelList=None):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        if length==None:
            length=eeg.num_datapoints-start
        ar = n.array(eeg.getData(start,length,stride,channelList))
        eeg.close()
        #return (base64.encodestring(ar.tostring()),ar.shape,str(ar.dtype))
        return base64.encodestring(pickle.dumps(ar,-1))

    def getChannelData(self,fn,channel):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        start=0
        length=eeg.num_datapoints-start
        stride=1
        ar = n.array(eeg.getData(start,length,stride,[channel]))
        eeg.close()
        return base64.encodestring(pickle.dumps(ar))

    def getChannelNames(self,fn):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        rv = eeg.channel_names
        eeg.close()
        return rv

    def getNumDatapoints(self,fn):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        rv = eeg.num_datapoints
        eeg.close()
        return rv

    def getShape(self,fn):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        rv = eeg.shape
        eeg.close()
        return rv

    def getSamplRate(self,fn):
        eeg = eegpy.F32(os.path.join(self._dirname,fn))
        rv = eeg.get_samplRate
        eeg.close()
        return rv
  
class F32Proxy(object):
    """Provides a reduced interface for accessing F32-files that are served via a F32Server instance.
    Supported methods:
     * getData
     * getChannelData
     * getChannelNames
     * getNumDatapoints
    Properties:
     * channel_names
     * num_channels
     * num_datapoints

    Example:
    from eegpy.computing.f32oip import F32Proxy
    fp = F32Proxy("hafner_bp.f32","http://10.47.8.201:51235")
    ndp = fp.num_datapoints
    ar = fp.getData(0,10000,1,[2,7,10]) #Get 10000 samples starting at sample 0 with stride 1 from channels 2,7,10
    ar = fp.getData(50000,1000) #Get 1000 samples starting at sample 50000 with stride 1 from all channels
    ar2 = fp.getChannelData(fp.channel_names.find("TL01")) #Get all data for channel TL01
    """
    def __init__(self,fn,*args,**kwargs):
        self._proxy = xmlrpclib.ServerProxy(*args,**kwargs)
        self._fn = fn

    def getData(self,start,length,stride=1,channelList=None):
        rv = self._proxy.getData(self._fn,start,length,stride,channelList)
        return pickle.loads(base64.decodestring(rv))

    def getChannelData(self,channel):
        rv = self._proxy.getChannelData(self._fn,channel)
        return pickle.loads(base64.decodestring(rv))

    def getChannelNames(self):
        return self._proxy.getChannelNames(self._fn)

    def getNumDatapoints(self):
        return self._proxy.getNumDatapoints(self._fn)

    def getNumChannels(self):
        return self._proxy.getShape(self._fn)[1]

    def getShape(self):
        return self._proxy.getShape(self._fn)

    num_channels = property(getNumChannels)
    num_datapoints = property(getNumDatapoints)
    shape = property(getShape)

if __name__ == "__main__":
    server = F32Server("/home/fell/thorsten/turmserver_neu/data/",("",port))
    server.serve_forever()
    
