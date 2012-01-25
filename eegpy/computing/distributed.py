#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module is supposed to provide sophisticated possibilites to use distributed computing 
for the analysis of large eeg-files like in epilepsy-research. 
It's meant to work in a decentralized way, every user can start jobs and use the power of the cluster"""

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
import pickle
import socket
from sys import stdout

import Pyro.core
import Pyro.errors

#import BaseHTTPServer
#import httpserver
#from pyftpdlib import ftpserver
from threading import Thread, Lock 
import os, os.path, time
from sys import stdout

port = 51235

class Bivariate(Pyro.core.ObjBase):
    """Pyro-Object. Is exposed via network-connection and offers functionality 
    for computation of bivariate measures."""
    _available = True
    
    def __init__(self):
        Pyro.core.ObjBase.__init__(self)
    
    def tell_available(self):
        return self._available
    
    def phase_coherence(self,x,y):
        self._available = False
        rv = eegpy.analysis.phases.phase_coherence(x, y)
        print rv
        self._available = True
        #time.sleep(2)
        return rv
    
    def average_distance(self,x,y):
        """Only for testing purposes. Trivial functionality."""
        self._available = False
        rv = (x-y).mean()
        print rv
        self._available = True
        return rv
        

class TestProxy:
    def __init__(self):
        test = Pyro.core.getProxyForURI("PYROLOC://localhost:%s/test" % port)
        test.tell_available()
        
class EEGPyDCClient:
    """Scans the network on the given port for a valid connection"""
    ips = ["10.47.13.%s"%x for x in range(1,50)] + ["10.47.8.%s"%x for x in range(200,220)] + ["10.47.12.%s"%x for x in range(1,200)]
    _responding = []
    _lock_responding = Lock()
    #_available = []
    #_lock_available = Lock()
    _searcher = None
    _workers = []
    _num_threads = 0
    _max_threads = 10
    _stop = False
    
    
    def __init__(self):
       self._searcher = Thread(target=self.search_servers)
       self._searcher.start()
       time.sleep(2)
       #self._checker = Thread(target=self.check_servers)
       #self._checker.start()
    
    def __del__(self):
        self._stop = True
       
    def shutdown(self):
        self._stop = True
            
    def search_servers(self):
        while not self._stop:
            for ip in self.ips:
                self._lock_responding.acquire()
                try:
                    sock = socket.socket(socket.AF_INET)
                    sock.settimeout(0.01)
                    sock.connect((ip,port))
                    sock.close()
                    proxy = Pyro.core.getProxyForURI("PYROLOC://%s:%i/bivariate" % (ip,port))
                    #proxy._setTimeout(3)
                    try:
                        if proxy.tell_available():
                            #print "%s: Server gefunden." % ip
                            pass
                        if not ip in self._responding:
                            self._responding.append(ip)
                    except Pyro.errors.ProtocolError, e:
                        #print "%s: Keine Verbindung" % ip, e
                        if ip in self._responding:
                            self._responding.remove(ip)
                except Exception, e:
                    #print "%s: Fehler beim Socket-Verbinden" % ip , e
                    if ip in self._responding:
                            self._responding.remove(ip)
                self._lock_responding.release()
            #print "r: ", self._responding
            #print "a: ", self._available
            time.sleep(10)
            
    def find_server(self,varname):
        for ip in self._responding:
            try:
                sock = socket.socket(socket.AF_INET)
                sock.settimeout(0.01)
                sock.connect((ip,port))
                sock.close()
                proxy = Pyro.core.getProxyForURI("PYROLOC://%s:%i/%s" % (ip,port,varname))
                try:
                    if proxy.tell_available():
                        #print "%s: Server gefunden." % ip
                        return ip
                except Pyro.errors.ProtocolError, e:
                    #print "%s: Keine Verbindung" % ip, e
                    pass
            except Exception, e:
                pass
        raise RuntimeError("Can't find any server")
    
    def collect_results(self, results):
        for i,w in enumerate(self._workers):
            if w.finished:
                results[w._result_indices] = w.result
                #print w.result
                del self._workers[i]
                
    def bivariate_measure(self,fn, measure_name, symmetric=True, channels = None, window_size=None, overlap=None):
        """Iterates through the eeg file and calculates a bivariate measure 
        for all channel-combinations for the given channels
        Works only if a method with the name "measure_name" exists for the remote object.
        If symmetric==True, only the necessary calculations are made and all redundant 
        calculations are omitted, otherwise all pairwise combinations are calculated."""
        try:
            f = open_eeg(fn)
        except Exception,e:
            print "Cannot open eeg-file"
        if channels==None:
            num_channels = f.numChannels
        else:
            num_channels = len(channels)
            
        results = n.ones((f.moving_windows_num(window_size,overlap),num_channels,num_channels),"d")
        
        i=0
        for data,pos in f.moving_windows(channels,window_size,overlap):
            print pos
            print i
            stdout.flush()
            for j in range(num_channels):
                if symmetric:
                    upper_k = j
                else:
                    upper_k = num_channels
                print upper_k
                stdout.flush()
                for k in range(upper_k):
                    while True:
                        if len(self._workers)<=self._max_threads and len(self._workers)<len(self._responding):
                            try:
                                #print "Point 1"
                                ip = self.find_server("bivariate")
                                #print "Point 2"
                                self._workers.append(WorkThread(ip,measure_name,data[:,j],data[:,k]))
                                self._workers[-1]._result_indices = (i,j,k)
                                self._workers[-1].start()
                                #print "Point 3"
                                #print "anzahl Worker:", len(self._workers)
                                break
                            except RuntimeError, e:
                                #Couldn't find a server
                                pass
                        self.collect_results(results)
                        #print "Point 4"
                            
                    #results[i,k,j] = results[i,j,k]

            i+=1
        while len(self._workers)>0:
            self.collect_results(results)
        if symmetric:
            for i in range(results.shape[0]):
                for j in range(results.shape[1]):
                    for k in range(j):
                        results [i,k,j] = results[i,j,k] 
        return results
    
    def univariate_measure(self,fn, measure_name, channels = None, window_size=None, overlap=None):
        """Iterates through the eeg file and calculates a bivariate measure 
        for all channel-combinations for the given channels
        Works only if a method with the name "measure_name" exists for the remote object.
        If symmetric==True, only the necessary calculations are made and all redundant 
        calculations are omitted, otherwise all pairwise combinations are calculated."""
        try:
            f = open_eeg(fn)
        except Exception,e:
            print "Cannot open eeg-file"
        if channels==None:
            num_channels = f.numChannels
        else:
            num_channels = len(channels)
            
        results = n.ones((f.moving_windows_num(window_size,overlap),num_channels),"d")
        
        i=0
        for data,pos in f.moving_windows(channels,window_size,overlap):
            print pos
            print i
            stdout.flush()
            for j in range(num_channels):
                while True:
                    if len(self._workers)<=self._max_threads and len(self._workers)<len(self._responding):
                        try:
                            #print "Point 1"
                            ip = self.find_server("univariate")
                            #print "Point 2"
                            self._workers.append(WorkThread(ip,measure_name,data[:,j]))
                            self._workers[-1]._result_indices = (i,j,k)
                            self._workers[-1].start()
                            #print "Point 3"
                            #print "anzahl Worker:", len(self._workers)
                            break
                        except RuntimeError, e:
                            #Couldn't find a server
                            pass
                    self.collect_results(results)
                    #print "Point 4"
                        
                #results[i,k,j] = results[i,j,k]

            i+=1
        while len(self._workers)>0:
            self.collect_results(results)
        if symmetric:
            for i in range(results.shape[0]):
                for j in range(results.shape[1]):
                    for k in range(j):
                        results [i,k,j] = results[i,j,k] 
        return results

    def phase_coherence(self,fn, measure_name,channels = None, window_size=None, overlap=None):
        return self.bivariate_measure(fn, "phase_coherence", True, channels, window_size, overlap)
    
class WorkThread(Thread):        
    """Increases the functionality of Thread. Calls the function given as argument "work" 
    and saves the result as property "result" """
    
    result = None
    finished = False
    _result_indices = None #Where to put the results
    _work = None
    _ip = None
    
    def __init__(self,ip,work,*args,**kwargs):
        Thread.__init__(self)
        self._ip = ip
        self._work = work
        self.myargs = args
        self.mykwargs = kwargs
        
    def run(self):
        try:
            proxy = Pyro.core.getProxyForURI("PYROLOC://%s:%i/%s" % (self._ip,port,"bivariate"))
            self.result = getattr(proxy, self._work)(*(self.myargs), **(self.mykwargs))
        except Exception, e:
            print e
            pass
        finally:
            self.finished = True
        

    
if __name__ == "__main__":
    #server = EEGPyDCServer()
    #server.serve_forever()
    Pyro.core.initServer()
    daemon=Pyro.core.Daemon(port=port)
    uri=daemon.connect(Bivariate(),"bivariate")
    
    print "The daemon runs on port:",daemon.port
    print "The object's uri is:",uri
    
    daemon.requestLoop()
