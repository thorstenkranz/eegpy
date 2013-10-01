#!/usr/bin/env python
# -*- coding: utf-8 -*-

import threading
import Queue

import eegpy 
from eegpy.misc import FATALERROR, debug
#from eegpy.events import EventTable

try:
    import numpy as np
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    
class EEG_file(object):
    _fn = None
    _channel_names = None
    _channel_units = None
    _shape = (0,0)
    _mode="r"
    
    def __init__(self, fn):
        pass
    
    def close(self):
        if debug:
            print "Method EEG_file.close() invoked"
        
    def __getitem__(self, item):
        """Wrapper for easy access to the getData methods. A Slice should be given."""
        if debug:
            print "Method EEG_file.__getitem() invoked"
    
    #TODO: Method for retrieving 
    def get_data_for_events(self,event_list,start_end = (-1000,2000),channels = None):
        """For a given list of events, this method collects the data from the eeg file 
        and returns them as numpy-array."""
        assert len(event_list)>0, "Event list is not a list"
        assert len(event_list)<10000, "Only up to 10000 timepoints support right now"
        assert start_end[0]<start_end[1], "Incorrect values for start and end"
        
        start = start_end[0]
        end = start_end[1]
        event_list = [int(x) for x in event_list]
        if (channels == None):
            channels = range(self.num_channels)
        channels.sort()
        rv = np.zeros((end-start,len(channels),len(event_list)),"d")
        arb = np.zeros((self.num_channels),np.bool)
        for c in channels:
            arb[c]=True
                     
        for i,t in enumerate(event_list):
            if t+start<0 or t+end>self.num_datapoints:
                raise IndexError("Cannot get data from %i to %i" % (t+start,t+end) )
            rv[:,:,i] = self[t+start:t+end,arb]
        return rv
    
    def moving_windows(self,channels = None, size = None, overlap = None, threads=False):
        """Generator for iterating through the file, moving-window-technique.
        threads: bool
          Use separate thread for reading."""

        class ThreadedWindowReader(threading.Thread):

            def __init__(self,eeg,ch_bs,size,overlap,qsize=10):
                threading.Thread.__init__(self)
                self.eeg = eeg
                self.ch_bs = ch_bs
                self.size=size
                self.overlap=overlap
                self.q = Queue.Queue(qsize)
                self.has_more = True
                self.is_running = False
            
            def run(self):
                start=0
                self.is_running = True
                while start+self.size<self.eeg.num_datapoints:
                    data = self.eeg[start:start+self.size,self.ch_bs]
                    self.q.put((data, start))
                    start = start+self.size-self.overlap
                self.is_running=False

            def get(self):
                if self.is_running:
                    data,start = self.q.get()
                else:
                    data,start = self.q.get_nowait()
                return data, start

        if channels==None:
            channels=range(self.num_channels)
        if size==None:
            size=4096
        if overlap==None:
            overlap=0
        start = 0
        ch_bs = np.zeros((self.num_channels),np.bool)
        for ch in channels:
            ch_bs[ch] = True
        if not threads:
            while start+size<self.num_datapoints:
                data = self[start:start+size,ch_bs]
                yield data, start
                start = start+size-overlap
        else:
            twr = ThreadedWindowReader(self,ch_bs,size,overlap)
            twr.start()
            try:
                while True:
                    yield twr.get()
            except Queue.Empty:
                print "Queue now empty"



    
    def moving_windows_num(self, size = None, overlap = None):
        """Helper functions. Calulates the number of windows needed to cover the whole file.
        This might be useful for allocating result-arrays prior to iteration"""
        if size==None:
            size=4096
        if overlap==None:
            overlap=0
        
        return (self.num_datapoints)/(size-overlap) #Not sure if it's right, need testing
    
    def get_shape(self):
        return self._shape

    def get_samplRate(self):
        return 1.0
    
    def get_num_datapoints(self):
        return 0
    
    def get_num_channels(self):
        return self.shape[1]
    
    def get_channel_names(self):
        return []
    
    def set_channel_names(self):
        pass
    
    def getChannelNames(self):
        return []
    
    def edit_channel_names(self):
        try:
            from eegpy.ui.widgets.edit_channel_names import ChannelNamesEditor
        except ImportError:
            print "Cannot load ChannelNamesEditor. Probably TraitsUi is not installed."
            raise NotImplementedError()
        ce = ChannelNamesEditor(channel_names = self.channel_names)
        ce.configure_traits()
        if len(self.channel_names) == len(ce.channel_names):
            self.channel_names = ce.channel_names
        else:
            raise ValueError("Number of channel-names may not be changed!")
     
    def get_filename(self):
        return self._fn    
    
    shape = property(get_shape)
    samplRate = property(get_samplRate)
    Fs = property(get_samplRate)
    num_datapoints = property(get_num_datapoints)
    numDatapoints = property(get_num_datapoints)
    num_channels = property(get_num_channels)
    channel_names = property(get_channel_names, set_channel_names)
    channelNames = property(get_channel_names)
    fn = property(get_filename)
    #TODO: Add functions for direct access to attributes like numChannels etc. 

class MemoryMappedBinaryDataFile(EEG_file): 
    _mm = None #memmap-object
    
    def __init__(self, filename, mode, shape, data_offset=0, dtype="f", Fs = 1000.0, 
                 channel_names = None):
        assert mode in ["r","r+","w+"], "Unsupported mode"
        self._fn = filename
        
        try:
            self.f = open(self._fn, mode)
        except:
            raise IOError, "The file could not be opened!"
        
        assert shape != None, "Shape must be given."
        self._mode = mode
#        if mode in ["r","r+"]:
        if any(np.array(shape) < 1):
            raise ValueError("Only values > 0 are allowed for number of channels and number of datapoints.")
        self._shape = shape
        if channel_names is not None:
            if len(channel_names) != shape[1]:
                raise ValueError("Lengths of channel names and memmap shape 1 should match.")
            self.set_channel_names(channel_names)
        else:
            self.set_channel_names([str(x+1) for x in range(shape[1])])  
        self._data_offset = data_offset
        self._Fs = Fs
        self._dtype = dtype
        self._mode = mode
        self._open_memmap()
        
    def _open_memmap(self):
        self._mm = np.memmap(self._fn,dtype=self._dtype,offset=self._data_offset,
                             shape=self._shape,mode=self._mode)
        
    def _reopen_memmap(self):
        mm_mode = "r" if self._mode == "r" else "r+"
        self._mm = np.memmap(self._fn,dtype=self._dtype,offset=self._data_offset,
                             shape=self._shape,mode=mm_mode)
          
    def __del__(self):
        self.close()
    
    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()
    
    def __str__(self):
        return "Memory-mapped data file at %s" % self._fn
 
    def close(self):
        self.f.close()
        if self._mm != None:
            self._mm.flush()
            del self._mm
            self._mm = None
    
    @property
    def is_open(self):
        if self._mm is None:
            return False
        if type(self._mm)==np.memmap:
            return True
        return False
                    
    def __getitem__(self,item):
        #for retrieving data from the memmap-object
        return self._mm.__getitem__(item)
    
    def __setitem__(self,key,item):
        return self._mm.__setitem__(key,item)
    
    def getData(self, start, length, stride=1, channelList = None):
        """For conformity with older classes.
        Now uses boolean indexing"""
        # No channelList given
        if (channelList == None):
            arb = np.ones((self.shape[1]),np.bool)
        else:
            arb = np.zeros((self.shape[1]),np.bool)
            channelList.sort()
            if (channelList[0]<0 or channelList[-1]>(self.numChannels-1)):
                raise ValueError("Channellist contains illegal channel-number.\nmin=%i, max=%%i" %(channelList[0],channelList[-1])) 
            for c in channelList:
                arb[c] = True

        tmpData = self[start:start+(length*stride):stride,:]
        rv = tmpData[:,arb]
        return rv
    
    def getOverviewData(self, numPoints, channelList = None):
        if(numPoints > self.numDatapoints):
            raise ValueError("To many datapoints for overview data were requested.")
        # Keine channelList übergeben
        if (channelList == None):
            channelList = range(self.numChannels)
        channelList.sort()
        if (channelList[0]<0 or channelList[-1]>(self.numChannels-1)):
            raise Exception, "Kanalliste enthaelt falsche Kanalnummern.\min=%i, max=%%i" %(channelList[0],channelList[-1]) 
        #Waehle geeignete Parameter und rufe damit getData auf
        #start=0, length=numPoints, stride=?
        stride = self.numDatapoints / numPoints
        return self.getData(0,numPoints,stride,channelList)
        
    def get_samplRate(self):
        return self._Fs
    
    def get_num_datapoints(self):
        return self._shape[0]
    
    def get_channel_names(self):
        return self._channel_names
    
    def getChannelNames(self):
        return self._channel_names
    
    def set_channel_names(self, cns):
        assert len(cns) == self.shape[1], "List of channelnames must contain exactly %i elements" % self.shape[1]
        assert self._mode in ["r+", "w+"], "Cannot set channel_names: F32 is read-only"
        self._channel_names = [str(x) for x in cns]
            
    @property
    def Fs(self):
        """Sampling rate of recording. 

        .. hint::
            The F32-format doesn' allow per-channel sampling rates."""
        return self.get_samplRate()

    samplRate = property(get_samplRate)
    num_datapoints = property(get_num_datapoints)
    numDatapoints = property(get_num_datapoints)
    channel_names = property(get_channel_names, set_channel_names)
    channelNames = property(get_channel_names)     

class EEGfiltered(object):
    """Object with included frequency-filters"""
    def __init__(self,filter_function):
        self._filter_function = filter_function

    @property
    def ff(self):
        return self._filter_function

    def __getitem__(self,item):
        """Calls method of other class, filters the return value.
        Make sure you also inherit from MemoryMappedBinaryDataFile!"""
        return self.ff(MemoryMappedBinaryDataFile.__getitem__(self,item))

if __name__=="__main__":
    print "Test speed for moving_windows with and without multithreading"
    import time
    eeg = eegpy.F32("/media/story/SchlafGed/iEEG/data/canseven_bp.f32")
    num_windows = eeg.moving_windows_num(size=4096,overlap=1000)
    start1 = time.time()
    for i,(w,start) in enumerate(eeg.moving_windows(size=4096,overlap=1000)):
        print "%i / %i" %(i,num_windows)
        time.sleep(0.01)
    t_normal = time.time()-start1

    start2 = time.time()
    for i,(w,start) in enumerate(eeg.moving_windows(size=4096,overlap=1000,threads=True)):
        print "%i / %i" %(i,num_windows)
        time.sleep(0.01)
    t_threads = time.time()-start2

    print t_normal, t_threads

