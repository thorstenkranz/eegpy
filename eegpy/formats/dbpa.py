#! /usr/bin/python
# -*- coding: utf8 -*-
import struct
import os
import StringIO
import tempfile

import numpy as n
#from scipy.signal import butter

from eegpy.formats.iobase import EEG_file

class DBPA(EEG_file):
    """Access to Sensorium DBPA datafiles
    """

    header = {"formatName":"Bonn F32 Data Format",
              "formatVer":"1.0.00",
              "sbChinfo":1024, # Weil header 1024 bytes hat
              "sbData":None,
              "numChannels":0,
              "numDatapoints":0,
              "samplRate":1.0,
              "tStart":0.0,
              "tEnd":0.0,
              "reservedString":"",
              "datatype":4096,
              "samplesize":17533,
              "minsample":0.0,
              "maxsample":0.0,
              "reserved2":""}
    reserved3 = ""
    fmtSample = "", # Format-String für einzelnes Sample
    #(self.formatName,self.formatVer,self.sbChinfo,self.sbData,self.numChannels,self.numDatapoints,self.samplRate,self.tStart,self.tEnd,self.reservedString, self.datatype, self.samplesize, self.minsample, self.maxsample,self.reserved2) = struct.unpack(fmtF32header,s)
    #self.fmtSample = "= %if" % self.numChannels
    #shape = property(self.get_shape)
    dataToWrite = None
    headerWritten = False
    _mm = None #memmap-object

    def __init__(self,filename,mode="r+",cNames=None,shape=None, Fs = 1000.0):
        assert mode in ["r","r+","w+"], "Unsupported mode"
        self._fn = filename

        try:
            self.f = open(self._fn, mode)
        except:
            raise IOError, "Die angegebene Datei konnte nicht geöffnet werden!"

        self._mode = mode
        if mode in ["r","r+"]:
            s = self.f.read(struct.calcsize(fmtF32header))
            (self.header["formatName"],self.header["formatVer"],self.header["sbChinfo"],self.header["sbData"],self.header["numChannels"],self.header["numDatapoints"],self.header["samplRate"],self.header["tStart"], self.header["tEnd"], self.header["reservedString"], self.header["datatype"], self.header["samplesize"], self.header["minsample"], self.header["maxsample"],self.header["reserved2"]) = struct.unpack(fmtF32header,s)

            self.numChannels=self.header["numChannels"]
            #print self.numChannels
            if self.numChannels > 500 or self.numChannels < 1:
                raise Exception, "The chosen file is either in an unrecognized format or corrupt!"
            if self.header["sbData"] != (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo)):
                print "Der Wert für das Anfangsbyte der Daten war falsch gesetzt. Korrigiere..."
                self.sbData = (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo))
            self.fmtSample = "= %if" % self.numChannels
            self._channel_names = []
            self._channel_units = []
            for i in range(0,self.numChannels):
                s = self.f.read(struct.calcsize(fmtF32channelinfo))
                (chName, chUnit, reserved3) = struct.unpack(fmtF32channelinfo,s)
                self._channel_names.append(chName)
                self._channel_units.append(chUnit)
            self._shape = (self.header["numDatapoints"],self.header["numChannels"])
            #FIXED: closinf self.f before opening memmap.
            self.f.close()
            self._mm = n.memmap(self._fn,dtype=n.float32,offset=self.header["sbData"],shape=self._shape,mode=mode)
        elif mode=="w+":
            assert shape != None, "Shape must be given."
            assert len(shape) == 2, "Only 2d is possible"
            self._shape = shape
            self.header["numChannels"] = self._shape[1]
            self.numChannels = self._shape[1]
            if cNames != None:
                assert len(cNames)>0, "Liste der Kanalnamen muss länger als 0 sein!"
                assert shape[1] == len(cNames)
                self._channel_names = [str(x) for x in cNames]
            self.header["sbData"] = int(self.header["sbChinfo"]+self.numChannels*struct.calcsize(fmtF32channelinfo))
            self.header["numChannels"] = self._shape[1]
            self.header["numDatapoints"] = self._shape[0]
            self.header["samplRate"] = float(Fs)
            self.f.close()
            self._mm = n.memmap(self._fn,dtype=n.float32,offset=self.header["sbData"],shape=self._shape,mode=mode)
            del self._mm
            self._mm = None
            self.f = open(self._fn, "r+")
            self.writeHeader()
            self.f.close()
            self._mm = n.memmap(self._fn,dtype=n.float32,offset=self.header["sbData"],shape=self._shape,mode="r+")

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()

    def __str__(self):
        rv = StringIO.StringIO()
        ks_header = self.header.keys()
        ks_header.sort()
        print >>rv, "eegpy F32-object\n----------------\n"
        print >>rv, "Header-Information:"
        for k in ks_header:
            print >>rv, "%s: %s" %(str(k),str(self.header[k]))
        print >>rv, "Channels:"
        for i,c in enumerate(self.channel_names):
            print >>rv, "%3i,%s" % (i,c)
        return rv.getvalue()

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
        if type(self._mm)==n.memmap:
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
        # Keine channelList übergeben
        arb = n.zeros((self.shape[1]),n.bool)
        if (channelList == None):
            arb = n.ones((self.shape[1]),n.bool)
        else:
            channelList.sort()
            if (channelList[0]<0 or channelList[-1]>(self.numChannels-1)):
                raise Exception, "Channellist contains illegal channel-number.\min=%i, max=%%i" %(channelList[0],channelList[-1])
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


    def writeHeader(self):
        """Writing the header"""
        #Preparation: chNames setzen, falls leer
        if self._channel_names == None:
            self._channel_names = ["%i"%(i+1) for i in range(self._shape[1])]
            self.numChannels = len(self._channel_names)
        #Preparation: chUnits setzen, falls leer
        if self._channel_units == None:
            self._channel_units = ["" for i in self._channel_names]

        self.f.seek(0)
        s=struct.pack(fmtF32header, self.header["formatName"],self.header["formatVer"],int(self.header["sbChinfo"]),self.header["sbData"],self.header["numChannels"],self.header["numDatapoints"],self.header["samplRate"],self.header["tStart"],self.header["tEnd"],self.header["reservedString"], self.header["datatype"], self.header["samplesize"], self.header["minsample"], self.header["maxsample"],self.header["reserved2"])
        self.f.write(s)
        for i in range(len(self._channel_names)):
            #print len(self._channel_names), i
            s=struct.pack(fmtF32channelinfo, self._channel_names[i], self._channel_units[i], self.reserved3)
            self.f.write(s)

        if(not (self.f.tell()==self.header["sbData"])):
            raise Exception, "Irgendwie ist die Position falsch, %i != %i" % (self.f.tell(),self.header["sbData"])
        self.headerWritten = True

    def get_samplRate(self):
        return self.header["samplRate"]

    def get_num_datapoints(self):
        return self.header["numDatapoints"]

    def get_channel_names(self):
        return self._channel_names

    def getChannelNames(self):
        return self._channel_names

    def set_channel_names(self, cns):
        assert len(cns) == self.shape[1], "List of channelnames must contain exactly %i elements" % self.shape[1]
        assert self._mode in ["r+", "w+"], "Cannot set channel_names: F32 is read-only"
        try:
            del self._mm
            self._mm = None
        except Exception, e:
            pass
        self._channel_names = [str(x) for x in cns]
        self.f = open(self._fn, "r+")
        self.writeHeader()
        self.f.close()
        self._mm = n.memmap(self._fn,dtype=n.float32,offset=self.header["sbData"],shape=self._shape,mode="r+")

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

class F32filtered(F32):
    """F32 read-only object with included frequency-filteres"""
    def __init__(self,filename,filter_function):
        F32.__init__(self,filename,"r+")
        self._filter_function = filter_function

    @property
    def ff(self):
        return self._filter_function

    def __getitem__(self,item):
        """Calls method of super-class, filters the return value"""
        return self.ff(F32.__getitem__(self,item))

class F32ReaderFiltered(F32Reader):
    """F32 read-only object with included frequency-filteres.
    Uses old F32Reader, without memmap."""
    def __init__(self,filename,filter_function):
        F32Reader.__init__(self,filename)
        self._filter_function = filter_function

    @property
    def ff(self):
        return self._filter_function

    def getData(self,start,length,stride=1,channelList = None):
        """Calls method of super-class, filters the return value"""
        return self.ff(F32Reader.getData(self,start,length,stride,channelList))

class F32filteredWithCache(F32):
    """F32 read-only object with included frequency-filteres"""
    def __init__(self,filename,btype='lp',fl=None,fh=None,border=2,filter_windowed=False,dirname=None):
        from eegpy.ui.eegpylab import freqfilt_eeg
        fd, tmpfn = tempfile.mkstemp(dir=dirname)
        freqfilt_eeg(filename,tmpfn,btype,fl,fh,border,filter_windowed)
        F32.__init__(self,tmpfn)

    def __del__(self):
        self.close()
        try:
            os.unlink(self._fn)
        except Exception,e:
            print "Cannot remove file %s"%self._fn, e



