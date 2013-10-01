#! /usr/bin/python
# -*- coding: utf8 -*-
import struct
import os
import StringIO
import tempfile

import numpy as np
#from scipy.signal import butter

from eegpy.formats.iobase import (MemoryMappedBinaryDataFile,
                    EEGfiltered)

fmtF32header = "= 21p 7p i i H i d d d 13p i H f f 931p"
#Erklaerung
# = Native byte order, standard size and alignment
# 21p: Format-Name
# 7p:  Format-Version
# i:   start-byte der Kanalinformationen
# i:   start-byte der Daten
# H:   Anzahl Kanaele unsigned short
# i:   Anzahl Datenpunkte
# d d d: Sampling-rate, Startzeit, Endzeit (Zeiten als Double, integraler Teil ist Anzahl Tage seit 12/30/1899, hinter komma dann bruchteil des Tages, der abgelaufen ist
# 13p: reserviert
# i:   datatype
# H:   samplesize
# f:   minsample
# f:   maxsample
# 931p:reserviert2, um header auf 1024b aufzublasen 

fmtPosOfNumDatapoints = "= 21p 7p i i H "
#Zeigt die Position von NumDatapoints im header an, damit diese nachträglich geschrieben werden kann.

fmtF32channelinfo = "= 256p 32p 224p"
# 256p: Kanalname
# 32p:  Einheit der Messgroesse (z.B. Volt)
# 224p: reserviert um Kanalinfo auf 512 bytes zu vergroessern

from _f32_old import F32Reader, F32Writer, F32WriterAdvanced
 
class F32(MemoryMappedBinaryDataFile):
    """Access to F32-datafiles
    
This is the new class used for both read and write access. 
It uses numpy.memmap for rapid read/write"""

    header = {"formatName":"Bonn F32 Data Format",
              "formatVer":"1.0.00",
              "sbChinfo":1024, # Because header has 1024 bytes
              "sbData":None, # Where the binary data start
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
        
        if mode in ["r","r+"]:
            s = self.f.read(struct.calcsize(fmtF32header))
            (self.header["formatName"],self.header["formatVer"],self.header["sbChinfo"],self.header["sbData"],self.header["numChannels"],self.header["numDatapoints"],self.header["samplRate"],self.header["tStart"], self.header["tEnd"], self.header["reservedString"], self.header["datatype"], self.header["samplesize"], self.header["minsample"], self.header["maxsample"],self.header["reserved2"]) = struct.unpack(fmtF32header,s)
            
            self.numChannels=self.header["numChannels"]
            #print self.numChannels
            if self.numChannels > 500 or self.numChannels < 1:
                raise Exception, "The chosen file is either in an unrecognized format or corrupt!"
            if self.header["sbData"] != (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo)):
                #print "Der Wert für das Anfangsbyte der Daten war falsch gesetzt. Korrigiere..."
                self.header["sbData"] = (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo))
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
            
            shape_from_header = (self.header["numDatapoints"], self.header["numChannels"])
            MemoryMappedBinaryDataFile.__init__(self, filename, mode, shape_from_header, 
                            data_offset = self.header["sbData"], 
                            dtype=np.float32, Fs = self.header["samplRate"],
                            channel_names=cNames)
        elif mode=="w+":
            assert shape != None, "Shape must be given."
            assert len(shape) == 2, "Only 2d is possible"
            if cNames != None:
                assert len(cNames)>0, "List of channel names must be longer than 0!"
                assert shape[1] == len(cNames) 
                self._channel_names = [str(x) for x in cNames]  
            self.header["sbData"] = int(self.header["sbChinfo"]+shape[1]*struct.calcsize(fmtF32channelinfo))   
            self.header["numChannels"] = shape[1]
            self.header["numDatapoints"] = shape[0]
            self.header["samplRate"] = float(Fs)
            self.f.close()            
            MemoryMappedBinaryDataFile.__init__(self, filename, mode, shape, 
                            data_offset = self.header["sbData"], 
                            dtype=np.float32, Fs = Fs,
                            channel_names=cNames)
            del self._mm
            self._mm = None
            self.f = open(self._fn, "r+")
            self.writeHeader()
            self.f.close()
            self._reopen_memmap()

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
 
    def set_channel_names(self, cns):
        assert len(cns) == self.shape[1], "List of channelnames must contain exactly %i elements" % self.shape[1]
        assert self._mode in ["r+", "w+"], "Cannot set channel_names: F32 is read-only"
        try:
            del self._mm
            self._mm = None
        except Exception:
            pass
        self._channel_names = [str(x) for x in cns]
        self.f = open(self._fn, "r+")
        self.writeHeader()
        self.f.close()
        self._mm = np.memmap(self._fn,dtype=np.float32,offset=self.header["sbData"],shape=self._shape,mode="r+")
        
class F32filtered(EEGfiltered, F32):
    """F32 read-only object with included frequency-filteres"""
    def __init__(self,filename,filter_function, **kw_args):
        F32.__init__(self,filename,"r+", **kw_args)
        EEGfiltered.__init__(self, filter_function)

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

__all__ = ["F32", "F32filtered", "F32Reader", "F32ReaderFiltered", "F32filteredWithCache", "F32Writer", "F32WriterAdvanced"]
    
