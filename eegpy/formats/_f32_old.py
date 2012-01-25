#! /usr/bin/python
# -*- coding: utf8 -*-
import struct
#import Gnuplot
import time
import numpy as n
#import matplotlib.cm as cm
#import pywt
#import pylab

#fn="/home/thorsten/f32/1.f32"    

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

class F32Reader:
    "Class for reading EEG-Files which are saved in the Bonn-F32-Format"
    def __init__(self, filename):
        # Initialisierung der Instanzvariablen
        self.fn = filename
        try:
            self.f = open(self.fn, "rb")
        except:
            raise IOError, "Die angegebene Datei wurde nicht gefunden!"
        self.formatName = ""
        self.formatVer = ""
        self.sbChinfo = 0
        self.sbData = 0
        self.numChannels = 0
        self.numDatapoints = 0
        self.samplRate = 0.0
        self.tStart = 0.0
        self.tEnd = 0.0
        self.reservedString = ""
        self.datatype = 0
        self.samplesize = 0
        self.minsample = 0.0
        self.maxsample = 0.0
        self.reserved2 = ""
        
        self.fmtSample = "" # Format-String für einzelnes Sample
        
        s = self.f.read(struct.calcsize(fmtF32header))
        (self.formatName,self.formatVer,self.sbChinfo,self.sbData,self.numChannels,self.numDatapoints,self.samplRate,self.tStart,self.tEnd,self.reservedString, self.datatype, self.samplesize, self.minsample, self.maxsample,self.reserved2) = struct.unpack(fmtF32header,s)
        if self.numChannels > 500 or self.numChannels < 1:
            raise Exception, "Die gewählte Datei ist entweder in einem falschen Format oder korupt!"
        if self.sbData != (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo)):
            print "Der Wert für das Anfangsbyte der Daten war falsch gesetzt. Korrigiere..."
            self.sbData = (struct.calcsize(fmtF32header)+self.numChannels*struct.calcsize(fmtF32channelinfo))
        self.fmtSample = "= %if" % self.numChannels
        self.channelNames = []
        self.channelUnits = []
        for i in range(0,self.numChannels):
            s = self.f.read(struct.calcsize(fmtF32channelinfo))
            (chName, chUnit, reserved3) = struct.unpack(fmtF32channelinfo,s)
            self.channelNames.append(chName)
            self.channelUnits.append(chUnit)
    
    def __del__(self):
        self.close()

    def showInfo(self):
        print "Dataformat:\t\t%s %s" % (self.formatName, self.formatVer)
        print "Anzahl Kanaele:\t\t%i" % self.numChannels
        print "Anzahl Datenpunkte:\t%i" % self.numDatapoints
        print "Formatstring Sample:\t\"%s\"" % self.fmtSample
    
    def getData(self, start, length, stride=1, channelList = None):
        if(start+stride*length>(self.numDatapoints) or start<0 or length < 1 or stride<1):
            raise Exception, "Falsche Parameter fuer start, length und stride"
        # Keine channelList übergeben
        if (channelList == None):
            channelList = range(self.numChannels)
        channelList.sort()
        if (channelList[0]<0 or channelList[-1]>(self.numChannels-1)):
            raise Exception, "Kanalliste enthaelt falsche Kanalnummern.\min=%i, max=%%i" %(channelList[0],channelList[-1]) 
        #print channelList
        #Search starting position in file
        startpos = self.sbData + struct.calcsize(self.fmtSample)*start
        self.f.seek(startpos)    #Seek from beginning of file
        #Daten auslesen
        #rvData = [[] for i in channelList]
        rvData = n.zeros((length,len(channelList)),"d")
        for i in range(length):
            s = self.f.read(struct.calcsize(self.fmtSample))
            #print self.fmtSample
            ld_sp = struct.unpack(self.fmtSample, s)
            for j in range(len(channelList)):
                #ld_spSelected.append(ld_sp[j])
                rvData[i][j]=ld_sp[channelList[j]]
            #rvData.append(ld_spSelected)
            self.f.seek((stride-1)*struct.calcsize(self.fmtSample),1)    #seek relativ to position
            
        return rvData
    
    def getOverviewData(self, numPoints, channelList = None):
        if(numPoints > self.numDatapoints):
            raise Exception, "So viele Datenpunkte sind nicht in der Datei vorhanden."
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

    def moving_windows(self,channels = None, size = None, overlap = None):
        """Generator for iterating through the file, moving-window-technique."""
        if channels==None:
            channels=range(self.numChannels)
        if size==None:
            size=4096
        if overlap==None:
            overlap=0
        start = 0
        while start+size<self.numDatapoints:
            data = self.getData(start,size,1,channels)
            yield data, start
            start = start+size-overlap
    
    def moving_windows_num(self, size = None, overlap = None):
        """Helper functions. Calulates the number of windows needed to cover the whole file.
        This might be useful for allocating result-arrays prior to iteration"""
        if size==None:
            size=4096
        if overlap==None:
            overlap=0
        
        return (self.numDatapoints)/(size-overlap) #Not sure if it's right, need testing
        
    def getChannelNames(self):
        return self.channelNames
    
    def getChannelName(self,i):
        return self.channelNames[i]

    def close(self):
        """Close filehandle and finalize"""
        try:
            self.f.close()
        except:
            pass

 
class F32Writer:
    "Class for writing EEG-Files which are saved in the Bonn-F32-Format"
    def __init__(self, filename):
        # Initialisierung der Instanzvariablen
        self.fn = filename
        try:
            self.f = open(self.fn, "wb")
        except:
            raise IOError, "Die angegebene Datei konnte nicht geöffnet werden!"
        self.formatName = "Bonn F32 Data Format"
        self.formatVer = "1.0.00"
        self.sbChinfo = 1024 # Weil header 1024 bytes hat
        self.sbData = None
        self.numChannels = 0
        self.numDatapoints = 0
        self.samplRate = 1.0
        self.tStart = 0.0
        self.tEnd = 0.0
        self.reservedString = ""
        self.datatype = 4096
        self.samplesize = 17533
        self.minsample = 0.0
        self.maxsample = 0.0
        self.reserved2 = ""
        self.reserved3 = ""
        self.fmtSample = "" # Format-String für einzelnes Sample
        #(self.formatName,self.formatVer,self.sbChinfo,self.sbData,self.numChannels,self.numDatapoints,self.samplRate,self.tStart,self.tEnd,self.reservedString, self.datatype, self.samplesize, self.minsample, self.maxsample,self.reserved2) = struct.unpack(fmtF32header,s)
        #self.fmtSample = "= %if" % self.numChannels
        self.channelNames = None
        self.channelUnits = None
        self.dataToWrite = None
        
        self.headerWritten = False;
        #for i in range(0,self.numChannels):
        #    s = self.f.read(struct.calcsize(fmtF32channelinfo))
        #    (chName, chUnit, reserved3) = struct.unpack(fmtF32channelinfo,s)
        #    self.channelNames.append(chName)
        #    self.channelUnits.append(chUnit)
    
    def showInfo(self):
        print "Dataformat:\t\t%s %s" % (self.formatName, self.formatVer)
        print "Anzahl Kanaele:\t\t%i" % self.numChannels
        print "Anzahl Datenpunkte:\t%i" % self.numDatapoints
        print "Formatstring Sample:\t\"%s\"" % self.fmtSample
    
    def setChannelNames(self, cn):
        "Setzt die Anzahl der Kanäle"
        if not self.dataToWrite == None:
            assert self.dataToWrite.shape[1] == len(cn)
        self.channelNames = []
        self.channelUnits = []
        for i in cn:
            self.channelNames.append(str(i))
            self.channelUnits.append("")
            
        self.numChannels = len(self.channelNames)
        #print nc
        self.sbData = int(self.sbChinfo+self.numChannels*struct.calcsize(fmtF32channelinfo))
    
    def setData(self, data):
        """Nimmt ein numpy-Array. Index1:Zeitpunkt, Index2: Kanal"""    
        if not self.channelNames == None:
            assert data.shape[1] == len(self.channelNames)
        self.dataToWrite = data
        self.numDatapoints = int(data.shape[0])
        self.minsample = data.min()
        self.maxsample = data.max()
    #def setNumChannels(self, nc):
    #    "Setzt die Anzahl der Kanäle"
    #    if(headerWritten):
    #        headerWritten=false
    #    self.numChannels = nc 
    
    def writeHeader(self):
        """Writing the header"""
        if self.dataToWrite != None:
            assert len(self.dataToWrite.shape)==2
        #Preparation: chNames setzen, falls leer
        if self.channelNames == None:
            self.channelNames = ["%i"%(i+1) for i in range(self.dataToWrite.shape[1])]
            self.numChannels = len(self.channelNames)
        #Preparation: chUnits setzen, falls leer
        if self.channelUnits == None:
            self.channelUnits = ["" for i in self.channelNames]
        
        self.sbData = int(self.sbChinfo+self.numChannels*struct.calcsize(fmtF32channelinfo))
        
        if self.dataToWrite != None:
            assert len(self.channelNames)==self.dataToWrite.shape[1]    
        
        self.f.seek(0)
        s=struct.pack(fmtF32header, self.formatName,self.formatVer,int(self.sbChinfo),self.sbData,self.numChannels,self.numDatapoints,self.samplRate,self.tStart,self.tEnd,self.reservedString, self.datatype, self.samplesize, self.minsample, self.maxsample,self.reserved2)
        self.f.write(s)
        for i in range(len(self.channelNames)):
            #print len(self.channelNames), i
            s=struct.pack(fmtF32channelinfo, self.channelNames[i], self.channelUnits[i], self.reserved3)
            self.f.write(s)
        
        if(not (self.f.tell()==self.sbData)):
            raise Exception, "Irgendwie ist die Position falsch, %i != %i" % (self.f.tell(),self.sbData)
        self.headerWritten = True
        
    def writeData(self, channelNames = None, units = None):
        assert len(self.dataToWrite.shape)==2        
        assert self.dataToWrite.shape[1]==len(self.channelNames)
        
        self.fmtSample = "%if" % self.numChannels
        
        self.dataToWrite.astype("f").tofile(self.f)
        #for i in range(self.numDatapoints):
        #    for j in range(self.dataToWrite.shape[1]):
                #print j[i]
        #        self.f.write(struct.pack("f",self.dataToWrite[i,j]))
    
        
    def write(self):
        assert len(self.dataToWrite.shape)==2        
        
        self.writeHeader()
        assert self.dataToWrite.shape[1]==len(self.channelNames)
        self.writeData()

class F32WriterAdvanced(F32Writer):
    """F32writer for sequential writing of data. Everytime the appendData method is called, 
    the data are appended to the file and, when the writer is closed, the numDatappoints-fiels 
    in the header is updated."""
    def __init__(self,filename,cNames=None,numCs=None):
        F32Writer.__init__(self,filename)
        if cNames != None:
            assert len(cNames)>0, "Liste der Kanalnamen muss länger als 0 sein!"
            self.channelNames = cNames
            self.numChannels = len(self.channelNames)
            self.writeHeader()
        elif numCs!=0:
            assert numCs>0, "Anzahl Kanäle muss größer als 0 sein!"
            self.channelNames = ["%i"%(i+1) for i in range(self.numChannels)]
            self.writeHeader()
        #self.channelNames = ["%i"%(i+1) for i in range(self.dataToWrite.shape[1])]
        #self.numChannels = len(self.channelNames)
        #Preparation: chUnits setzen, falls leer
        #if self.channelUnits == None:
        #    self.channelUnits = ["" for i in self.channelNames]
        pass
    
    def __del__(self):
        self.writeNumDatapoints()
        #print "numDatapoints geschrieben..."
    
    def writeHeader(self):
        """Writing the header"""
        #Preparation: chNames setzen, falls leer
        if self.channelNames == None:
            if self.numChannels!=None:
                self.channelNames = ["%i"%(i+1) for i in range(self.numChannels)]
        else:
            self.numChannels = len(self.channelNames)
        #Preparation: chUnits setzen, falls leer
        if self.channelUnits == None:
            self.channelUnits = ["" for i in self.channelNames]
        
        self.sbData = int(self.sbChinfo+self.numChannels*struct.calcsize(fmtF32channelinfo))
        
        if self.dataToWrite != None:
            assert len(self.channelNames)==self.dataToWrite.shape[1]    
        
        self.f.seek(0)
        s=struct.pack(fmtF32header, self.formatName,self.formatVer,int(self.sbChinfo),self.sbData,self.numChannels,self.numDatapoints,self.samplRate,self.tStart,self.tEnd,self.reservedString, self.datatype, self.samplesize, self.minsample, self.maxsample,self.reserved2)
        self.f.write(s)
        for i in range(len(self.channelNames)):
            #print len(self.channelNames), i
            s=struct.pack(fmtF32channelinfo, self.channelNames[i], self.channelUnits[i], self.reserved3)
            self.f.write(s)
        
        if(not (self.f.tell()==self.sbData)):
            raise Exception, "Irgendwie ist die Position falsch, %i != %i" % (self.f.tell(),self.sbData)
        self.headerWritten = True
        
    def appendData(self,data):
        assert len(data.shape)==2, "Daten sind kein 2d-Array"
        assert self.headerWritten, "Header ist noch nicht geschrieben"
        assert data.shape[1] == self.numChannels, "Daten haben nicht die gleiche Anzahl Kanäle wie Datei!"
        
        #write data
        self.f.seek(0,2)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                #print j[i]
                self.f.write(struct.pack("f",data[i,j]))
        self.numDatapoints += data.shape[0]
    
    def writeNumDatapoints(self):
        self.f.seek(struct.calcsize(fmtPosOfNumDatapoints))
        self.f.write(struct.pack("=i", self.numDatapoints))
 
        
        
if __name__=='__main__':
    w=F32writer_advanced("test.f32",["c1","c2","c3"])
    w.appendData(n.zeros((10000,3),"d"))
    w.appendData(n.ones((10000,3),"d")*100)
    w.appendData(n.zeros((10000,3),"d"))
    w.appendData(n.ones((10000,3),"d")*100)
    del w
    
    

