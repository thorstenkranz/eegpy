# -*- coding: utf-8 -*-
#Dieses Modul dient zur Konvertierung von EEG-Daten im Brain Vision Exchange Format in das f32-format.

from f32 import F32Writer
import os.path
import sys

class BvHeaderParser:
    """This class is used to extract the informations of the Header-file that comes with the """
    def __init__(self,fn):
        self.myfile = open(fn,"r")
        self.contents = self.myfile.readlines()
        self.declareVars()
        self.parse()
        
            
    def declareVars(self):
        """Predeclares the instance-parameters to ensure they exist"""
        self.datafile = None
        self.markerfile = None
        self.dataformat = None
        self.dataorientation = None
        self.datatype = None
        self.numCh = None
        self.numDatapoints = None
        self.samplRate = None
        self.binaryformat = None
        self.channelNames = []
        
    def parse(self):
        """This method parses the header-file and saves the values as instance-parameters"""
        #Allgemeine Informationen auslesen
        for line in self.contents:
            if line.startswith("DataFile"):
                self.datafile = line.split("=")[1].replace("\n","").replace("\r","")
            if line.startswith("MarkerFile"):
                self.markerfile = line.split("=")[1].replace("\n","").replace("\r","")
            if line.startswith("DataFormat"):
                self.dataformat = line.split("=")[1].replace("\n","").replace("\r","")
            if line.startswith("DataOrientation"):
                self.dataorientation = line.split("=")[1].replace("\n","").replace("\r","")
            if line.startswith("DataType"):
                self.datatype = line.split("=")[1].replace("\n","").replace("\r","")
            if line.startswith("NumberOfChannels"):
                self.numCh = int(line.split("=")[1].replace("\n","").replace("\r",""))
            if line.startswith("DataPoints"):
                self.numDatapoints = int(line.split("=")[1].replace("\n","").replace("\r",""))
            if line.startswith("SamplingInterval"):
                self.samplRate = 1000000/int(line.split("=")[1].replace("\n","").replace("\r",""))
            if line.startswith("BinaryFormat"):
                self.binaryformat = line.split("=")[1].replace("\n","").replace("\r","")
        #Kanal-Namen extrahieren
        inChannelInfo=False
        #channelNames = []
        for line in self.contents:
            if line.startswith("["):
                inChannelInfo = False
            if line.startswith("[Channel Infos]"):
                inChannelInfo = True
            if inChannelInfo:
                if line.startswith("Ch"):
                    self.channelNames.append(line.split("=")[1].split(",")[0])
        #print self.channelNames
        self.numCh = len(self.channelNames)
                    
class BvToF32writer(F32Writer):
    def __init__(self, filename):
        F32Writer.__init__(self,filename)
        self.infile=None
        #self.parser = BvHeaderParser(self.fn)
    
    def setHeaderData(self, p):     
        """Sets all relevant header-data, retrieved from the parser-object p"""
        self.setChannelNames(p.channelNames)
        self.numDatapoints = p.numDatapoints
        self.samplRate = p.samplRate
        
    def writeData(self,infilename):
        """Takes the Name of the binary .dat file and just copies all of its content to the new f32-file"""
        if not self.headerWritten:
            return False
        
        try:
            self.infile = open(infilename, "rb")
        except:
            raise IOError, "Die angegebene Datei "+ infilename+" konnte nicht geöffnet werden!"
        #Für Statistik: Länge ermitteln
        self.infile.seek(0,2)
        fileLength = self.infile.tell()
        self.infile.seek(0)
        s=self.infile.read(8192)
        #print s
        nextMessageAt = 0.1
        while (len(s)>0):
            self.f.write(s)
            s=self.infile.read(8192)
            ratio = float(self.f.tell()) / float(fileLength)
            if ratio > nextMessageAt:
                #print ratio, fileLength
                print "%i Prozent ..." % int(nextMessageAt*100) 
                nextMessageAt +=0.1 
        self.f.close()
        self.infile.close()
 
class BvMarkerParser:   
    """This class is used to extract the informations of the marker-file that comes with the datafile"""
    def __init__(self,fn):
        self.fn = fn
        self.myfile = open(fn,"r")
        self.contents = self.myfile.readlines()
        self.markerDescriptions = []
        self.markers = []
        self.markersWithDescription = {}
        self.triggers =  {}
        self.allEvents = []
        self.parse()   
        self.myfile.close()
        
    def parse(self):
        """This method parses the marker-file and saves the values as adequate instance-parameters"""
        #Marker extrahieren
        inMarkerInfo=False
        #channelNames = []
        for line in self.contents:
            if line.startswith("["):
                inMarkerInfo = False
            if line.startswith("[Marker Infos]"):
                inMarkerInfo = True
            if inMarkerInfo:
                if line.startswith("Mk"):
                    #print line
                    ldAsList = line.split("=")[1].split(",")
                    self.markers.append( {"type":ldAsList[0] , "description":ldAsList[1], "pos":ldAsList[2]} )
                    #self.markers.append(line.split("=")[1].split(",")[0])
                    #if not self.markersWithDescription.has_key(ldAsList[1]):
                    #    self.markersWithDescription[ldAsList[1]] = []
                    #self.markersWithDescription[ldAsList[1]].append(self.markersWithDescription[ldAsList[0]])
        
        for m in self.markers:
            if not self.markersWithDescription.has_key(m["description"]):
                self.markersWithDescription[m["description"]] = []
            #if m["type"] == "Trigger":
            self.markersWithDescription[m["description"]].append(int(m["pos"]))
            
            if m["type"]=="Trigger":
                trigNum = int(m["description"][2:])
                if not self.triggers.has_key(trigNum):
                    self.triggers[trigNum] = []
                self.triggers[trigNum].append(int(m["pos"]))
            
            self.allEvents.append(int(m["pos"]))
                
        #print self.channelNames
        #self.numCh = len(self.channelNames)       