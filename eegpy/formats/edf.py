#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy 
from eegpy.misc import FATALERROR, debug
from eegpy.formats.iobase import EEG_file

try:
    import scipy
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


__biosig_docstring = """BIOSIG Toolbox for Python
    Copyright (C) 2005-2006 by Martin Hieden <martin.hieden@gmx.at> 
                 and Alois Schloegl <a.schloegl@ieee.org>

    $Id: biosig.py,v 1.1.1.1 2006/03/24 18:17:06 schloegl Exp $
    This function is part of the "BioSig for Python" repository 
    (biosig4python) at http://biosig.sf.net/ 


    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


FUNCTIONS
===================
    sopen
    sclose
    sread
    swrite
    seof
    srewind
    sseek
    stell


HEADER
===================
    TYPE:        type of file format
    VERSION:    GDF version number 
    FileName:    name of opened file
    HeadLen:    length of header in bytes
    NS:        number of channels
    SPR:        samples per block (when different sampling rates are used, this is the LCM(CHANNEL[..].SPR)
    NRec:        number of records/blocks -1 indicates length is unknown.
    Dur:        Duration of each block in seconds expressed in the fraction Dur[0]/Dur[1]
    SampleRate:    Sampling rate
    IPaddr:        IP address of recording device (if applicable)
    T0:        starttime of recording


    data        last data read
    ------------------------------
        block:        data block  
        size:        size {rows, columns} of data block    


    Patient:    Patient information
    -----------------------------------
        Name:        not recommended because of privacy protection 
        Id:        identification code as used in hospital 
        Weight:        weight in kilograms [kg] 0:unkown, 255: overflow 
        Height:        height in centimeter [cm] 0:unkown, 255: overflow 
        Birthday:    Birthday of Patient
        Age:        Age of Patient
        Headsize:    circumference, nasion-inion, left-right mastoid in millimeter
        Sex
        Handedness
        Smoking
        AlcoholAbuse
        DrugAbuse
        Medication
        Impairment
            Visual


    ID:        recording identification
    ----------------------------------------
        Technician
        Hospital
        Equipment:    identfies this software


    LOC:        location of recording according to RFC1876
    ----------------------------------------------------------
        VertPre
        HorizPre
        Size
        Version
        Latitude:    in degrees
        Longitude:    in degrees
        Altitude:    in metres


    ELEC:        position of electrodes; see also HDR.CHANNEL[k].XYZ
    -------------------------------------------------------------------
        REF:        XYZ position of reference electrode
        GND:        XYZ position of ground electrode


    EVENT:        EVENTTABLE
    --------------------------
        SampleRate:    for converting POS and DUR into seconds 
        N:        number of events
        TYP:        defined at http://cvs.sourceforge.net/viewcvs.py/biosig/biosig/t200/eventcodes.txt?view=markup
        POS:        starting position [in samples]
        DUR:        duration [in samples]
        CHN:        channel number; 0: all channels 


    FLAG:        flags
    ---------------------
        OVERFLOWDETECTION:    overflow & saturation detection 0: OFF, !=0 ON
        UCAL:            UnCalibration  0: scaling  !=0: NO scaling - raw data return 


    FILE:        File specific data
    ----------------------------------
        FID:        file handle 
        POS:        current reading/writing position in samples 
        OPEN:        0: closed, 1:read, 2: write
        LittleEndian:    not in use


    AS:        internal variables
    ----------------------------------
        PID:        patient identification
        RID:        recording identification 
        spb:        total samples per block
        bpb:        total bytes per block
        bi:        not in use
        Header1:    not in use
        rawdata:    raw data block 


    CHANNEL[k]:    channel specific data
    -------------------------------------
        Label:        Label of channel 
        Transducer:    transducer e.g. EEG: Ag-AgCl electrodes
        PhysDim:    physical dimension
        PhysDimCode:    code for physical dimension
        PreFilt:    pre-filtering
    
        LowPass:    lowpass filter
        HighPass:    high pass
        Notch:        notch filter
        XYZ:        electrode position
        Impedance:    in Ohm
    
        PhysMin:    physical minimum
        PhysMax:    physical maximum
        DigMin:        digital minimum
        DigMax:        digital maximum
    
        GDFTYP:        data type
        SPR:        samples per record (block)
        bpr:        bytes per record (block)
    
        OnOff:        1: include, 0: exclude in sread
        Cal:        gain factor 
        Off:        bias"""


try:
################
# MODUL IMPORT #
################
    
    try:
        import scipy
        from scipy import NaN
    except ImportError:
        raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    import math
    import struct
    import datetime
    
###############
# GLOBAL DATA #
###############
    
    # not used so far:    
    #    __FileFormat = ('ACQ', 'BKR', 'BDF', 'CNT', 'DEMG', 'EDF', 'EVENT', 'FLAC', 'GDF', 'MFER', 'NEX1', 'PLEXON') 
    
    __HANDEDNESS = ('Unknown', 'Right', 'Left', 'Equal')
    __GENDER = ('Unknown', 'Male', 'Female')
    __SCALE = ('Unknown', 'No', 'Yes', 'Corrected')
    
    __GDFTYP_NAME = [None]
    try:
        from scipy import int8
        from scipy import uint8
        from scipy import int16
        from scipy import uint16
        from scipy import int32
        from scipy import uint32
        __GDFTYP_NAME.append(int8)
        __GDFTYP_NAME.append(uint8)
        __GDFTYP_NAME.append(int16)
        __GDFTYP_NAME.append(uint16)
        __GDFTYP_NAME.append(int32)
        __GDFTYP_NAME.append(uint32)
    except NameError:
        raise FATALERROR('Standard datatypes are not supported by this SciPy version!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    try:
        from scipy import int64 as int64
        __GDFTYP_NAME.append(int64)
    except NameError:
        __GDFTYP_NAME.append(None)
        print 'int64 is not supported by this SciPy version.'
        print 'Please visit www.scipy.org or numeric.scipy.org for more information.'
    try:
        from scipy import uint64 as uint64
        __GDFTYP_NAME.append(uint64)
    except NameError:
        __GDFTYP_NAME.append(None)
        print 'uint64 is not supported by this SciPy version.'
        print 'Please visit www.scipy.org or numeric.scipy.org for more information.'
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    __GDFTYP_NAME.append(None)
    try:
        from scipy import float32 as float32
        __GDFTYP_NAME.append(float32)
    except NameError:
        raise FATALERROR('Standard datatypes are not supported by this SciPy version!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')
    try:
        from scipy import float64 as float64
        __GDFTYP_NAME.append(float64)
    except NameError:
        __GDFTYP_NAME.append(None)
        print 'float64 is not supported by this SciPy version.'
        print 'Please visit www.scipy.org or numeric.scipy.org for more information.'
    #try:
    #    from scipy import float128
    #    __GDFTYP_NAME.append(float128)
    #except NameError:
    #    __GDFTYP_NAME.append(None)
    #    print 'float128 is not supported by this SciPy version.'
    #    print 'Please visit www.scipy.org or numeric.scipy.org for more information.'
    
    __GDFTYP_BYTE = (1, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 0, 0, 0, 0, 0, 4, 8, 16)
    
    
#########################
# FATAL ERROR EXCEPTION #
#########################

except FATALERROR, e:
    print ''
    print 'FATAL ERROR:'
    print '============'
    print e.value
    print '============'
    print ''



###########################
# Class Implementation    #
###########################


class EDFReader(EEG_file):
    """Functionality of biosig4python in a class. Now will try to make it usable like an array."""
    HDR = None
    
    def __init__(self, FileName, MODE = 'r', HDR = None):
        """input:  string FileName, [r,w] MODE, HDR_TYPE HDR
        output: HDR_TYPE
        Opens a file for reading.
        Writing is not yet implemented.
        Supported dataformats: BDF, EDF, GDF, GDF2"""
        self.HDR = HDR
        if self.HDR == None:
            self.HDR = HDR_TYPE()
            
        try:
            
            # READ
            if MODE == 'r':
                self.HDR.FILE.FID = open(FileName, 'rb')
                version = self.HDR.FILE.FID.read(8)
                
                # BDF & EDF
                if version[1:] == 'BIOSEMI' or version[:] == '0       ':
                    if version[1:] == 'BIOSEMI':
                        self.HDR.TYPE = 'BDF'
                        self.HDR.VERSION = -1
                    else:
                        self.HDR.TYPE = 'EDF'
                        self.HDR.VERSION = 0
                    
                    #FIXED HEADER
                    self.HDR.AS.PID = self.HDR.FILE.FID.read(80)
                    self.HDR.AS.RID = self.HDR.FILE.FID.read(80)
                    tm = self.HDR.FILE.FID.read(16)
                    if int(tm[6:8]) < 85:
                        tm = '20' + tm[6:8] + tm[3:5] + tm[0:2] + tm[8:10] + tm[11:13] + tm[14:16] + '00'
                    else:
                        tm = '19' + tm[6:8] + tm[3:5] + tm[0:2] + tm[8:10] + tm[11:13] + tm[14:16] + '00'
                    self.HDR.T0 = self.__gdf_time2py_time(tm)
                    self.HDR.HeadLen = int(self.HDR.FILE.FID.read(8))
                    reserved = self.HDR.FILE.FID.read(44)    #44bytes reserved
                    if reserved[0:4] == 'EDF+':
                        pid = self.HDR.AS.PID.split(' ')    #PID split up in 'patient identification code', 'sex', 'birthday', 'patient name' and 'patient classification'
                        if len(pid) >= 4:
                            self.HDR.Patient.Id = pid[0]
                            if pid[1][0].lower() == 'f':
                                self.HDR.Patient.Sex = self.__GENDER[2]
                            else:
                                self.HDR.Patient.Sex = self.__GENDER[1]
                            bday = pid[2].split('-')
                            if len(bday) >= 3:
                                if bday[1].lower() == 'jan':
                                    bday[1] = '01'
                                elif bday[1].lower() == 'feb':
                                    bday[1] = '02'
                                elif bday[1].lower() == 'mar':
                                    bday[1] = '03'
                                elif bday[1].lower() == 'apr':
                                    bday[1] = '04'
                                elif bday[1].lower() == 'may':
                                    bday[1] = '05'
                                elif bday[1].lower() == 'jun':
                                    bday[1] = '06'
                                elif bday[1].lower() == 'jul':
                                    bday[1] = '07'
                                elif bday[1].lower() == 'aug':
                                    bday[1] = '08'
                                elif bday[1].lower() == 'sep':
                                    bday[1] = '09'
                                elif bday[1].lower() == 'oct':
                                    bday[1] = '10'
                                elif bday[1].lower() == 'nov':
                                    bday[1] = '11'
                                elif bday[1].lower() == 'dec':
                                    bday[1] = '12'
                                self.HDR.Patient.Birthday = self.__gdf_time2py_time(bday[2] + bday[1] + bday[0] + '12000000')
                            self.HDR.Patient.Name = pid[3]
                    
                    self.HDR.NRec = int(self.HDR.FILE.FID.read(8))
                    self.HDR.Dur = scipy.array([float(self.HDR.FILE.FID.read(8)),1.])
                    self.HDR.NS = int(self.HDR.FILE.FID.read(4))
                    
                    #VARIABLE HEADER
                    self.HDR.SPR = 1
                    vh = self.HDR.FILE.FID.read(self.HDR.HeadLen-256)
                    for k in range(self.HDR.NS):
                        self.HDR.CHANNEL.append(CHANNEL_TYPE())
                        i = k*16
                        self.HDR.CHANNEL[k].Label = vh[i:i+16]
                        i = i + (self.HDR.NS-k)*16 + k*80
                        self.HDR.CHANNEL[k].Transducer = vh[i:i+80]
                        i = i + (self.HDR.NS-k)*80 + k*8
                        self.HDR.CHANNEL[k].PhysDim = vh[i:i+8]
                        i = i + (self.HDR.NS-k)*8 + k*8
                        self.HDR.CHANNEL[k].PhysMin = float(vh[i:i+8])
                        i = i + (self.HDR.NS-k)*8 + k*8
                        self.HDR.CHANNEL[k].PhysMax = float(vh[i:i+8])
                        i = i + (self.HDR.NS-k)*8 + k*8
                        self.HDR.CHANNEL[k].DigMin = float(vh[i:i+8])
                        i = i + (self.HDR.NS-k)*8 + k*8
                        self.HDR.CHANNEL[k].DigMax = float(vh[i:i+8])
                        i = i + (self.HDR.NS-k)*8 + k*80
                        self.HDR.CHANNEL[k].PreFilt = vh[i:i+80]
                        i = i + (self.HDR.NS-k)*80 + k*8
                        self.HDR.CHANNEL[k].SPR = long(vh[i:i+8])
                        self.HDR.SPR = self.__lcm(self.HDR.SPR, self.HDR.CHANNEL[k].SPR)    # least common SPR
                        self.HDR.AS.spb = self.HDR.AS.spb + self.HDR.CHANNEL[k].SPR
                        i = i + (self.HDR.NS-k)*8 + k*32
                        if self.HDR.TYPE == 'BDF':
                            self.HDR.CHANNEL[k].GDFTYP = 255 + 24
                            self.HDR.CHANNEL[k].bpr = 3 * self.HDR.CHANNEL[k].SPR
                            self.HDR.AS.bpb = self.HDR.AS.bpb + self.HDR.CHANNEL[k].bpr
                        else:
                            self.HDR.CHANNEL[k].GDFTYP = 3
                            self.HDR.CHANNEL[k].bpr = self.__GDFTYP_BYTE[self.HDR.CHANNEL[k].GDFTYP] * self.HDR.CHANNEL[k].SPR
                            self.HDR.AS.bpb = self.HDR.AS.bpb + self.HDR.CHANNEL[k].bpr
                        
                        #reserved
                        i = i +(self.HDR.NS-k)*32
                            
                        self.HDR.CHANNEL[k].Cal = (self.HDR.CHANNEL[k].PhysMax - self.HDR.CHANNEL[k].PhysMin) / (self.HDR.CHANNEL[k].DigMax - self.HDR.CHANNEL[k].DigMin)
                        self.HDR.CHANNEL[k].Off = self.HDR.CHANNEL[k].PhysMin - self.HDR.CHANNEL[k].Cal * self.HDR.CHANNEL[k].DigMin
                            
                        if i != len(vh):
                            raise FATALERROR('Error reading variable Header!\nSignal: ' + str(k))
                    
                    #Finished reading Header
                    self.HDR.FLAG.OVERFLOWDETECTION = 0    #EDF does not support automated overflow and saturation detection
                    self.HDR.FILE.OPEN = 1
                    self.HDR.FILE.FID.seek(self.HDR.HeadLen)
                    self.HDR.FILE.POS = 0
                    self.HDR.SampleRate = float(self.HDR.SPR) * self.HDR.Dur[1] / self.HDR.Dur[0]
                        
                    self._shape = (self.HDR.NRec*self.HDR.SPR, self.HDR.NS)
                    if debug:
                        print 'Finished reading header'

                
                # GDF
                # fromstring(...).tolist()[0] is necessary because of wrong type (..._arrtype) returned by fromstring(...)[0]
                elif version[:3] == 'GDF':
                    HDR.TYPE = 'GDF'
                    HDR.VERSION = float(version[4:])
                    
                    # GDF 1.x
                    if HDR.VERSION < 1.9:
                        #FIXED HEADER
                        HDR.AS.PID = HDR.FILE.FID.read(80)
                        pid = HDR.AS.PID.split(' ', 2)    #PID split up in 'patient identification code', 'patient name' and 'patient classification'
                        if len(pid) >= 2:
                            self.HDR.Patient.Id = pid[0]
                            self.HDR.Patient.Name = pid[1]
                        
                        self.HDR.AS.RID = self.HDR.FILE.FID.read(80)
                        self.HDR.T0 = self.__gdf_time2py_time(self.HDR.FILE.FID.read(16))
                        self.HDR.HeadLen = scipy.fromstring(self.HDR.FILE.FID.read(8), int64).tolist()[0]
                        self.HDR.ID.Equipment = scipy.fromstring(self.HDR.FILE.FID.read(8), uint8)
                        self.HDR.ID.Hospital = scipy.fromstring(self.HDR.FILE.FID.read(8), uint8)
                        self.HDR.ID.Technician = scipy.fromstring(self.HDR.FILE.FID.read(8), uint8)
                        self.HDR.FILE.FID.seek(20,1)    #20bytes reserved
                        self.HDR.NRec = scipy.fromstring(self.HDR.FILE.FID.read(8), int64).tolist()[0]
                        self.HDR.Dur = scipy.fromstring(self.HDR.FILE.FID.read(8), uint32)
                        self.HDR.NS = scipy.fromstring(self.HDR.FILE.FID.read(4), uint32).tolist()[0]
                        
                        #VARIABLE HEADER
                        self.HDR.SPR = 1
                        vh = self.HDR.FILE.FID.read(self.HDR.HeadLen-256)
                        for k in range(self.HDR.NS):
                            self.HDR.CHANNEL.append(CHANNEL_TYPE())
                            i = k*16
                            self.HDR.CHANNEL[k].Label = vh[i:i+16]
                            i = i + (self.HDR.NS-k)*16 + k*80
                            self.HDR.CHANNEL[k].Transducer = vh[i:i+80]
                            i = i + (self.HDR.NS-k)*80 + k*8
                            self.HDR.CHANNEL[k].PhysDim = vh[i:i+8]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].PhysMin = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].PhysMax = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].DigMin = scipy.fromstring(vh[i:i+8], int64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].DigMax = scipy.fromstring(vh[i:i+8], int64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*80
                            self.HDR.CHANNEL[k].PreFilt = vh[i:i+80]
                            i = i + (self.HDR.NS-k)*80 + k*4
                            self.HDR.CHANNEL[k].SPR = scipy.fromstring(vh[i:i+4], uint32).tolist()[0]
                            self.HDR.SPR = self.__lcm(self.HDR.SPR, self.HDR.CHANNEL[k].SPR)    # least common SPR
                            self.HDR.AS.spb = self.HDR.AS.spb + self.HDR.CHANNEL[k].SPR
                            i = i + (self.HDR.NS-k)*4 + k*4
                            self.HDR.CHANNEL[k].GDFTYP = scipy.fromstring(vh[i:i+4], uint32).tolist()[0]
                            self.HDR.CHANNEL[k].bpr = self.__GDFTYP_BYTE[self.HDR.CHANNEL[k].GDFTYP] * self.HDR.CHANNEL[k].SPR
                            self.HDR.AS.bpb = self.HDR.AS.bpb + self.HDR.CHANNEL[k].bpr
                            i = i + (self.HDR.NS-k)*4 + k*32
                            #reserved
                            i = i +(self.HDR.NS-k)*32
                            
                            self.HDR.CHANNEL[k].Cal = (self.HDR.CHANNEL[k].PhysMax - self.HDR.CHANNEL[k].PhysMin) / (self.HDR.CHANNEL[k].DigMax - self.HDR.CHANNEL[k].DigMin)
                            self.HDR.CHANNEL[k].Off = self.HDR.CHANNEL[k].PhysMin - self.HDR.CHANNEL[k].Cal * self.HDR.CHANNEL[k].DigMin
                            
                            if i != len(vh):
                                raise FATALERROR('Error reading variable Header!\nSignal: ' + str(k))
                        
                        #EVENT TABLE
                        etp = self.HDR.HeadLen + self.HDR.NRec*self.HDR.AS.bpb
                        self.HDR.FILE.FID.seek(etp)
                        etmode = self.HDR.FILE.FID.read(1)
                        if etmode != '':
                            etmode = scipy.fromstring(etmode, uint8).tolist()[0]
                            sr = scipy.fromstring(self.HDR.FILE.FID.read(3), uint8)
                            self.HDR.EVENT.SampleRate = sr[0]
                            for i in range(1,len(sr)):
                                self.HDR.EVENT.SampleRate = self.HDR.EVENT.SampleRate + sr[i]*2**i
                            
                            self.HDR.EVENT.N = scipy.fromstring(self.HDR.FILE.FID.read(4), uint32).tolist()[0]
                            self.HDR.EVENT.POS = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*4), uint32)
                            self.HDR.EVENT.TYP = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*2), uint16)
                            
                            if etmode == 3:
                                self.HDR.EVENT.CHN = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*2), uint16)
                                self.HDR.EVENT.DUR = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*4), uint32)
                        
                        #Finished reading Header
                        self.HDR.FILE.OPEN = 1
                        self.HDR.FILE.FID.seek(self.HDR.HeadLen)
                        self.HDR.FILE.POS = 0
                        self.HDR.SampleRate = float(self.HDR.SPR) * self.HDR.Dur[1] / self.HDR.Dur[0]
                        if self.HDR.EVENT.SampleRate == 0:
                            self.HDR.EVENT.SampleRate = self.HDR.SampleRate
                        
                        self._shape = (self.HDR.NRec*self.HDR.SPR, self.HDR.NS)
                        if debug:
                            print 'Finished reading header'
                    
                    # GDF 2.0
                    else:
                        #FIXED HEADER
                        self.HDR.AS.PID = self.HDR.FILE.FID.read(66)
                        pid = self.HDR.AS.PID.split(' ', 2)    #PID split up in 'patient identification code', 'patient name' and 'patient classification'
                        if len(pid) >= 2:
                            self.HDR.Patient.Id = pid[0]
                            self.HDR.Patient.Name = pid[1]
                        
                        self.HDR.FILE.FID.seek(10,1)    #10bytes reserved
                        sadm = scipy.fromstring(self.HDR.FILE.FID.read(1), uint8).tolist()[0]    #Smoking / Alcohol abuse / drug abuse / medication
                        self.HDR.Patient.Smoking = self.__SCALE[sadm%4]
                        self.HDR.Patient.AlcoholAbuse = self.__SCALE[(sadm>>2)%4]
                        self.HDR.Patient.DrugAbuse = self.__SCALE[(sadm>>4)%4]
                        self.HDR.Patient.Medication = self.__SCALE[(sadm>>6)%4]
                        
                        self.HDR.Patient.Weight = scipy.fromstring(self.HDR.FILE.FID.read(1), uint8).tolist()[0]
                        if self.HDR.Patient.Weight == 0 or self.HDR.Patient.Weight == 255:
                            self.HDR.Patient.Weight = NaN
                        self.HDR.Patient.Height = scipy.fromstring(self.HDR.FILE.FID.read(1), uint8).tolist()[0]
                        if self.HDR.Patient.Height == 0 or self.HDR.Patient.Height == 255:
                            self.HDR.Patient.Height = NaN
                        ghi = scipy.fromstring(self.HDR.FILE.FID.read(1), uint8).tolist()[0]    #Gender / Handedness / Visual Impairment
                        self.HDR.Patient.Sex = self.__GENDER[ghi%4]
                        self.HDR.Patient.Handedness = self.__HANDEDNESS[(ghi>>2)%4]
                        self.HDR.Patient.Impairment.Visual = self.__SCALE[(ghi>>4)%4]
                        
                        self.HDR.AS.RID = self.HDR.FILE.FID.read(64)
                        rl = self.HDR.FILE.FID.read(16)    #Recording location (Lat, Long, Alt)
                        vhsv = scipy.fromstring(rl[0:4], uint8)
                        if vhsv[3] == 0:
                            self.HDR.LOC.VertPre = 10 * int(vhsv[0]>>4) + int(vhsv[0]%16)
                            self.HDR.LOC.HorzPre = 10 * int(vhsv[1]>>4) + int(vhsv[1]%16)
                            self.HDR.LOC.Size = 10 * int(vhsv[2]>>4) + int(vhsv[2]%16)
                        else:
                            self.HDR.LOC.VertPre = 29
                            self.HDR.LOC.HorizPre = 29
                            self.HDR.LOC.Size = 29
                        self.HDR.LOC.Version = 0
                        self.HDR.LOC.Latitude = float(scipy.fromstring(rl[4:8], uint32).tolist()[0]) / 3600000
                        self.HDR.LOC.Longitude = float(scipy.fromstring(rl[8:12], uint32).tolist()[0]) / 3600000
                        self.HDR.LOC.Altitude = float(scipy.fromstring(rl[12:16], int32).tolist()[0]) / 100
                        
                        self.HDR.T0 = self.__gdf2_time2py_time(scipy.fromstring(self.HDR.FILE.FID.read(8), uint64).tolist()[0])
                        self.HDR.Patient.Birthday = self.__gdf2_time2py_time(scipy.fromstring(self.HDR.FILE.FID.read(8), uint64).tolist()[0])
                        if self.HDR.Patient.Birthday != datetime.datetime(1,1,1,0,0):
                            today = datetime.datetime.today()
                            self.HDR.Patient.Age = today.year - self.HDR.Patient.Birthday.year
                            today = today.replace(year=self.HDR.Patient.Birthday.year)
                            if today < self.HDR.Patient.Birthday:
                                self.HDR.Patient.Age -= 1
                        else:
                            Age = NaN
                        
                        self.HDR.HeadLen = scipy.fromstring(self.HDR.FILE.FID.read(2), uint16).tolist()[0]*256
                        self.HDR.FILE.FID.seek(6,1)    #6bytes reserved
                        self.HDR.ID.Equipment = scipy.fromstring(self.HDR.FILE.FID.read(8), uint8)
                        self.HDR.IPaddr = scipy.fromstring(self.HDR.FILE.FID.read(6), uint8)
                        self.HDR.Patient.Headsize = scipy.fromstring(self.HDR.FILE.FID.read(6), uint16)
                        self.HDR.Patient.Headsize = scipy.asarray(self.HDR.Patient.Headsize, float32)
                        self.HDR.Patient.Headsize = scipy.ma.masked_array(self.HDR.Patient.Headsize, scipy.equal(self.HDR.Patient.Headsize, 0), NaN).filled()
                        
                        self.HDR.ELEC.REF = scipy.fromstring(self.HDR.FILE.FID.read(12), float32)
                        self.HDR.ELEC.GND = scipy.fromstring(self.HDR.FILE.FID.read(12), float32)
                        self.HDR.NRec = scipy.fromstring(self.HDR.FILE.FID.read(8), int64).tolist()[0]
                        self.HDR.Dur = scipy.fromstring(self.HDR.FILE.FID.read(8), uint32)
                        self.HDR.NS = scipy.fromstring(self.HDR.FILE.FID.read(2), uint16).tolist()[0]
                        self.HDR.FILE.FID.seek(2,1)    #2bytes reserved
                        
                        #VARIABLE HEADER
                        self.HDR.SPR = 1
                        vh = self.HDR.FILE.FID.read(self.HDR.HeadLen-256)
                        for k in range(self.HDR.NS):
                            self.HDR.CHANNEL.append(CHANNEL_TYPE())
                            i = k*16
                            self.HDR.CHANNEL[k].Label = vh[i:i+16]
                            i = i + (self.HDR.NS-k)*16 + k*80
                            self.HDR.CHANNEL[k].Transducer = vh[i:i+80]
                            i = i + (self.HDR.NS-k)*80 + k*6
                            self.HDR.CHANNEL[k].PhysDim = vh[i:i+6]
                            i = i + (self.HDR.NS-k)*6 + k*2
                            self.HDR.CHANNEL[k].PhysDimCode = scipy.fromstring(vh[i:i+2], uint16).tolist()[0]
                            i = i + (self.HDR.NS-k)*2 + k*8
                            self.HDR.CHANNEL[k].PhysMin = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].PhysMax = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].DigMin = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*8
                            self.HDR.CHANNEL[k].DigMax = scipy.fromstring(vh[i:i+8], float64).tolist()[0]
                            i = i + (self.HDR.NS-k)*8 + k*68
                            self.HDR.CHANNEL[k].PreFilt = vh[i:i+68]
                            i = i + (self.HDR.NS-k)*68 + k*4
                            self.HDR.CHANNEL[k].LowPass = scipy.fromstring(vh[i:i+4], float32).tolist()[0]
                            i = i + (self.HDR.NS-k)*4 + k*4
                            self.HDR.CHANNEL[k].HighPass = scipy.fromstring(vh[i:i+4], float32).tolist()[0]
                            i = i + (self.HDR.NS-k)*4 + k*4
                            self.HDR.CHANNEL[k].Notch = scipy.fromstring(vh[i:i+4], float32).tolist()[0]
                            i = i + (self.HDR.NS-k)*4 + k*4
                            self.HDR.CHANNEL[k].SPR = scipy.fromstring(vh[i:i+4], uint32).tolist()[0]
                            self.HDR.SPR = self.__lcm(self.HDR.SPR, self.HDR.CHANNEL[k].SPR)    # least common SPR
                            self.HDR.AS.spb = self.HDR.AS.spb + self.HDR.CHANNEL[k].SPR
                            i = i + (self.HDR.NS-k)*4 + k*4
                            self.HDR.CHANNEL[k].GDFTYP = scipy.fromstring(vh[i:i+4], uint32).tolist()[0]
                            self.HDR.CHANNEL[k].bpr = self.__GDFTYP_BYTE[self.HDR.CHANNEL[k].GDFTYP] * self.HDR.CHANNEL[k].SPR
                            self.HDR.AS.bpb = self.HDR.AS.bpb + self.HDR.CHANNEL[k].bpr
                            i = i + (self.HDR.NS-k)*4 + k*12
                            self.HDR.CHANNEL[k].XYZ = scipy.fromstring(vh[i:i+12], float32)
                            i = i + (self.HDR.NS-k)*12 + k*1
                            self.HDR.CHANNEL[k].Impedance= pow(2,float(scipy.fromstring(vh[i:i+1], uint8)[0])/8)
                            i = i + (self.HDR.NS-k)*1 + k*19
                            #reserved
                            i = i +(self.HDR.NS-k)*19
                            
                            self.HDR.CHANNEL[k].Cal = (self.HDR.CHANNEL[k].PhysMax - self.HDR.CHANNEL[k].PhysMin) / (self.HDR.CHANNEL[k].DigMax - self.HDR.CHANNEL[k].DigMin)
                            self.HDR.CHANNEL[k].Off = self.HDR.CHANNEL[k].PhysMin - self.HDR.CHANNEL[k].Cal * self.HDR.CHANNEL[k].DigMin
                            
                            if i != len(vh):
                                print 'Error reading variable Header!'
                                break
                        
                        #EVENT TABLE
                        etp = self.HDR.HeadLen + self.HDR.NRec*self.HDR.AS.bpb
                        self.HDR.FILE.FID.seek(etp)
                        etmode = self.HDR.FILE.FID.read(1)
                        if etmode != '':
                            etmode = scipy.fromstring(etmode, uint8).tolist()[0]
                            
                            if self.HDR.VERSION < 1.94:
                                sr = scipy.fromstring(self.HDR.FILE.FID.read(3), uint8)
                                self.HDR.EVENT.SampleRate = sr[0]
                                for i in range(1,len(sr)):
                                    self.HDR.EVENT.SampleRate = self.HDR.EVENT.SampleRate + sr[i]*2**i
                                self.HDR.EVENT.N = scipy.fromstring(self.HDR.FILE.FID.read(4), uint32).tolist()[0]
                            else:
                                ne = scipy.fromstring(self.HDR.FILE.FID.read(3), uint8)
                                self.HDR.EVENT.N = ne[0]
                                for i in range(1,len(ne)):
                                    self.HDR.EVENT.N = self.HDR.EVENT.N + ne[i]*2**i
                                self.HDR.EVENT.SampleRate = scipy.fromstring(self.HDR.FILE.FID.read(4), float32).tolist()[0]
    
                            self.HDR.EVENT.POS = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*4), uint32)
                            self.HDR.EVENT.TYP = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*2), uint16)
                            
                            if etmode == 3:
                                self.HDR.EVENT.CHN = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*2), uint16)
                                self.HDR.EVENT.DUR = scipy.fromstring(self.HDR.FILE.FID.read(self.HDR.EVENT.N*4), uint32)
                        
                        #Finished reading Header
                        self.HDR.FILE.OPEN = 1
                        self.HDR.FILE.FID.seek(self.HDR.HeadLen)
                        self.HDR.FILE.POS = 0
                        self.HDR.SampleRate = float(self.HDR.SPR) * self.HDR.Dur[1] / self.HDR.Dur[0]
                        if self.HDR.EVENT.SampleRate == 0:
                            self.HDR.EVENT.SampleRate = self.HDR.SampleRate
                        
                        self._shape = (self.HDR.NRec*self.HDR.SPR, self.HDR.NS)
                        if debug:
                            print 'Finished reading header'
                
                # other File Formats
                else:
                    print 'This file format is not implemented yet!'
                
            # WRITE
            elif MODE == 'w':
                print 'Writing is not implemented yet!'
            
        except FATALERROR, e:
            print ''
            print 'FATAL ERROR:'
            print '============'
            print e.value
            print '============'
            print ''
        # End of __init__
    
    def close(self):
        """Closes the according file.
        
        output: [-1, 0]
        Returns 0 if successfull, -1 otherwise."""
        
        if self.HDR.FILE.OPEN != 0:
            self.HDR.FILE.FID.close()
            self.HDR.FILE.FID = 0
            self.HDR.FILE.OPEN = 0
            return 0
        return -1
        # End of SCLOSE
        
    #functions needed to implement the sigfile-Interface
    
    def getData(self, start, length, stride=1, channelList = None):
        """Wrapper for function _sread. Conform to the sigfile-interface.
        
        The perfomance of this function is poor for large strides. 
        It will be rewritten later, to read directly from the file."""
        
        start=int(start)
        length=int(length)
        stride=int(stride)
        if (channelList == None):
            channelList = range(self.HDR.NS)
        channelList.sort()    
        
        #Assertions...
        assert not (start+stride*length>(self.numDatapoints) or start<0 or length < 1 or stride<1), "Wrong parameters for start, length und stride"
        assert not (channelList[0]<0 or channelList[-1]>(self.numChannels-1)), "Channel-list contains illegal channel-numbers: min=%i, max=%i" %(channelList[0],channelList[-1])
        
        d={} # This is a dictionary. It's keys will correspond to the blocks that have to be read, the values are lists containing the indices of the Datapoints to take in each of these blocks.
        
        for i in range(length):
            x=start+i*stride
            if not d.has_key(x/self.HDR.SPR):
                d[x/self.HDR.SPR] = []
            d[x/self.HDR.SPR].append(x%self.HDR.SPR)
        ks = d.keys()
        ks.sort()
        
        rv = scipy.ones((length,len(channelList)), float64) * NaN    # right sized output matrix filled with NaN's
        idx=0
        for k in ks:
            tmpData = self._sread(1,k,channelList)
            for i, pos in enumerate(d[k]):
                #print k, i, pos, idx
                rv[idx,:] = tmpData[pos,:]
                idx += 1
                   
        #print self.HDR.NS, self.HDR.SPR 
        #print d
        return rv
    
    def getOverviewData(self, numPoints=1000, channelList = None):
        assert not (numPoints > self.numDatapoints), "Number of requested datapoints exceeds number of datapoints in file."
        # Keine channelList übergeben
        if (channelList == None):
            channelList = range(self.numChannels)
        channelList.sort()
        #Waehle geeignete Parameter und rufe damit getData auf
        #start=0, length=numPoints, stride=?
        stride = self.numDatapoints / numPoints
        return self.getData(0,numPoints,stride,channelList)
    
    def getChannelNames(self):
        return [self.HDR.CHANNEL[i].Label for i in range(self.numChannels)]
    
    def _sread(self, length = -1, start = 0, channelList = None):
        """input: HDR_TYPE HDR, int length, int start
        output: array
        Reads and returns data according to length and start.
        length is the number of blocks to read. Use -1 to read all blocks until the end.
        start is the block to begin reading with.
        Use HDR.CHANNEL[k].OnOff to exclude single channels from reading."""
        
        raw = ''            # const input string
        channel = []            # const channels to read
        ns = 0                # const number of signals to read
        
        i = 0                # index for writing in output matrix
        bstart = 0            # start index (in bytes) of next value on the input string
        bend = 0            # length (in bytes) of the next value on the input string and next bstart
        row = 0                # row in output matrix
        column = 0            # column in output matrix
        ch = 0                # actual channel
        value = 0            # actual value to insert into output matrix
        
        if channelList == None:
            for ch in range(self.HDR.NS):
                if self.HDR.CHANNEL[ch].OnOff != 0:
                    channel.append(ch)
        else:
            channel = channelList
            
        ns = len(channel)
        
        seek = 0
        if start >= 0:
            seek = self._sseek(start, -1)
        
        if seek == 0:
            if length == -1:
                length = self.HDR.NRec    # read all blocks
            length = max(min(length, self.HDR.NRec - self.HDR.FILE.POS),0)        # number of blocks to read
            self.HDR.data.block = scipy.ones((self.HDR.SPR * length,ns), float64) * NaN    # right sized output matrix filled with NaN's
            raw = self.HDR.FILE.FID.read(self.HDR.AS.bpb * length)            # read input string from file
            
            for i in range(scipy.size(self.HDR.data.block)):
                row = i / (ns * self.HDR.SPR) * self.HDR.SPR + i % self.HDR.SPR
                column = (i / self.HDR.SPR) % ns
                ch = channel[column]
                
                # OnOff calculations
                # not all channels are being read and it's the first row in a block
                if (ns != self.HDR.NS) and (row % self.HDR.SPR == 0):
                    # first channel is missing
                    if (column == 0) and (channel[0] != 0):
                        # first block -> take only the leading missing channels in account
                        # other blocks -> take trailing and leading missing channels in account
                        if i != 0:
                            # for all channels being left out between channel[-1] (last wanted channel) and self.HDR.NS (last possible channel)
                            for leftout in range(channel[-1] + 1, self.HDR.NS):
                                bstart += self.HDR.CHANNEL[leftout].bpr
                        # for all channels being left out between 0 and channel[column]
                        for leftout in range(ch):
                            bstart += self.HDR.CHANNEL[leftout].bpr
                    
                    # channels in between are missing
                    elif ch != channel[column - 1] + 1:
                        # for all channels being left out between channel[column - 1] and channel[column]
                        for leftout in range(channel[column - 1] + 1, ch):
                            bstart += self.HDR.CHANNEL[leftout].bpr
                
                
                # reading from input string, calculating right value and inserting in output matrix
                if (row * self.HDR.CHANNEL[ch].SPR) % self.HDR.SPR == 0:
                    # reading and calculating of value
                    # BDF
                    if self.HDR.CHANNEL[ch].GDFTYP == 255 + 24:
                        bend = bstart + 3
                        value = scipy.fromstring(raw[bstart:bend], uint8).tolist()
                        if value[2] >= 128:    # minus
                            value = value[0] + value[1] * 2**8 + (value[2] - 256) * 2**16
                        else:            # plus
                            value = value[0] + value[1] * 2**8 + value[2] * 2**16
                    # EDF, GDF
                    elif self.HDR.CHANNEL[ch].GDFTYP < 18:
                        bend = bstart + __GDFTYP_BYTE[self.HDR.CHANNEL[ch].GDFTYP]
                        value = scipy.fromstring(raw[bstart:bend], __GDFTYP_NAME[self.HDR.CHANNEL[ch].GDFTYP]).tolist()[0]
                        
                    else:
                        raise FATALERROR('Error SREAD: datatype ' + self.HDR.CHANNEL[ch].GDFTYP + ' not supported!')
                        
                    # calculating new input string position
                    bstart = bend
                    
                    # overflow and saturation detection
                    if (self.HDR.FLAG.OVERFLOWDETECTION != 0) and ((value <= self.HDR.CHANNEL[ch].DigMin) or (value >= self.HDR.CHANNEL[ch].DigMax)):
                        value = NaN
                        
                    # calibration
                    elif self.HDR.FLAG.UCAL == 0:
                        value = value * self.HDR.CHANNEL[ch].Cal + self.HDR.CHANNEL[ch].Off
    
                    # inserting in output matrix
                    self.HDR.data.block[row, column] = value
                
                
                # inserting same value multiply times (self.HDR.CHANNEL[].SPR != self.HDR.SPR)
                else:
                    self.HDR.data.block[row, column] = value
            
            
            self.HDR.data.size = self.HDR.data.block.shape
            self.HDR.FILE.POS = self.HDR.FILE.POS + length
        
            if eegpy.debug:
                print 'Finished reading data'
        return self.HDR.data.block
        # End of SREAD
    
        
        ###########
        # SWRITE  #
        ###########
    def _swrite(self, ptr, nelem):
        """Writing is not implemented yet!"""
        print 'Writing is not implemented yet!'
        #End of SWRITE
        
        
        ###########
        #  SEOF   #
        ###########
    def _seof(HDR):
        """input: HDR_TYPE HDR
        output: [-1, 0]
        Indicates if end of data is reached.
        Returns 0 after the last block, -1 otherwise."""
        
        if self.HDR.FILE.POS >= self.HDR.NRec:
            return True
        return False
        #End of SEOF
        
        
        ###########
        # SREWIND #
        ###########
    def _srewind(self):
        """input: HDR_TYPE HDR
        output: None
        Sets the data pointer back to the beginning.
        No return value."""
        
        if self.HDR.FILE.OPEN != 0:
            self.HDR.FILE.FID.seek(self.HDR.HeadLen)
            self.HDR.FILE.POS = 0
        #End of SREWIND
        
        
        ###########
        #  SSEEK  #
        ###########
    def _sseek(self, offset = 1, whence = 0):
        """input: HDR_TYPE HDR, int offset, [-1, 0, 1] whence
        output: [-1, 0]
        Sets the position pointer to the desired position.
        offset is measured in blocks
        whence:    -1 -> from beginning of data
             0 -> from actual position
             1 -> from end of data
        If an error occurs, the data pointer is not moved and function returns -1, 0 otherwise."""
        
        if self.HDR.FILE.OPEN != 0:
            if whence < 0:
                pos = self.HDR.HeadLen + offset * self.HDR.AS.bpb
            elif whence == 0:
                pos = self.HDR.HeadLen + (self.HDR.FILE.POS + offset) * self.HDR.AS.bpb
            elif whence > 0:
                pos = self.HDR.HeadLen + (self.HDR.NRec + offset) * self.HDR.AS.bpb
            
            if pos < self.HDR.HeadLen or pos > self.HDR.HeadLen + self.HDR.NRec * self.HDR.AS.bpb:
                return -1
            else:
                self.HDR.FILE.FID.seek(pos)
            
            self.HDR.FILE.POS = (pos - self.HDR.HeadLen) / self.HDR.AS.bpb
            return 0
        return -1
        #End of SSEEK
        
    
        ###########
        #  STELL  #
        ###########
    def _stell(self):
        """input: HDR_TYPE HDR
        output: int
        Returns the actual position of the data pointer in blocks.
        If an error occurs, function returns -1."""
        
        if self.HDR.FILE.OPEN != 0:
            pos = self.HDR.FILE.FID.tell()
            if pos == self.HDR.FILE.POS * self.HDR.AS.bpb + self.HDR.HeadLen:
                return self.HDR.FILE.POS
        return -1
        #End of STELL
        
    
    ######################
    # INTERNAL FUNCTIONS #
    ######################
    
    def __gcd(self, a, b):
        while (b != 0):
            t = b
            b = a % b
            a = t
        return a
    
    def __lcm(self, a, b):
        return (abs(a*b)/self.__gcd(a,b))
    
    def __py_time2t_time(self, t):
        delta = t - datetime.datetime(1970,1,1)
        return (delta.days * 86400 + float(delta.seconds + delta.microseconds * pow(10,-6)))
    
    def __t_time2py_time(self, t):
        return (datetime.datetime(1970,1,1) + datetime.timedelta(t / 86400))
    
    def __gdf_time2py_time(self, t):
        if t[14:16] == '  ':
            t[14:16] = '00'
        #month = int(t[4:6])    
        #if month not in range(1,13):
        #    month=1
        try:
            return (datetime.datetime(int(t[0:4]),month,int(t[6:8]),int(t[8:10]),int(t[10:12]),int(t[12:14]),int(t[14:16])*pow(10,4)))
        except Exception, e:
            return datetime.datetime(1900,1,1)
    
    def __py_time2gdf_time(self, t):
        return (t.strftime('%Y%m%d%H%M%S') + str(t.microsecond/pow(10,4)))
    
    def __gdf2_time2py_time(self, t):
        if t == 0:
            t = 367 * pow(2,32)
        return (datetime.datetime(1,1,1) + datetime.timedelta(t * pow(2,-32) - 367))    # don't know if there's an error in datetime but the 367 days shift gives the right result
    
    def __py_time2gdf2_time(self, t):
        delta = t - datetime.datetime(1,1,1)
        if t == datetime.datetime(1,1,1):
            delta = datetime.timedelta(-367)
        return int((delta.days + float(delta.seconds + delta.microseconds * pow(10,-6)) / 86400 + 367) * pow(2,32))
    
    def __gdf2_time2t_time(self, t):
        return ((t * pow(2,-32) - 719529) * 86400)
    
    def __t_time2gdf2_time(self, t):
        return int((float(t) / 86400 + 719529) * pow(2,32))
    
    def get_shape(self):
        return self._shape

    def get_samplRate(self):
        return self.header["samplRate"]
    
    def get_num_datapoints(self):
        return self.HDR.NRec*self.HDR.SPR
    
    def get_channel_names(self):
        return self._channel_names
    
    def getChannelNames(self):
        return self._channel_names
    
    def set_channel_names(self, cns):
        assert len(cns) == self.shape[1], "List of channelnames must contain exactly %i elements" % self.shape[1]
        assert self._mode in ["r+", "w+"], "Cannot set channel_names: F32 is read-only"
        pass
            
    samplRate = property(get_samplRate)
    Fs = property(get_samplRate)
    num_datapoints = property(get_num_datapoints)
    numDatapoints = property(get_num_datapoints)
    channel_names = property(get_channel_names, set_channel_names)
    channelNames = property(get_channel_names)
    
################
# UNTERKLASSEN #
################

class CHANNEL_TYPE:
    OnOff = 1    #
    Label = ''    # Label of channel 
    Transducer = ''    # transducer e.g. EEG: Ag-AgCl electrodes
    PhysDim = ''    # physical dimension
    PhysDimCode = 0    # code for physical dimension
    PreFilt = ''    # pre-filtering

    LowPass = 0    # lowpass filter
    HighPass = 0    # high pass
    Notch = 0    # notch filter
    XYZ = 0        # electrode position
    Impedance = 0    # in Ohm

    PhysMin = 0    # physical minimum
    PhysMax = 0    # physical maximum
    DigMin = 0    # digital minimum
    DigMax = 0    # digital maximum

    GDFTYP = 0    # data type
    SPR = 0        # samples per record (block)
    bpr = 0        # bytes per record (block)

    Cal = 1        # gain factor 
    Off = 0        # bias 

class DATA_TYPE:
    size = scipy.array([0,0])    # size {rows, columns} of data block    
    block =  scipy.array([])    # data block  
    
class IMPAIRMENT_TYPE:
    Visual = 'Unknown'

# Patient specific information
class PATIENT_TYPE:
    Name = ''    # not recommended because of privacy protection 
    Id = ''        # identification code as used in hospital 
    Weight = 0    # weight in kilograms [kg] 0:unkown, 255: overflow 
    Height = 0    # height in centimeter [cm] 0:unkown, 255: overflow 
    Birthday = 0    # Birthday of Patient
    Age = 0
    Headsize = 0    # circumference, nasion-inion, left-right mastoid in millimeter; 
    Sex = 'Unknown'     
    Handedness = 'Unknown'    
    # Patient classification
    Smoking = 'Unknown'
    AlcoholAbuse = 'Unknown'
    DrugAbuse = 'Unknown'
    Medication = 'Unknown'
    Impairment = IMPAIRMENT_TYPE()

class ID_TYPE:
    Technician = ''     
    Hospital = ''
    Equipment = ''    # identfies this software

# location of recording according to RFC1876
class LOC_TYPE:
    VertPre = 0
    HorizPre = 0
    Size = 0
    Version = 0
    Latitude = 0    # in degrees
    Longitude = 0    # in degrees
    Altitude = 0    # in metres

# position of electrodes; see also HDR.CHANNEL[k].XYZ
class ELEC_TYPE:
    REF = ''    # XYZ position of reference electrode
    GND = ''    # XYZ position of ground electrode

# EVENTTABLE 
class EVENT_TYPE:
    SampleRate = 0    # for converting POS and DUR into seconds 
    N = 0        # number of events
    TYP = ''    # defined at http://cvs.sourceforge.net/viewcvs.py/biosig/biosig/t200/eventcodes.txt?view=markup
    POS = 0        # starting position [in samples]
    DUR = 0        # duration [in samples]
    CHN = 0        # channel number; 0: all channels 

# flags
class FLAG_TYPE:
    OVERFLOWDETECTION = 1    # overflow & saturation detection 0: OFF, !=0 ON
    UCAL = 0    # UnCalibration  0: scaling  !=0: NO scaling - raw data return 

# File specific data
class FILE_TYPE:
    FID = ''    # file handle 
    POS = 0        # current reading/writing position in samples 
    OPEN = 0    # 0: closed, 1:read, 2: write
    LittleEndian = ''    # 

# internal variables (not public) 
class AS_TYPE:
    PID = ''    # patient identification
    RID = ''    # recording identification 
    spb = 0        # total samples per block
    bpb = 0        # total bytes per block
    bi = 0
    Header1 = ''
    rawdata = scipy.array([])    # raw data block 

class HDR_TYPE:
    TYPE = ''    # type of file format
    VERSION = 0    # GDF version number 
    FileName = ''

    HeadLen = 0    # length of header in bytes
    NS = 0        # number of channels
    SPR = 0        # samples per block (when different sampling rates are used, this is the LCM(CHANNEL[..].SPR)
    NRec = -1    # number of records/blocks -1 indicates length is unknown.    
    Dur = ''    # Duration of each block in seconds expressed in the fraction Dur[0]/Dur[1]
    SampleRate = 0    # Sampling rate
    IPaddr = ''    # IP address of recording device (if applicable)    
    T0 = ''        # starttime of recording
    
    data = DATA_TYPE()
    Patient = PATIENT_TYPE()
    ID = ID_TYPE()
    LOC = LOC_TYPE()
    ELEC = ELEC_TYPE()
    EVENT = EVENT_TYPE()
    FLAG = FLAG_TYPE()
    FILE = FILE_TYPE()
    AS = AS_TYPE()
    CHANNEL = []
    