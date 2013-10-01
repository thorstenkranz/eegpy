#! /usr/bin/python
# -*- coding: utf8 -*-
import os
import tempfile

from eegpy.formats.iobase import (MemoryMappedBinaryDataFile,
                                  EEGfiltered)

class DBPA(MemoryMappedBinaryDataFile):
    """Access to Sensorium DBPA datafiles"""
    _eventChannels = [0]
    _responseChannels = [1] 
    _dataChannels =  range(2,119)
    _eyeChannels = range(119,122)

    def __init__(self,filename,mode="r+",cNames=None,shape=None, Fs = 1000.0, numChannels=122):
        """Documentation needed"""
        fileSize = os.stat(filename).st_size
        numDatapoints = int(fileSize/numChannels/4)
        MemoryMappedBinaryDataFile.__init__(self, filename, mode, 
                    shape = (numDatapoints, numChannels),
                    data_offset=0,
                    dtype=">f", Fs = Fs, channel_names=cNames)

class DBPAfiltered(EEGfiltered, DBPA):
    """F32 read-only object with included frequency-filteres"""
    def __init__(self,filename,filter_function, **kw_args):
        DBPA.__init__(self,filename,"r+", **kw_args)
        EEGfiltered.__init__(self, filter_function)
        
class DBPAfiltered(DBPA):
    """DBPA read-only object with included frequency-filters"""
    def __init__(self,filename,filter_function):
        DBPA.__init__(self,filename,"r+")
        self._filter_function = filter_function

    @property
    def ff(self):
        return self._filter_function

    def __getitem__(self,item):
        """Calls method of super-class, filters the return value"""
        return self.ff(DBPA.__getitem__(self,item))

class DBPAFilteredWithCache(DBPA):
    """DBPA read-only object with included frequency-filteres"""
    def __init__(self,filename,btype='lp',fl=None,fh=None,border=2,filter_windowed=False,dirname=None):
        from eegpy.ui.eegpylab import freqfilt_eeg
        fd, tmpfn = tempfile.mkstemp(dir=dirname)
        freqfilt_eeg(filename,tmpfn,btype,fl,fh,border,filter_windowed)
        DBPA.__init__(self,tmpfn)

    def __del__(self):
        self.close()
        try:
            os.unlink(self._fn)
        except Exception,e:
            print "Cannot remove file %s"%self._fn, e

if __name__ == "__main__":
    fn = "/home/thorsten/eyecalibrate.000.bh.dat"
    dbpa = DBPA(fn)
