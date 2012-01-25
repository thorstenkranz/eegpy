#!/usr/bin/env python
# -*- coding: utf-8 -*-

from eegpy.misc import FATALERROR

import pickle
import sys
import tempfile

try:
    import pylab as P
except Exception,e:
    print "helper.py: Cannot import pylab. I continue, let's see what happens"

try:
    import numpy as n
    from numpy.fft import *
    from scipy.interpolate import splrep,splev,interp1d
    from scipy.signal import resample
    from scipy.optimize import fmin, fmin_powell, brute
    from scipy.stats import fprob, ksprob
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pyosd
except ImportError:
    pyosd=None

def load(fn):
    """Load a saved object from disk.
    Uses the pickle-module."""
    try:
        rv = pickle.load(open(fn,"rb"))
        return rv
    except Exception, e:
        print "Cannot load object from disk, ", e 

def fconv(x,h):
    """Implements the fast convolution"""
    Ly=len(x)+len(h)-1
    Ly2=int(2**(nextpow2(Ly)))       # Find smallest power of 2 that is > Ly
    X=fft(x, Ly2)                # Fast Fourier transform
    H=fft(h, Ly2)               # Fast Fourier transform
    Y=X*H                      # 
    y=(ifft(Y, Ly2))#.real      # Inverse fast Fourier transform
    y=y[0:Ly]               # Take just the first N elements
    return y
    #y=y/max(abs(y));           # Normalize the output
    
def nextpow2(x):
    return n.ceil(n.log(x)/n.log(2))

def is_power_of_2(x):
    """Check if x i a power of two. 
    Cp.: http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
    """
    if not int(x) == x:
        return False
    x = int(x)
    return bool(x) and not bool(x & (x - 1))

def factorial(x):
    x = int(x)
    return n.arange(1,x+1).prod()

def findFirstMax(x,thres = 200,search_width=500):
    #Findet in einem 1d oder 2d-Array (wie getData es liefert) die Stelle des ersten Imaging-Artefakt-Extremums
    assert len(x.shape) in [1,2]
    if len(x.shape)==1:
        diff = n.diff(abs(x))
    else:
        diff = n.diff(abs(x).mean(axis=1))
    import pylab
    #pylab.ion()
    #pylab.plot(diff)
    #pylab.plot(x)
    #pylab.show()
    #raw_input("Test")
    for i in range(search_width):
        if diff[i] < -thres:
            return i
    return -1

def demean(sig):  
    #Returns an array X with X.mean()==0
    X = n.copy(sig)
    X -= X.mean()
    return X

def taper(sig, tlip = 0.125, tapertype = 'sinus'):
    # taper for wavelet and other purposes
    # tlip = Taper Length in prozent/100 : für jede Seite
    data = n.copy(sig)
    if tapertype == 'sinus':
        nd = int(len(data)*tlip)
        sw = (n.pi/2) / nd
        hc = n.arange(0,(n.pi/2),sw)
        halbwelle = n.sin(hc)
        for i in range(nd):
            data[i] *= halbwelle[i]
            data[-i] *= halbwelle[i]
    return data

def upsample(x,factor=10):
    """Performs upsampling using splines. Now works for 1d, 2d and 3d arrays."""
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for the upsampling-factor")
    if len(x.shape) == 1: #1d 
        #y = n.zeros((factor*x.shape[0]),"d")
        spl = splrep(n.arange(x.shape[0]),x)
        y = splev(n.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 2: #1d 
        y = n.zeros((factor*x.shape[0],x.shape[1]),"d")
        for i in range(x.shape[1]):
            spl = splrep(n.arange(x.shape[0]),x[:,i])
            y[:,i] = splev(n.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 3: #1d 
        y = n.zeros((factor*x.shape[0],x.shape[1],x.shape[2]),"d")
        for i in range(x.shape[1]):
            for j in range(x.shape[2]):
                spl = splrep(n.arange(x.shape[0]),x[:,i,j])
                y[:,i,j] = splev(n.arange(0,x.shape[0],1.0/factor),spl)
    else:
        raise ValueError("Only 1d, 2d and 3d arrays are accepted for upsampling up to now.")
    return y 

def tmp_memmap(*args,**kwargs):
    tmp_fn = tempfile.mktemp()
    y=n.memmap(tmp_fn,*args,**kwargs)
    return y, tmp_fn

def upsample_to_memmap(x,factor=10):
    """Performs upsampling using splines and returns a memmap. Now works for 1d, 2d and 3d arrays."""
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for the upsampling-factor")
    if len(x.shape) == 1: #1d 
        #y = n.zeros((factor*x.shape[0]),"d")
        spl = splrep(n.arange(x.shape[0]),x)
        #f = interp1d(n.arange(x.shape[0]),x)
        y,tmp_fn = tmp_memmap(dtype=n.double,shape=(x.shape[0]*factor),mode="w+")
        y[:] = splev(n.arange(0,x.shape[0],1.0/factor),spl)
        #y[:-factor] = f(n.arange(0,x.shape[0],1.0/factor)[:-factor])
        #y[:] = resample(x,x.shape[0]*factor)
    #TODO: make kind of interpolation consistent for different dimensionalities
    elif len(x.shape) == 2: #1d 
        y,tmp_fn = tmp_memmap(dtype=n.double,shape=(x.shape[0]*factor,x.shape[1]),mode="w+")
        for i in range(x.shape[1]):
            spl = splrep(n.arange(x.shape[0]),x[:,i])
            y[:,i] = splev(n.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 3: #1d 
        y,tmp_fn = tmp_memmap(dtype=n.double,shape=(x.shape[0]*factor,x.shape[1],x.shape[2]),mode="w+")
        for i in range(x.shape[1]):
            for j in range(x.shape[2]):
                spl = splrep(n.arange(x.shape[0]),x[:,i,j])
                y[:,i,j] = splev(n.arange(0,x.shape[0],1.0/factor),spl)
    else:
        raise ValueError("Only 1d, 2d and 3d arrays are accepted for upsampling up to now.")
    return y, tmp_fn
    
def downsample(x,factor=10,start=0):
    "Downsampling using slices"
    if start<-x.shape[0] or start>x.shape[0]-1:
        raise ValueError("Illegal value for the start value")
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for the upsampling-factor")
    y_shape = list(x.shape)
    y_shape[0] = x.shape[0]/factor
    y = n.zeros(y_shape,"d")
    offset = 0
    while start<0:
        offset+=1
        start+=factor
    tmp = x[:x.shape[0]-offset*factor:factor,...]
    y[offset:offset+tmp.shape[0],...] = tmp[:]
    return y

def find_max_overlap(x,y,max_shift=10):
    """Finds the value for a shift of y for which the mean squared error becomes minimal."""
    assert max_shift<len(x), "max_shift too large."
    assert len(x) == len(y), "works only for equal-size arrays"
    _len = len(x)
    rv = 0
    rv_error = 1000000000
    for i in range(-max_shift,max_shift+1):
        if i<0:
            err = ((x[-i:]-y[:i])**2).sum()/(_len-abs(i))
            #print "Fall1", i,
        elif i==0:
            err = ((x-y)**2).sum()/_len
            #print "Fall2", i,
        else:
            err = ((x[:-i]-y[i:])**2).sum()/(_len-i)
            #print "Fall3", i,
        #print err
        if err<rv_error:
            rv = i
            rv_error=err
            #print "Fehler: ", err, i
        #print "Fehler: ", err
        #if abs(i)<=3:
            #print i, err, "   ",
    #print rv, rv_error
    return rv

def find_maxs(x,num_maxs=35,width=86):
    """In a 1d-array x, find the indices of the num_maxs tallest local maxima.
    Assume width to be the minimal distance between the maxima."""
    x = x.copy()
    maxs = []
    for i in range(num_maxs):
        #P.plot(x)
        #raw_input()
        max = x.argmax()
        maxs.append(max)
        #print "Max", max
        startdel=0
        enddel=-1
        if max-width/2>0:
            startdel=max-width/2
        if max+width/2<x.shape[0]:
            enddel=max+width/2
        #print "sd,ed", startdel,enddel
        x[startdel:enddel] = n.zeros(x[startdel:enddel].shape,"d")
    maxs.sort()
    return maxs

def find_all_maxs(x,ratio=0.5,width=86,nth_global=1,callback=None):
    """In a 1d-array x, find the indices of all local maxima that are at least "ratio" the siz of the global maximum.
    Assume width to be the minimal distance between the maxima.
    New: parameter nth_global. If >1, do cutoff with respect to the nth global maximum. 
    Used to coke with artifactual maxima. """
    x = x.copy()
    #Remove borders as they might disturb the result
    x[:width] = n.zeros((width),"d")
    x[-width:] = n.zeros((width),"d")
    maxs = []
    last_ratio=1
    n_global = nth_global
    global_max=10e50
    while last_ratio>ratio or n_global>0:
        #print nth_global, n_global, global_max, ratio, last_ratio
        if len(maxs)%1000==999:
            #print len(maxs)+1
            if not callback is None:
                callback(len(maxs))
        if n_global > 0:
            n_global -= 1
            global_max=x.max()
        max = x.argmax()
        maxs.append(max)
        last_ratio = x[max]/global_max
        startdel=0
        enddel=-1
        if max-width/2>0:
            startdel=max-width/2
        if max+width/2<x.shape[0]:
            enddel=max+width/2
        x[max-width/2:max+width/2] = n.zeros(x[max-width/2:max+width/2].shape,"d")
        #print max, len(maxs), last_ratio
        #sys.stdout.flush()
    maxs.sort()
    return maxs



def os_display(msg,timeout=100,wait_for_enter=True):
    if pyosd!=None:
        opt = {"font":"-microsoft-verdana-bold-*-*-*-180-*-*-*-*-*-*-*",
               "colour": "#FF0000",
               "shadow": 2,
               "offset": 300,
               "pos": pyosd.POS_TOP,
               "timeout": timeout}
        osd = pyosd.osd(font=opt["font"] ,colour=opt["colour"],timeout=opt["timeout"],shadow=opt["shadow"],offset=opt["offset"], pos=opt["pos"])
        osd.set_horizontal_offset(10)
        osd.display(msg)
        if wait_for_enter:
            raw_input("Press <enter> to go on")
        else:
            osd.wait_until_no_display()
    else:
        print "pyosd is not available, but needed for displaying on-screen messages."
        
class ProgressBar: 
    def __init__(self, finalcount, progresschar=None):
        import sys
        self.finalcount=finalcount
        self.blockcount=0
        #
        # See if caller passed me a character to use on the
        # progress bar (like "*").  If not use the block
        # character that makes it look like a real progress
        # bar.
        #
        if not progresschar: self.block="="#chr(178)
        else:                self.block=progresschar
        #
        # Get pointer to sys.stdout so I can use the write/flush
        # methods to display the progress bar.
        #
        self.f=sys.stdout
        #
        # If the final count is zero, don't start the progress gauge
        #
        if not self.finalcount : return
        #self.f.write('\n0----------------- % Progress -------------------1\n')
        #self.f.write('    1    2    3    4    5    6    7    8    9    0\n')
        #self.f.write('----0----0----0----0----0----0----0----0----0----0\n')
        return

    def progress(self, count):
        #
        # Make sure I don't try to go off the end (e.g. >100%)
        #
        count=min(count, self.finalcount)
        #
        # If finalcount is zero, I'm done
        #
        if self.finalcount:
            percentcomplete=int(round(100*count/self.finalcount))
            if percentcomplete < 1: percentcomplete=1
        else:
            percentcomplete=100
            
        #print "percentcomplete=",percentcomplete
        blockcount=int(percentcomplete/2)
        
        #print "blockcount=",blockcount
        if blockcount > self.blockcount:
            line="0"+self.block*blockcount+(50-blockcount)*"_"+"1"
            #for i in range(self.blockcount,blockcount):
            self.f.write("\r")
            self.f.write(line)
            self.f.flush()
                
        if percentcomplete == 100: self.f.write("\n")
        self.blockcount=blockcount
        return
    
#def simple_progress(fraction,notify_step=0.01):

def prepend_append_zeros(x,num_z=30):
    if not len(x.shape) in [1,2]:
        raise ValueError("Must be 2d")
    if len(x.shape) == 1:
        rv = n.zeros((x.shape[0]+2*num_z))
        rv[num_z:num_z+x.shape[0]] = x[:]    
    else:
        rv = n.zeros((x.shape[0]+2*num_z,x.shape[1]))
        rv[num_z:num_z+x.shape[0],:] = x[:]
    return rv

def fprobi(df1,df2,prob):
    """Calculate the critical F-value for given degrees of freedom and p-value.
    Uses optimization. Not optimal, but works"""
    def err(p,fp_ref):
        #print p,fp_ref
        return abs(fp_ref-fprob(df1,df2,p[0]))
    p = [1.68]
    return fmin(err,p,n.array([prob]),disp=False)[0]

def ksprobi(n1,n2,prob):
    """Calculate the critical k-value for given p-value.
    Uses optimization. Not optimal, but works"""
    def err(p,fp_ref):
        #print p,fp_ref
        #print "p, fp_ref:", p, fp_ref
        return abs(fp_ref-ksprob((en+0.12+0.11/en)*p[0]))
        #return abs(fp_ref-ksprob(ksprob((en+0.12+0.11/en)*p)))

    en = n.sqrt(n1*n2/float(n1+n2))
    p = [0.3]
    return fmin(err,p,n.array([prob]),maxiter=1000,disp=False)[0]
    #return brute(err,[n.s_[0:0.5:0.5/100]],args=[prob])[0]
        
def cm2in(length):
    """Convert cm to inches, used for pylab"""
    return length*.393700787
