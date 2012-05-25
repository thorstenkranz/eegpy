#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""General-purpose utilities and helpers"""

from eegpy.misc import FATALERROR, EegpyBase

import pickle
import sys
import tempfile
import re


try:
    import numpy as np
    from numpy.fft import fft, ifft
    from scipy.interpolate import splrep,splev,interp1d
    from scipy.signal import resample, detrend
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

    Uses the pickle-module.

    :param fn: The filename
    :type fn: string

    :returns: Unpickled object, usually subclass of type EegpyBase
    """
    try:
        rv = pickle.load(open(fn,"rb"))
        return rv
    except Exception, e:
        print "Cannot load object from disk, ", e 

def fconv(x,h):
    """Implements the fast convolution
    
    :param x: Array one
    :type x: array
    :param h: Array two
    :type h: array

    :returns: convolved array    
    """
    Ly=len(x)+len(h)-1
    Ly2=int(2**(nextpow2(Ly)))       # Find smallest power of 2 that is > Ly
    X=fft(x, Ly2)                # Fast Fourier transform
    H=fft(h, Ly2)               # Fast Fourier transform
    Y=X*H                      # 
    y=(ifft(Y, Ly2))#.real      # Inverse fast Fourier transform
    y=y[0:Ly]               # Take just the first N elements
    if x.dtype==np.double:
        return y.real
    else:
        return y
    #y=y/max(abs(y));           # Normalize the output
    
def nextpow2(x):
    """
    Determine next power of 2
    
    Calculate smallest possible n so that
    
    .. math::

        n \in N,\; 2^n>x

    :param x: Number for which the next power of 2 should be found.
    :type x: float or int
    """
    if x<0:
        raise ValueError("x must be non-negative")
    return np.ceil(np.log(x)/np.log(2))

def is_power_of_2(x):
    """Check if x is a power of two. 

    Cp.: http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
    """
    try:
        int(x)
    except:
        raise ValueError("Cannot convert to int")
    if not int(x) == x:
        return False
    x = int(x)
    return bool(x) and not bool(x & (x - 1))

def factorial(x):
    """Calculate the factorial of x"""
    x = int(x)
    return np.arange(1,x+1).prod()

def findFirstMax(x,thres = 200,search_width=500):
    #Findet in einem 1d oder 2d-Array (wie getData es liefert) die Stelle des ersten Imaging-Artefakt-Extremums
    assert len(x.shape) in [1,2]
    if len(x.shape)==1:
        diff = np.diff(abs(x))
    else:
        diff = np.diff(abs(x).mean(axis=1))
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

def demean(x,axis=0,bp=0): 
    """
    Demean array along one axis.
    
    :param x: array to demean. Remains unchanged.
    :type x: array
    :param axis: the axis along which to demean. Default is 0, not -1 as in scipy!
    :type axis: int
    :param bp: A sequence of break points. If given, an individual fit is
               performed for each part of `x` between two break points.
               Break points are specified as indices into `x`.
    :type bp: array_like of ints, optional
    
    :returns: Demeaned array.
    """
    return detrend(x, axis, "constant", bp)

def taper(sig, tlip = 0.125, tapertype = 'sinus'):
    # taper for wavelet and other purposes
    # tlip = Taper Length in prozent/100 : für jede Seite
    data = np.copy(sig)
    if tapertype == 'sinus':
        nd = int(len(data)*tlip)
        sw = (np.pi/2) / nd
        hc = np.arange(0,(np.pi/2),sw)
        halbwelle = np.sin(hc)
        for i in range(nd):
            data[i] *= halbwelle[i]
            data[-i] *= halbwelle[i]
    return data

def upsample(x,factor=10):
    """
    Performs upsampling using splines. 
    
    Now works for 1d, 2d and 3d arrays.
    
    Uses slices to perform simple downsampling. No temporal filtering is being 
    performed to cope with aliasing. See :py:meth:`downsample` for upsampling.

    :param x: array to upsample. Remains unchanged.
    :type x: array
    :param factor: factor by which the upsampling is being perfomred.
    :type factor: int

    :returns: upsampled array
    """
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for the upsampling-factor")
    if len(x.shape) == 1: #1d 
        #y = np.zeros((factor*x.shape[0]),"d")
        spl = splrep(np.arange(x.shape[0]),x)
        y = splev(np.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 2: #1d 
        y = np.zeros((factor*x.shape[0],x.shape[1]),"d")
        for i in range(x.shape[1]):
            spl = splrep(np.arange(x.shape[0]),x[:,i])
            y[:,i] = splev(np.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 3: #1d 
        y = np.zeros((factor*x.shape[0],x.shape[1],x.shape[2]),"d")
        for i in range(x.shape[1]):
            for j in range(x.shape[2]):
                spl = splrep(np.arange(x.shape[0]),x[:,i,j])
                y[:,i,j] = splev(np.arange(0,x.shape[0],1.0/factor),spl)
    else:
        raise ValueError("Only 1d, 2d and 3d arrays are accepted for upsampling up to now.")
    return y 

def tmp_memmap(*args,**kwargs):
    tmp_fn = tempfile.mktemp()
    y=np.memmap(tmp_fn,*args,**kwargs)
    return y, tmp_fn

def upsample_to_memmap(x,factor=10):
    """Performs upsampling using splines and returns a memmap. Now works for 1d, 2d and 3d arrays."""
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for the upsampling-factor")
    if len(x.shape) == 1: #1d 
        #y = np.zeros((factor*x.shape[0]),"d")
        spl = splrep(np.arange(x.shape[0]),x)
        #f = interp1d(np.arange(x.shape[0]),x)
        y,tmp_fn = tmp_memmap(dtype=np.double,shape=(x.shape[0]*factor),mode="w+")
        y[:] = splev(np.arange(0,x.shape[0],1.0/factor),spl)
        #y[:-factor] = f(np.arange(0,x.shape[0],1.0/factor)[:-factor])
        #y[:] = resample(x,x.shape[0]*factor)
    #TODO: make kind of interpolation consistent for different dimensionalities
    elif len(x.shape) == 2: #1d 
        y,tmp_fn = tmp_memmap(dtype=np.double,shape=(x.shape[0]*factor,x.shape[1]),mode="w+")
        for i in range(x.shape[1]):
            spl = splrep(np.arange(x.shape[0]),x[:,i])
            y[:,i] = splev(np.arange(0,x.shape[0],1.0/factor),spl)
    elif len(x.shape) == 3: #1d 
        y,tmp_fn = tmp_memmap(dtype=np.double,shape=(x.shape[0]*factor,x.shape[1],x.shape[2]),mode="w+")
        for i in range(x.shape[1]):
            for j in range(x.shape[2]):
                spl = splrep(np.arange(x.shape[0]),x[:,i,j])
                y[:,i,j] = splev(np.arange(0,x.shape[0],1.0/factor),spl)
    else:
        raise ValueError("Only 1d, 2d and 3d arrays are accepted for upsampling up to now.")
    return y, tmp_fn
    
def downsample(x,factor=10,start=0):
    """
    Downsampling using slices.
    
    Uses slices to perform simple downsampling. No temporal filtering is being 
    performed to cope with aliasing. See :py:meth:`upsample` for upsampling.

    :param x: array to downsample. Remains unchanged.
    :type x: array
    :param factor: factor by which the downsampling is being perfomred.
    :type factor: int
    :param start: at which timepoints start slicing
    :type start: int

    :returns: downsampled array
    """
    if start<-x.shape[0] or start>x.shape[0]-1:
        raise ValueError("Illegal value for start")
    factor = int(round(factor))
    if factor<1 or factor>1000000:
        raise ValueError("Illegal value for factor")
    y_shape = list(x.shape)
    y_shape[0] = x.shape[0]/factor
    y = np.zeros(y_shape,"d")
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
        x[startdel:enddel] = np.zeros(x[startdel:enddel].shape,"d")
    maxs.sort()
    return maxs

def find_all_maxs(x,ratio=0.5,width=86,nth_global=1,callback=None):
    """In a 1d-array x, find the indices of all local maxima that are at least "ratio" the siz of the global maximum.
    Assume width to be the minimal distance between the maxima.
    New: parameter nth_global. If >1, do cutoff with respect to the nth global maximum. 
    Used to coke with artifactual maxima. """
    x = x.copy()
    #Remove borders as they might disturb the result
    x[:width] = np.zeros((width),"d")
    x[-width:] = np.zeros((width),"d")
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
        x[max-width/2:max+width/2] = np.zeros(x[max-width/2:max+width/2].shape,"d")
        #print max, len(maxs), last_ratio
        #sys.stdout.flush()
    maxs.sort()
    return maxs



def os_display(msg,timeout=100,wait_for_enter=True):
    """Display a simple onscreen message, e.g., to signal completion of a long running task
    
    :param msg: The message to be diplayed
    :type msg: string
    :param timeout: How long to show the message. If wait_for_enter, this setting is ignored.
    :type timeout: int
    :param wait_for_enter: Show message until <Enter> is has been pressed. 
    :type wait_for_enter: bool
    """
    if pyosd!=None:
        opt = {"font":#"-*-fixed-*-*-*-*-150-*-*-*-c-90-*-*",
                "-adobe-helvetica-bold-r-normal-*-*-320-*-*-p-*-*",
               "colour": "#FF0000",
               "shadow": 2,
               "offset": 300,
               "pos": pyosd.POS_TOP,
               "timeout": timeout}
        osd = pyosd.osd(font=opt["font"] ,colour=opt["colour"],timeout=opt["timeout"],shadow=opt["shadow"],offset=opt["offset"], pos=opt["pos"])
        #osd = pyosd.osd()
        osd.set_horizontal_offset(10)
        osd.display(msg)
        if wait_for_enter:
            raw_input("Press <enter> to go on")
        else:
            osd.wait_until_no_display()
    else:
        print "pyosd is not available, but needed for displaying on-screen messages."
        
class ProgressBar(EegpyBase): 
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
    """
    Make copy of an array and prepend AND append zeros.
    
    :param x: array to be extended with zeros. Stays unchanged. Only 1d- and 2d-arrays supported.
    :type x: array
    :param num_z: number of zeros appended on each side. 
    :type num_z: int

    .. hint::

        .. plot:: 

            import matplotlib.pyplot as p
            import numpy as np
            from eegpy.helper import prepend_append_zeros

            array1 = np.random.random((10,10))
            array2 = prepend_append_zeros(array1,3)

            p.subplot(121)
            p.imshow(array1,interpolation="nearest")
            p.title("(10,10)")
            p.subplot(122)
            p.imshow(array2,interpolation="nearest")
            p.title("num_z = 3")
            p.colorbar()
    """
    if not len(x.shape) in [1,2]:
        raise ValueError("Must be 2d")
    if len(x.shape) == 1:
        rv = np.zeros((x.shape[0]+2*num_z))
        rv[num_z:num_z+x.shape[0]] = x[:]    
    else:
        rv = np.zeros((x.shape[0]+2*num_z,x.shape[1]))
        rv[num_z:num_z+x.shape[0],:] = x[:]
    return rv

def fprobi(df1,df2,prob):
    """Calculate the critical F-value for given degrees of freedom and p-value.
    Uses optimization. Not optimal, but works"""
    def err(p,fp_ref):
        #print p,fp_ref
        return abs(fp_ref-fprob(df1,df2,p[0]))
    p = [1.68]
    return fmin(err,p,np.array([prob]),disp=False)[0]

def ksprobi(n1,n2,prob):
    """Calculate the critical k-value for given p-value.
    Uses optimization. Not optimal, but works"""
    def err(p,fp_ref):
        #print p,fp_ref
        #print "p, fp_ref:", p, fp_ref
        return abs(fp_ref-ksprob((en+0.12+0.11/en)*p[0]))
        #return abs(fp_ref-ksprob(ksprob((en+0.12+0.11/en)*p)))

    en = np.sqrt(n1*n2/float(n1+n2))
    p = [0.3]
    return fmin(err,p,np.array([prob]),maxiter=1000,disp=False)[0]
    #return brute(err,[np.s_[0:0.5:0.5/100]],args=[prob])[0]
        
def cm2in(length):
    """Convert cm to inches, used for pylab"""
    return length*.393700787

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.

        Example:

        .. sourcecode:: ipython

            In [1]: from eegpy.helper import sort_nicely
            In [2]: l = ["image9.jpg","image10.jpg"]
            In [3]: l.sort()
            In [4]: l
            Out[4]: ['image10.jpg', 'image9.jpg']
            In [5]: sort_nicely(l)
            In [6]: l
            Out[6]: ['image9.jpg', 'image10.jpg']
        
    """
    def tryint(s):
        try:
            return int(s)
        except:
            return s
        
    def alphanum_key(s):
        """ Turn a string into a list of string and number chunks.
            "z23a" -> ["z", 23, "a"]
        """
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]

    l.sort(key=alphanum_key)
