# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR, debug
from eegpy.helper import demean

try:
    import numpy as n
    from scipy.signal import hilbert, detrend
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')


def phase_coherence(x,y):
    """Diese Funktion berechnet aus den Arrays die mittlere Phasenkohärenz R, 
    ein Mass zwischen 0 und 1 (zirkuläre Statistik).
    If the arrays are 1d, calculate the mean over time.
    If they are 2d (second dimensions usually trials), calculate mean over this dimension."""
    assert x.shape == y.shape
    if len(x.shape)==1:
        X = hilbert(demean(x))
        Y = hilbert(demean(y))
        
        try:
            X = X/abs(X)
            Y = Y/abs(Y)
        except ArithmeticException, ae:
            print "Error: Division by zero", e
        except:
            print "Anderer Fehler! (Vielleicht Index?)"
        
        pd = n.exp(n.log(X)-n.log(Y))
        return abs(pd.mean())
    elif len(x.shape)==2:
        X = n.zeros((x.shape),"D")
        Y = n.zeros((y.shape),"D")
        for i in range(x.shape[1]):
            X[:,i] = hilbert(demean(x[:,i]))
            Y[:,i] = hilbert(demean(y[:,i]))
        try:
            for i in range(X.shape[1]):
                X[:,i] /= abs(X[:,i])
                Y[:,i] /= abs(Y[:,i])
        except ArithmeticException, ae:
            print "Error: Division by zero", e
        
        pd = n.exp(n.log(X)-n.log(Y))
        return abs(pd.mean(axis=1))
        

#def phase_coherence_fetptrials

def phases(data,continuous=False):
    """Berechnung der Phasen einer Funktion mit Hilfe der Hilbert-Transformation"""
    transf = hilbert(data)
    #phasen = zeros(len(data),"d")
    phasen = n.arctan2(transf.imag,transf.real)
    #for i in range(len(data)):
    #    phasen[i] = n.arctan(#atan(transf.imag[i]/transf.real[i])
    #    if transf.real[i] < 0:# and transf.imag[i] > 0:
    #        phasen[i] = phasen[i]+pi
    #    phasen[i] = phasen[i]+pi/2
        #if transf.imag[i] < 0 and transf.real[i] > 0:
        #    phases[i] = phases[i]-pi/2
    if continuous==False:
        return phasen
    else:
        return make_phases_continuous(phasen)

def make_phases_continuous(phases):
    """Makes an 1d-array of real valued phases continuous"""
    if len(phases.shape) == 1:
        diff = n.diff(phases)
        PSP = n.arange(phases.shape[0])[diff>n.pi]+1
        NSP = n.arange(phases.shape[0])[diff<-n.pi]+1
        for k in PSP:
            phases[k:] -= 2*n.pi
        for k in NSP:
            phases[k:] += 2*n.pi
        return phases
    elif len(phases.shape) == 2:
        for i in range(phases.shape[1]):
            diff = n.diff(phases[:,i])
            PSP = n.arange(phases.shape[0])[diff>n.pi]+1
            NSP = n.arange(phases.shape[0])[diff<-n.pi]+1
            for k in PSP:
                phases[k:,i] -= 2*n.pi
            for k in NSP:
                phases[k:,i] += 2*n.pi
        return phases
    else:
        raise ValueError("Only 1d and 2d arrays supported")
    
def phases_from_complex(wts, continuous=False, do_detrend=False):
    """Calculates phases from 1d or 2d wavelet/hilbert arrays, dim0 is time"""
    if len(wts.shape) == 1:
        #1d
        phasen = n.arctan2(wts.imag,wts.real)
        if not (continuous or do_detrend):
            return phasen
        else:
            phasen = make_phases_continuous(phasen)
            if do_detrend:
                phasen = detrend(phasen,axis=0)
            return phasen
    elif len(wts.shape) == 2:
        #2d
        phasen = n.arctan2(wts.imag,wts.real)
        if not (continuous or do_detrend):
            return phasen
        else:
            phasen = make_phases_continuous(phasen)
            if do_detrend:
                phasen = detrend(phasen,axis=0)
            return phasen
        
    else:
        raise ValueError("Only 1d and 2d arrays supported")


if __name__ == "__main__":
    import pylab as p
    from eegpy.analysis.wavelet import wt_analyze
    eeg = eegpy.F32("/media/story/SchlafGed/iEEG/data/canseven_bp.f32")
    #p.plot(eeg[10000:20000,:])
    data = eeg[21000:30000,20]
    freqs=n.arange(1,15,5)
    wts = wt_analyze(data,freqs=freqs)
    for i in range(wts.shape[1]):
        p.subplot(511)
        p.plot(phases_from_complex(wts[:,i],continuous=False))
        p.subplot(512)
        phase_cont = phases_from_complex(wts[:,i],continuous=True)
        p.plot(phase_cont)
        xs = n.arange(phase_cont.shape[0])
        pfacs = p.polyfit(xs[0:500],phase_cont[0:500],1)
        print pfacs
        p.plot(xs,p.polyval(pfacs,xs),"k--",)
        p.subplot(513)
        p.plot(phase_cont-p.polyval(pfacs,xs))
        p.subplot(514)
        pfacs = ((freqs[i]*2*np.pi)/1000.0 , 0)
        print pfacs
        tmp = phase_cont-p.polyval(pfacs,xs)
        tmp -= tmp[0:500].mean()
        p.plot(tmp)
        p.subplot(515)
        p.plot(phases_from_complex(wts[:,i],do_detrend=True))

