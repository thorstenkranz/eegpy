#!/usr/bin/env python
# -*- coding: utf-8 -*-

from eegpy.misc import FATALERROR
import sys
import pylab as p

try:
    import numpy as np
    from numpy.fft import *
    from scipy.integrate import odeint
    from scipy.optimize import anneal, fmin, fmin_powell, brute
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pyosd
except ImportError:
    pyosd=None

def expand_array(ar, mode="per"):
    if mode == "per":
        bigar = np.zeros(np.array(ar.shape)+2)
        # Mitte
        bigar[1:-1,1:-1] = ar
        # Ecken
        bigar[0,0] = ar[-1,-1]
        bigar[-1,0] = ar[0,-1]
        bigar[0,-1] = ar[-1,0]
        bigar[-1,-1] = ar[0,0]
        #Kanten
        bigar[0,1:-1] = ar[-1,:]
        bigar[1,1:-1] = ar[0,:]
        bigar[1:-1,0] = ar[:,-1]
        bigar[1:-1,1] = ar[:,0]
        return bigar
    else:
        raise ValueError("expand_array: mode must one of ['per']")

class CNN2d(object):
    """My CNN implementation. 
    Only 2d arrays, only 3x3 templates, only periodic border, only sigmoidal output, ...
    
    Names:
        input: u
        state: x
        output: y
        Feedforward-template: a
        Feedback-template: b
        bias:z
    """

    def __init__(self,input=None,state=None,ff_template=None,fb_template=None,bias=0):
        """Do initialisation"""
        #Input and state
        if input==None:
            self.u = None
        else:
            self.u = np.array(input)
        if state==None:
            if self.u == None:
                raise ValueError("At least one of input or state must be set to determine shape")
            else:
                self.x = np.zeros(self.u.shape)
        else:
            self.x = np.array(state)
        if self.u == None and self.x != None:
            self.u = np.zeros(self.x.shape)
        self.bigu = expand_array(self.u)

        # Templates
        self.a = np.array(ff_template)
        assert self.a.shape == (3,3)
        self.b = np.array(fb_template)
        assert self.b.shape == (3,3)
        # Bias
        self.z = float(bias)
        #TODO: Do some sanity checks

    def sigmoid(self,x):
        return 2.0/(1.0+np.exp(-x)) - 1

    @property
    def output(self):
        return self.sigmoid(self.x)

    def ode(self, x, t):
        """The ode describing the CNN"""
        x = x.reshape((self.x.shape))
        #self.x = x
        #print t, self.output.mean()
        bigx = expand_array(x)
        #print x.shape, bigx.shape
        bigy = self.sigmoid(bigx)
        dydt = (-1) * x + self.z
        #print dydt.mean()
        for os1 in range(self.a.shape[0]): # os1 = offset1
            for os2 in range(self.a.shape[1]): # os2 = offset2
                dydt += self.a[os1,os2] * self.bigu[os1:os1+dydt.shape[0],os2:os2+dydt.shape[1]]
                dydt += self.b[os1,os2] * bigy[os1:os1+dydt.shape[0],os2:os2+dydt.shape[1]]
        #print dydt.shape
        return dydt.flatten()

    def integrate(self, ts = [50,100,150,200]):
        """Integrate the CNN and return values at times ts"""
        vals = odeint(self.ode,self.x.flatten(),ts,)
        self.x[:,:] = vals[-1,:].reshape(self.x.shape)
        return vals

class CNNTrainer(object):
    """Training CNNs"""
    def __init__(self, a0=None,b0=None,z0=None,errfunc="mean_abs_diff",opt_method="anneal",t_sim=50):
        """Takes lists of arrays and performs training """
        for var in ["a0","b0","z0","errfunc","opt_method","t_sim"]:
            exec("self.%s = %s"%(var,var))
        self.do_pp = False

    def __call__(self, inputs,states,refs):
        if self.a0 == None:
            a0 = np.random.random((3,3))*2-1
        if self.b0 == None:
            b0 = np.random.random((3,3))*2-1
        if self.z0 == None:
            z0 = np.random.random()*2-1

        params = np.array(list(a0.flatten())+list(b0.flatten())+[z0])
        print "Initial parameters:", params
        assert params.shape[0] == 19
        
        #TODO: assert all arrays have same dimensions

        # Do optimization
        if self.opt_method == "anneal":
            #Do the optimization 100 times
            par, fopt, it, fcalls, t5, t6, t7 = anneal(self.error,params, lower=-2, upper=2, full_output=True)
        elif self.opt_method == "simplex":
            #Do the optimization 100 times
            #best_fopt = 10e50
            #for i in range(2):
            #    print "\nStep", i
            #    params = np.random.random((19))*0.2-0.1
            #    print "Initial guess", params
            #    self.do_pp = True
                #par_tmp, ftol, it, fcalls,wf, allv_tmp = fmin(self.error,params,xtol=-1,ftol=-1,maxiter=2000,maxfun=10e10,full_output=True,retall=True)
            par, fopt, it, fcalls,wf, allv_tmp = fmin(self.error,params,maxiter=2000,maxfun=10e10,full_output=True,retall=True)
            #    print fopt,# par_tmp
            #    print "Final:", par_tmp
            #    print "Veränderungen:", par_tmp-params
            #    if fopt < best_fopt:
            #        best_fopt=fopt
            #        par=par_tmp
            #        allv = allv_tmp
            #        print "besser!",
            #par, ftol, it, fcalls,wf, allv = fmin(self.error,params,xtol=-1,ftol=-1,maxiter=2000,maxfun=10e10,full_output=True,retall=True)
            #print "params:", params
            #par, ftol, it, fcalls,wf, allv = fmin(self.error,params,full_output=True,retall=True)
            #print par, ftol, it, fcalls, wf, allv
            #print np.array(allv).shape
                    #allv = np.array(allv)
                    #p.figure(3)
                    #for i in range(allv.shape[1]):
                    #    p.plot(allv[:,i])
                    #p.show()
                    #time.sleep(3)
        elif self.opt_method == "brute":
            par = brute(self.error,[slice(-2,2,2j) for i in range(19)])
            print par
            #par = par[0]

        elif self.opt_method == "powell":
            par = fmin_powell(self.error,params)
            print "Veränderungen:", par-params
        else:
            raise ValueError("Optimization method not known")

        par = np.array(par)
        return par[:9].reshape((3,3)), par[9:18].reshape((3,3)), par[18]
        
    def error(self, p):
        if self.do_pp:
            #print "in error:",  p
            self.do_pp=False
        err = 0
        #print "parameters[0:5]:", p[:5], "Anzahl NaN", (~np.isfinite(p)).sum()
        if (~np.isfinite(p)).sum()>0:
            print "Fehler", 10000
            return 10000
            #raise ValueError("No NaN allowed for CNN parameters")
        for i_n in range(len(inputs)):
            cnn = CNN2d(inputs[i_n],states[i_n],ff_template=p[0:9].reshape((3,3)),fb_template=p[9:18].reshape((3,3)),bias=p[18])
            #print "a,b,z", cnn.a, cnn.b, cnn.z
            #print "Vor integration: ", cnn.output.mean()
            cnn.integrate([0,self.t_sim/2,self.t_sim])
            #print "a,b,z", cnn.a, cnn.b, cnn.z
            #print "CNN output", i_n, cnn.output.mean()
            if self.errfunc=="mean_abs_diff":
                err += abs(refs[i_n]-cnn.output).mean()
            elif self.errfunc=="sum_sqr_diff":
                err += ((refs[i_n]-cnn.output)**2).sum()
            #print err
        print "Fehler", err
        return err
            
    
if __name__ == "__main__":
#if True: # __name__ == "__main__":
    import time
    from eegpy.filter.smoothing import smooth
    #fb = [[0,0,0,],[0,2.0,1.0],[0,0,0]]
    #ff = [[0,0,0,],[0,1.0,0.0],[0,0,0]]
    #bias = 0
    ##myc = CNN2d(np.diag(np.ones((10))),np.zeros((10,10)),ff,fb,1)
    #t0 = time.time()
    #myc = CNN2d(np.diag(np.ones((20))),np.zeros((20,20)),ff,fb,bias)
    #rv = myc.integrate(ts=np.arange(1,26,1))
    #t1 = time.time()-t0
    #print "odeint:", t1, " sec."
    #t0 = time.time()
    #myc2 = CNN2d_rk4(np.diag(np.ones((20))),np.zeros((20,20)),ff,fb,bias)
    #rv2 = []
    #for i in range(25):
    #    rv2.append(myc2.integrate())
    #t2 = time.time()-t0
    #print "rk4:", t2, " sec."
    #p.ioff()
    #p.figure(1)
    #for i in range(25):
    #    p.subplot(2,25,i+1)
    #    p.imshow(myc.sigmoid(rv[i].reshape(myc.x.shape)),interpolation="nearest",vmin=-1,vmax=1)
    #    p.subplot(2,25,i+1+25)
    #    p.imshow(rv2[i],interpolation="nearest",vmin=-1,vmax=1)
    #p.show()
       # time.sleep(0.2)

    inputs = []
    states = []
    refs = []
    for i in range(5):
        #inputs.append(np.diag(np.ones((10))*np.random.random()))
        inputs.append(np.random.random((10,10)))
        #states.append(np.zeros((10,10)))
        states.append(inputs[-1])
        refs.append((-1)*np.ones((10,10)))
    for i in range(5):
        inputs.append(np.random.random((10,10)))
        #inputs.append(np.diag(np.ones((10))*np.random.random()))
        inputs[-1][:] = smooth(inputs[-1],9)#+0.1
        #states.append(np.zeros((10,10)))
        states.append(inputs[-1])
        refs.append(np.ones((10,10))*(1))
    ct = CNNTrainer(t_sim=40,errfunc="mean_abs_diff", opt_method="simplex")
    #ct = CNNTrainer(t_sim=40,errfunc="sum_sqr_diff", opt_method="powell")
    #ct = CNNTrainer(t_sim=40,errfunc="mean_abs_diff", opt_method="anneal")
    a,b,z = ct(inputs,states,refs)
    #p.ioff()
    #print "Test"
    p.hot()
    p.figure(1)
    for i in range(10):
        print i,
        cnn = CNN2d(inputs[i],states[i],a,b,z)
        try:
            o = cnn.integrate(np.arange(1,21,1))
            #p.figure(i+2)
            #for i2 in range(20):
            #    p.subplot(4,5,i2+1)
            #    #p.axes(frameon=False)
            #    p.imshow(cnn.sigmoid(o[i2].reshape((10,10))),interpolation="nearest",vmin=-1,vmax=1)
            #    p.show()
            #assert 1==0
            time.sleep(1)
            p.ioff()
            p.figure(1)
            p.subplot(3,10,i+1)
            p.imshow(inputs[i],interpolation="nearest",vmin=-1,vmax=1)
            p.subplot(3,10,i+1+10)
            p.imshow(refs[i],interpolation="nearest",vmin=-1,vmax=1)
            p.subplot(3,10,i+1+20)
            p.imshow(cnn.output,interpolation="nearest",vmin=-1,vmax=1)
            p.title("%.1f"%cnn.output.mean())
            print "Fertig"
        except Exception, e:
            pass

    p.show()


