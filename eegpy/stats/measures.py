#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.stats import ks_2samp, ksprob, f_oneway, scoreatpercentile, fprob
import numpy as np

def pairwise_ks(*args,**kwargs):#, combine = np.median):
    """For each feature, do a Kolm.-Smirn.-Test for every pair of groups. 
    :Parameters:
      *args: arrays
        Some 1d-arrays, each representing one group
      combine
    Combine the results somehow.
    """
    
    #combine = np.median
    combine = lambda x: scoreatpercentile(x,80)
    if kwargs.has_key("combine"):
        combine = kwargs["combine"]

    num_grps = len(args)
    
    ks = np.zeros((num_grps*(num_grps+1)/2),"d")
    ps = ks.copy()
    
    offset=0
    for i in range(num_grps):
        #print i
        for j in range(i):
            ks[offset], ps[offset] = ks_2samp(args[i],args[j])
            offset+=1
    
    #Many different ways of combining are imaginable, e.g. mean,median, scoreatpercentile...
    result = combine(ks)
    return result

def f_oneway_repeated_measures(M):
    """Calculate One-Way ANOVA for repeated measures. 
    Models the difference between 'subjects' as random effect.
    Example code from Roger Lew:
    ---------------
    import numpy as np
    from scipy.stats import fprob

    # M contains subjects as rows and conditions as columns
    M=np.array([[21,22,8,6,6],
                [20,19,10,4,4],
                [17,15,5,4,5],
                [25,30,13,12,17],
                [30,27,13,8,6],
                [19,27,8,7,4],
                [26,16,5,2,5],
                [17,18,8,1,5],
                [26,24,14,8,9]],dtype='float')

    mu=np.mean(M)
    SS_total=np.sum([[(v-mu)**2 for v in row] for row in M])
    SS_subjects=np.sum([ M.shape[1]*(np.mean(row)-mu)**2 for row in M])
    SS_conditions=np.sum([ M.shape[0]*(np.mean(row)-mu)**2 for row in M.T])
    SS_error=SS_total-SS_subjects-SS_conditions

    df_total=M.size-1
    df_conditions=M.shape[1]-1
    df_subjects=M.shape[0]-1
    df_error=df_total-df_subjects-df_conditions

    F=(SS_conditions/df_conditions)/(SS_error/df_error)
    p=fprob(df_conditions,df_error,F)
    print F,p
    ------------------
    :Parameters:
      M: array-like, 2d
        The array containing the data, 2d, subjects as rows, 
        conditions as columns
    """
    mu=np.mean(M)
    SS_total=np.sum([[(v-mu)**2 for v in row] for row in M])
    SS_subjects=np.sum([ M.shape[1]*(np.mean(row)-mu)**2 for row in M])
    SS_conditions=np.sum([ M.shape[0]*(np.mean(row)-mu)**2 for row in M.T])
    SS_error=SS_total-SS_subjects-SS_conditions

    df_total=M.size-1
    df_conditions=M.shape[1]-1
    df_subjects=M.shape[0]-1
    df_error=df_total-df_subjects-df_conditions

    F=(SS_conditions/df_conditions)/(SS_error/df_error)
    p=fprob(df_conditions,df_error,F)
    return F,p

def ks_2samp(data1, data2):
    """
    Computes the Kolmogorov-Smirnof statistic on 2 samples.

    This is a two-sided test for the null hypothesis that 2 independent samples
    are drawn from the same continuous distribution.

    TAK: taken from scipy.stats. Extended to return the value for which
    the KS-statistic is maximal.

    Parameters
    ----------
    a, b : sequence of 1-D ndarrays
        two arrays of sample observations assumed to be drawn from a continuous
        distribution, sample sizes can be different


    Returns
    -------
    D : float
        KS statistic
    p-value : float
        two-tailed p-value
    x_D: float
        value of the control-variable for which the KS statistic is maximal


    Notes
    -----

    This tests whether 2 samples are drawn from the same distribution. Note
    that, like in the case of the one-sample K-S test, the distribution is
    assumed to be continuous.

    This is the two-sided test, one-sided tests are not implemented.
    The test uses the two-sided asymptotic Kolmogorov-Smirnov distribution.

    If the K-S statistic is small or the p-value is high, then we cannot
    reject the hypothesis that the distributions of the two samples
    are the same.

    Examples
    --------

    >>> from scipy import stats
    >>> import numpy as np
    >>> from scipy.stats import ks_2samp

    >>> #fix random seed to get the same result
    >>> np.random.seed(12345678);

    >>> n1 = 200  # size of first sample
    >>> n2 = 300  # size of second sample

    different distribution
    we can reject the null hypothesis since the pvalue is below 1%

    >>> rvs1 = stats.norm.rvs(size=n1,loc=0.,scale=1);
    >>> rvs2 = stats.norm.rvs(size=n2,loc=0.5,scale=1.5)
    >>> ks_2samp(rvs1,rvs2)
    (0.20833333333333337, 4.6674975515806989e-005)

    slightly different distribution
    we cannot reject the null hypothesis at a 10% or lower alpha since
    the pvalue at 0.144 is higher than 10%

    >>> rvs3 = stats.norm.rvs(size=n2,loc=0.01,scale=1.0)
    >>> ks_2samp(rvs1,rvs3)
    (0.10333333333333333, 0.14498781825751686)

    identical distribution
    we cannot reject the null hypothesis since the pvalue is high, 41%

    >>> rvs4 = stats.norm.rvs(size=n2,loc=0.0,scale=1.0)
    >>> ks_2samp(rvs1,rvs4)
    (0.07999999999999996, 0.41126949729859719)

    """
    data1, data2 = map(np.asarray, (data1, data2))
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    n1 = len(data1)
    n2 = len(data2)
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    data_all = np.concatenate([data1,data2])
    cdf1 = np.searchsorted(data1,data_all,side='right')/(1.0*n1)
    cdf2 = (np.searchsorted(data2,data_all,side='right'))/(1.0*n2)
    d = np.max(np.absolute(cdf1-cdf2))
    x_d = data_all[np.argmax(np.absolute(cdf1-cdf2))]
    #Note: d absolute not signed distance
    en = np.sqrt(n1*n2/float(n1+n2))
    try:
        prob = ksprob((en+0.12+0.11/en)*d)
    except:
        prob = 1.0
    return d, prob, x_d


if __name__ == "__main__":
    #from mvpa.datasets import Dataset
    import pylab as p
    #ngrps = 15
    #npergrp = 20
   # 
    #testdata = np.random.normal(size=(ngrps*npergrp,1200))
    ##First 200: All differ
    #for i in range(ngrps):
    #    testdata[i*npergrp:(i+1)*npergrp,:200] += np.ones((npergrp,200))*(i-ngrps/2)*3
    ##200:400: only 1 is different
    #testdata[:npergrp,200:400] += np.ones((20,200))*3
    ##400:600: no difference
   # 
    ##600:800: smaller difference
    #for i in range(ngrps):
    #    testdata[i*npergrp:(i+1)*npergrp,600:800] += np.ones((npergrp,200))*(i-ngrps/2)*1
    ##800:1000: even smaller difference
    #for i in range(ngrps):
    #    testdata[i*npergrp:(i+1)*npergrp,800:1000] += np.ones((npergrp,200))*(i-ngrps/2)*0.1
    ##800:1000: every second grp has a small difference
    #for i in range(ngrps):
    #    if i%2==0:
    #        testdata[i*npergrp:(i+1)*npergrp,1000:] += np.ones((npergrp,200))*(i-ngrps/2)*0.5
    #    
    #p.subplot(411)
    #for i in range(ngrps):
    #    p.plot(testdata[i*npergrp:(i+1)*npergrp,:].mean(axis=0))
   # 
    #k1s = np.zeros((testdata.shape[1]),"d")
    #k2s = np.zeros((testdata.shape[1]),"d")
    #fs = np.zeros((testdata.shape[1]),"d")
    #for i in range(testdata.shape[1]):
    #    grps = [testdata[j*npergrp:(j+1)*npergrp,i]for j in range(ngrps)]
    #    k1s[i] = pairwise_ks(*grps)
    #    #k2s[i] = pairwise_ks(*grps)
    #    k2s[i] = pairwise_ks(*grps,combine=np.mean)
    #    fs[i] = f_oneway(*grps)[0]
   # 
    #p.subplot(412)
    #p.plot(k1s)
    #p.subplot(413)
    #p.plot(k2s)
    #p.subplot(414)
    #p.plot(fs)
    #fwm = PairwiseKS(lambda x: scipy.stats.scoreatpercentile(x,20))
    #sens = fwm(dataset)
    #p.subplot(414)
    #p.plot(sens)
    
    
    
    #fsel = EffectAndDifferenceFeatureSelection(n_select=0.30)
    
    #selection = fsel(dataset)
    #print fsel._selector.felements
    #print selection[0].samples.shape
    #p.plot(fsel.sensitivity)
    #p.plot(fsel.selected_ids,selection, "rD")
    #p.plot(selection, "rD")
    #p.show()

    #Test ks_2samp
    for offset in np.arange(0,10.0,1.):
        ar1 = np.random.normal(0,1,(10000))
        ar2 = np.random.normal(offset,1,(20000))
        d,prob,x_d = ks_2samp(ar1,ar2)
        print offset,d,prob,x_d
        p.subplot(10,1,int(offset/1.)+1)
        p.hist(ar1,50,normed=True,ec="b",fc="none")
        p.hist(ar2,50,normed=True,ec="g",fc="none")
        p.axvline(x_d,color="r")

