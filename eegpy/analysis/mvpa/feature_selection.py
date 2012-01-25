#!/usr/bin/env python
# -*- coding: utf-8 -*-

from random import shuffle

from mvpa.measures.base import FeaturewiseDatasetMeasure, CombinedFeaturewiseDatasetMeasure
from mvpa.measures.anova import OneWayAnova
from mvpa.misc.transformers import FirstAxisMean, SecondAxisSumOfAbs, Absolute
from mvpa.misc.state import StateVariable, ClassWithCollections 
from mvpa.featsel.base import FeatureSelection, SensitivityBasedFeatureSelection, CombinedFeatureSelection
from mvpa.featsel.helpers import FixedNElementTailSelector, FractionTailSelector
from mvpa.clfs.svm import LinearCSVMC
from mvpa.suite import Dataset

import numpy as np
from scipy.stats import ttest_1samp, fprob, ttest_ind, ks_2samp, f_oneway

"""Define some custom classes for Feature-Selection with PyMVPA"""

class TTestFeaturewiseMeasure(FeaturewiseDatasetMeasure):
    """For each feature, test for the null-hypothesis that the mean is equal to a certain value"""
    
    def __init__(self, val=0, combiner=SecondAxisSumOfAbs, **kwargs):
        """Initialize 
   
          :Parameters:
            val : float
              Real-valued number for the null-hypothesis.
            combiner : Functor 
              The combiner is only applied if the computed featurewise dataset 
              measure is more than one-dimensional. This is different from a 
              `transformer`, which is always applied. By default, the sum of 
              absolute values along the second axis is computed. 
        """
        FeaturewiseDatasetMeasure.__init__(self,combiner,**kwargs)
        self.__val = val
        
    def __call__(self, dataset, labels=None):
        """Actually calculate the p-values."""
        # number of groups 
        if labels is None: 
            labels = dataset.labels
        
        alldata = dataset.samples
        #print alldata.shape
        ps = ttest_1samp(alldata,self.__val,axis=0)[1]
        #print ps.shape, ps
        return ps
        
    @property
    def val(self):
        """Return value for null-hypothesis"""
        return self.__val
    
class OneWayAnovaProbs(FeaturewiseDatasetMeasure):
    """For each feature, test for the null-hypothesis that the mean is equal to a certain value"""
        
    def __call__(self, dataset, labels=None):
        """Actually calculate the p-values."""
        f = OneWayAnova()(dataset)
        # number of groups 
        if labels is None: 
            labels = dataset.labels
        # Calculate degrees of freedom
        ul = np.unique(labels)
        na = len(ul) 
        bign = float(dataset.nsamples)
        dfbn = na-1 
        dfwn = bign - na
        # Now propabilities
        ps = fprob(dfbn,dfwn,f)
        return ps
    
class EffectAndDifferenceFeatureSelection(FeatureSelection):
    """Selects those features that show a significant effect over all samples and differ between groups"""
    def __init__(self,val=0,n_select=0.05):
        FeatureSelection.__init__(self)
        self.__ttfm = TTestFeaturewiseMeasure(val)
        self._combinedSA = CombinedFeaturewiseDatasetMeasure(analyzers=[self.__ttfm,OneWayAnovaProbs()], combiner = np.prod)
        if n_select<1 and n_select>0:
            self._selector=FractionTailSelector(n_select, mode='select', tail='lower')
        elif n_select>1:
            self._selector=FixedNElementTailSelector(int(n_select), mode='select', tail='lower')
        else:
            raise ValueError("n_select must be > 0, is %f"%n_select)
        self.__fs = SensitivityBasedFeatureSelection(
                        self._combinedSA,
                        self._selector,
                        enable_states = ["sensitivity", "selected_ids"],
                    )
        
    def __call__(self,dataset,testdataset=None):
        print "calling EADFS"
        result =  self.__fs(dataset)
        print "Sensitivity", self.__fs.sensitivity
        print result
        return result
    
    #Expose some properties
    sensitivity = property(fget=lambda self:self.__fs.sensitivity)
    sensitivity_analyzer = property(fget=lambda self:self.__fs.sensitivity_analyzer)
    selected_ids = property(fget=lambda self:self.__fs.selected_ids)
    fs = property(fget=lambda self:self.__fs)

class PairwiseKS(FeaturewiseDatasetMeasure):
    """For each feature, do a Kolm.-Smirn.-Test for every pair of labels. 
    Combine the results somehow.
    """
    
    def __init__(self, combine = lambda x: np.median(x,axis=1), **kwargs):
        """Initialize 
        """
        FeaturewiseDatasetMeasure.__init__(self,**kwargs)
        self._combine = combine
        
    def __call__(self, dataset, labels=None):
        """Actually calculate the k-values and do some kind of averaging."""
        # number of groups 
        if labels is None: 
            labels = dataset.labels
        
        ul = np.unique(labels)
        alldata = dataset.samples
        
        ks = np.zeros((alldata.shape[1],(len(ul)*(len(ul)-1))/2),"d")
        ps = ks.copy()
        
        offset=0
        for i in range(len(ul)):
            #print i
            for j in range(i):
                data1 = alldata[labels==ul[i]]
                data2 = alldata[labels==ul[j]]
                for k in range(data1.shape[1]):
                    ks[k,offset], ps[k,offset] = ks_2samp(data1[:,k],data2[:,k])
                offset+=1
        
        #Many different ways of combining are imaginable, e.g. mean,median, scoreatpercentile...
        result = self._combine(ks)
        return result

class GroupByGroupFeatureSelection(FeatureSelection):
    """Feature selection.

    For each feature, calculate 'separability from rest' for each group. 
    Then choose num_features features, balanced for all groups.
    """

    sensitivity = StateVariable(enabled=False)

    def __init__(self,
                 test_type = "ks",
                 num_features = 10,
                 **kwargs
                 ):
        """Initialize feature selection

        :Parameters:
          test_type : 'anova', 'ttest', 'ks', 'svm'
            type of statistical test to perform
          num_features : 
            Given a sensitivity map it has to return the ids of those
            features that should be kept.

        """

        # base init first
        FeatureSelection.__init__(self, **kwargs)

        self.__test_type = test_type
        """Test-type, One of 'anova', 'ttest', 'ks', 'svm'"""

        self.__num_features = num_features
        """Number of features to include"""


    def untrain(self):
        #if __debug__:
        #    debug("FS_", "Untraining sensitivity-based FS: %s" % self)
        #self.__sensitivity_analyzer.untrain()
        pass


    def __call__(self, dataset, testdataset=None):
        """Select the most important features

        :Parameters:
          dataset : Dataset
            used to compute sensitivity maps
          testdataset: Dataset
            optional dataset to select features on

        Returns a tuple of two new datasets with selected feature
        subset of `dataset`.
        """

        tt = self.__test_type
        
        ul = np.unique(dataset.labels)
        
        sensitivity = np.zeros((len(ul),dataset.nfeatures),"d")
        data = dataset.samples.reshape(dataset.samples.shape[0],-1)
        sorted_ids = []
        for i in range(len(ul)):
            d1 = data[dataset.labels==ul[i]]
            d2 = data[dataset.labels!=ul[i]]
            if tt == 'anova':
                for j in range(data.shape[1]):
                    sensitivity[i,j] = f_oneway(d1[:,j],d2[:,j])[0]
            elif tt == 'ttest':
                sensitivity[i,:] = abs(ttest_ind(d1,d2,axis=0)[0])
            elif tt == 'ks':
                for j in range(data.shape[1]):
                    sensitivity[i,j] = ks_2samp(d1[:,j],d2[:,j])[0]
            elif tt == 'svm':
                dataset2 = dataset.copy()
                dataset2.labels[dataset.labels==ul[i]] = [0]*len(dataset2.labels[dataset.labels==ul[i]])
                dataset2.labels[dataset.labels!=ul[i]] = [1]*len(dataset2.labels[dataset.labels!=ul[i]])
                sensana =  LinearCSVMC().getSensitivityAnalyzer(transformer=Absolute)
                sensitivity[i,:] = sensana(dataset2)
            #Sorted list of feature_ids for each  group
            ids = range(dataset.nfeatures)
            ids.sort(cmp=lambda x,y: int(np.sign(sensitivity[i,y]-sensitivity[i,x])))
            #print "ids", i, ids, sensitivity.mean()
            #print ["%.2f"%x for x in sensitivity[i,ids]]
            sorted_ids.append(ids)
            
        self.sensitivity = sensitivity

        # Select features to preserve
        sid = []
        grp_order = range(len(ul))
        shuffle(grp_order)
        offsets = np.zeros((len(ul)),np.int32)
        #print "Sensitivity: min:%.2f max:%.2f"%(sensitivity.min(),sensitivity.max()),
        while len(sid) < self.__num_features:
            for i in grp_order:
                if len(sid) < self.__num_features:
                    while offsets[i]<dataset.nfeatures:
                        if not sorted_ids[i][offsets[i]] in sid:
                            #print "Chose feature", sorted_ids[i][offsets[i]], "from group", i, "offset", offsets[i]
                            #print "Grp.", i, "feat.", sorted_ids[i][offsets[i]], "sens.", sensitivity[i,sorted_ids[i][offsets[i]]] 
                            sid.append(sorted_ids[i][offsets[i]])
                            offsets[i]+=1
                            break
                        else:
                            offsets[i]+=1
        #sid.sort()
        print sid, len(sid), "selected ids"
        selected_ids = sid 

        #if __debug__:
        #    debug("FS_", "Sensitivity: %s Selected ids: %s" %
        #          (sensitivity, selected_ids))

        # Create a dataset only with selected features
        wdataset = dataset.selectFeatures(selected_ids)

        if not testdataset is None:
            wtestdataset = testdataset.selectFeatures(selected_ids)
        else:
            wtestdataset = None

        # Differ from the order in RFE when actually error reported is for
        results = (wdataset, wtestdataset)

        # WARNING: THIS MUST BE THE LAST THING TO DO ON selected_ids
        selected_ids.sort()
        self.selected_ids = np.array(selected_ids)

        # dataset with selected features is returned
        return results

if __name__ == "__main__":
    from mvpa.datasets import Dataset
    import pylab as p
    ngrps = 15
    npergrp = 20
    
    testdata = np.random.normal(size=(ngrps*npergrp,1200))
    #First 200: All differ
    for i in range(ngrps):
        testdata[i*npergrp:(i+1)*npergrp,:200] += np.ones((npergrp,200))*(i-ngrps/2)*3
    #200:400: only 1 is different
    testdata[:npergrp,200:400] += np.ones((20,200))*3
    #400:600: no difference
    
    #600:800: smaller difference
    for i in range(ngrps):
        testdata[i*npergrp:(i+1)*npergrp,600:800] += np.ones((npergrp,200))*(i-ngrps/2)*1
    #800:1000: even smaller difference
    for i in range(ngrps):
        testdata[i*npergrp:(i+1)*npergrp,800:1000] += np.ones((npergrp,200))*(i-ngrps/2)*0.1
    #800:1000: every second grp has a small difference
    for i in range(ngrps):
        if i%2==0:
            testdata[i*npergrp:(i+1)*npergrp,1000:] += np.ones((npergrp,200))*(i-ngrps/2)*0.5
            
    #Ein paar Ausreisser
    for i in [250,500,700,750]:
        testdata[(ngrps/2)*npergrp:(ngrps/2+1)*npergrp,i] += 50
        
    p.subplot(411)
    for i in range(ngrps):
        p.plot(testdata[i*npergrp:(i+1)*npergrp,:].mean(axis=0))
    #Labels
    labels=[]
    for i in range(ngrps):
        labels+=[i]*npergrp
    dataset = Dataset(samples=testdata, labels=labels)
    print testdata.shape, dataset.samples.shape
    
    fwm = PairwiseKS(lambda x: np.median(x,axis=1))
    sens = fwm(dataset)
    p.subplot(412)
    p.plot(sens)
    fwm = PairwiseKS(lambda x: np.mean(x,axis=1))
    sens = fwm(dataset)
    p.subplot(413)
    p.plot(sens)
    fwm = OneWayAnova()
    sens = fwm(dataset)
    #p.subplot(414)
    #p.plot(sens)
    
    fs = GroupByGroupFeatureSelection("svm",200,enable_states=["selected_ids","sensitivity"])
    fs(dataset)
    p.subplot(414)
    p.plot(fs.selected_ids,np.zeros((fs.selected_ids.shape),"d"),"|r",ms=4)
    for i in range(fs.sensitivity.shape[0]):
        p.plot(fs.sensitivity[i,:])
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
    p.show()
    
    
