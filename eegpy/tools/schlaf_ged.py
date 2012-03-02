#!/usr/bin/env python
# -*- coding: utf-8 -*-

import urllib2
from eegpy.misc import EegpyBase
import numpy as np
#from nitime.utils import zscore

def get_good_channels(patient_name):
    """Reads from my website the list of good channels for a patient"""
    address = "http://www.thorstenkranz.de/SchlafGed_good_channels/good_channels.php?patient=%s"
    opener = urllib2.build_opener()
    page = opener.open(address%patient_name)
    rv = []
    try:
        rv = eval(page.read())
    except SyntaxError,se:
        if patient_name.endswith("_bp"):
            try:
                page = opener.open(address%patient_name[:-3])
                rv = eval(page.read())
            except SyntaxError,se:
                raise ValueError("No good channels. Probably patient not found")
        else:
            raise ValueError("No good channels. Probably patient not found")
    return rv

class PredictNightsResult(EegpyBase):
    """Object for analyzing results from classifier-predictions for nights"""
    def __init__(self,num_classifiers,num_wins_n1,num_wins_n2):
        self.log = ""
        self.cv_confusion = None # Confusion matrix from CrossValidation on the experiment-data
        self.num_classifiers = num_classifiers
        self.num_wins_n1 = num_wins_n1
        self.num_hits_n1 = 0
        self.vote_n1 = np.zeros((num_wins_n1),"d")
        self.predictions_n1 = np.zeros((num_wins_n1),"d")
        self.artifacts_n1 = np.ones((num_wins_n1),"d")
        self.num_wins_n2 = num_wins_n2
        self.num_hits_n2 = 0
        self.vote_n2 = np.zeros((num_wins_n1),"d")
        self.predictions_n2 = np.zeros((num_wins_n1),"d")
        self.artifacts_n2 = np.ones((num_wins_n1),"d")

    def num_hits_per_label(self,night=1,thres=None):
        """Returns hits for each label. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        #Assign local variables according to night chosen
        if night==1:
            predictions = self.predictions_n1
            vote = self.vote_n1
            artifacts = self.artifacts_n1
            num_hits = self.num_hits_n1
            num_wins = self.num_wins_n1
        elif night==2:
            predictions = self.predictions_n2
            vote = self.vote_n2
            artifacts = self.artifacts_n2
            num_hits = self.num_hits_n2
            num_wins = self.num_wins_n2
        #up = np.unique() # unique predictions
        #choose threshold automatically?
        if thres==None:
            thres==self.num_classifiers
        rv = {}
        for pred in predictions[vote>thres]:
            if not rv.has_key(pred):
                rv[pred]=0
            rv[pred]+=1
        return rv

    def num_hits(self,night=1,thres=None):
        """Returns total number of hits. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        hits_per_label = self.num_hits_per_label(night,thres)
        return np.sum(hits_per_label.values())
        
class PredictNightsResultSingleClassifier(EegpyBase):
    """Object for analyzing results from classifier-predictions for nights.
    In contrast to PredictNightsResult, used for only ONE classifier.
    It stores the probabilities for each class for each prediction"""
    def __init__(self,num_classes,num_wins_n1,num_wins_n2):
        self.log = ""
        self.cv_confusion = None # Confusion matrix from CrossValidation on the experiment-data
        self.num_classes = num_classes
        self.num_wins_n1 = num_wins_n1
        #self.num_hits_n1 = 0
        self.vote_n1 = np.zeros((num_wins_n1, num_classes),"d")
        self.predictions_n1 = np.zeros((num_wins_n1),"d")
        self.artifacts_n1 = np.ones((num_wins_n1),"d")
        self.num_wins_n2 = num_wins_n2
        #self.num_hits_n2 = 0
        self.vote_n2 = np.zeros((num_wins_n2,num_classes),"d")
        self.predictions_n2 = np.zeros((num_wins_n2),"d")
        self.artifacts_n2 = np.ones((num_wins_n2),"d")

    def num_hits_per_label(self,night=1,thres=None,zscored=True):
        """Returns hits for each label. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        #Assign local variables according to night chosen
        if night==1:
            predictions = self.predictions_n1
            vote = self.vote_n1
            #artifacts = self.artifacts_n1
            #num_hits = self.num_hits_n1
            #num_wins = self.num_wins_n1
        elif night==2:
            predictions = self.predictions_n2
            vote = self.vote_n2
            #artifacts = self.artifacts_n2
            #num_hits = self.num_hits_n2
            #num_wins = self.num_wins_n2
        #up = np.unique() # unique predictions
        #choose threshold automatically?
        if not zscored:
            #take max value
            vote = vote.max(axis=1)
        else:
            #zscoring of only max-value
            vote = (vote.max(axis=1) - vote.mean(axis=1)) / vote.std(axis=1)
        
        #if thres==None:
        #    thres==self.num_classifiers
        rv = {}
        for pred in predictions[vote>thres]:
            if not rv.has_key(pred):
                rv[pred]=0
            rv[pred]+=1
        return rv

    def num_hits(self,night=1,thres=None,zscored=True):
        """Returns total number of hits. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        hits_per_label = self.num_hits_per_label(night,thres,zscored)
        real_hits = [hits_per_label[k] for k in hits_per_label.keys() if k != 0.0]
        return np.sum(real_hits)

class PredictNightsResultWithSurrogates(EegpyBase):
    """Object for analyzing results from classifier-predictions for nights.
    In contrast to PredictNightsResult, used for only ONE classifier.
    It stores the probabilities for each class for each prediction
    Offers infrastructure to analyze data for surrogates as well"""
    def __init__(self,num_classes,num_surrogates,num_wins_n1,num_wins_n2,art_thr=7.0,art_thr_diff=7.5):
        self.log = ""
        self.cv_confusions = None # Confusion matrix from CrossValidation on the experiment-data
        self.num_classes = num_classes
        self.num_wins_n1 = num_wins_n1
        self.num_surrogates = num_surrogates
        #self.num_hits_n1 = 0
        self.vote_n1 = np.zeros((num_wins_n1,num_classes, num_surrogates+1),"d")
        self.predictions_n1 = np.zeros((num_wins_n1, num_surrogates+1),"d")
        self.artifacts_n1 = np.ones((num_wins_n1,3),"d")
        self.num_wins_n2 = num_wins_n2
        #self.num_hits_n2 = 0
        self.vote_n2 = np.zeros((num_wins_n2,num_classes, num_surrogates+1),"d")
        self.predictions_n2 = np.zeros((num_wins_n2, num_surrogates+1),"d")
        self.artifacts_n2 = np.ones((num_wins_n2,3),"d")
        self.art_thr = art_thr
        self.art_thr_diff = art_thr_diff

    def no_artifact(self,night):
        if night==1:
            an = self.artifacts_n1
        elif night==2:
            an = self.artifacts_n2
        else:
            raise ValueError("Invalid value for night")
        return (an[:,1]<self.art_thr) * (an[:,2]<self.art_thr_diff)

    def num_hits_per_label(self,night=1,thres=0.7,zscored=True):
        """Returns hits for each label. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        #Assign local variables according to night chosen
        if night==1:
            predictions = self.predictions_n1
            vote = self.vote_n1
            #artifacts = self.artifacts_n1
            #num_hits = self.num_hits_n1
            #num_wins = self.num_wins_n1
        elif night==2:
            predictions = self.predictions_n2
            vote = self.vote_n2
            #artifacts = self.artifacts_n2
            #num_hits = self.num_hits_n2
            #num_wins = self.num_wins_n2
        #up = np.unique() # unique predictions
        #choose threshold automatically?
        if not zscored:
            #take max value
            vote = vote.max(axis=1)
        else:
            #zscoring of only max-value
            vote = (vote.max(axis=1) - vote.mean(axis=1)) / vote.std(axis=1)
        
        #if thres==None:
        #    thres==self.num_classifiers
        rv = [{} for i in range(self.num_surrogates+1)]
        for i_surr in range(self.num_surrogates+1):
            tmp_pred = predictions[self.no_artifact(night),i_surr]
            tmp_vote = vote[self.no_artifact(night),i_surr]
            for pred in tmp_pred[tmp_vote>thres]:
                if not rv[i_surr].has_key(pred):
                    rv[i_surr][pred]=0
                rv[i_surr][pred]+=1
        return rv

    def num_hits(self,night=1,thres=0.7,zscored=True):
        """Returns total number of hits. 
        Threshold for this decision is thres. If None, take number of classifiers"""
        hits_per_label = self.num_hits_per_label(night,thres,zscored)
        rv = np.zeros((len(hits_per_label)))
        for i, d in enumerate(hits_per_label):
            keys = d.keys()
            keys.sort()
            rv = [np.sum(d[k]) for k in keys]
        return rv

    @property
    def num_surr(self):
        return self.num_surrogates

if __name__=="__main__":
    import shelve
    shv = shelve.open("/media/story/SchlafGed/iEEG/predict_nights_results/canseven_nacht_20100131.shv")
    pn_result = PredictNightsResult(shv["vote_n2"].max(),shv["vote_n1"].shape[0],shv["vote_n2"].shape[0])
    pn_result.num_hits_n1 = shv["num_hits_n1"]
    pn_result.vote_n1 = shv["vote_n1"]
    pn_result.predictions_n1 = shv["predictions_n1"]
    pn_result.artifacts_n1 = shv["artifacts_n1"]
    pn_result.num_hits_n2 = shv["num_hits_n2"]
    pn_result.vote_n2 = shv["vote_n2"]
    pn_result.predictions_n2 = shv["predictions_n2"]
    pn_result.artifacts_n2 = shv["artifacts_n2"]




