#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import mvpa
from mvpa.clfs.libsvmc.svm import SVM
from mvpa.clfs.meta import FeatureSelectionClassifier
from mvpa.featsel.base import SensitivityBasedFeatureSelection
from mvpa.featsel.helpers import FixedNElementTailSelector, FractionTailSelector
import mvpa.misc.transformers
from mvpa.misc.transformers import Absolute
from mvpa.measures.anova import OneWayAnova

import numpy as N
np = N

"""My classifier-warehouse"""
__all__ = ["calc_cm_in_group"]

def calc_cm_in_group(cm, groups):
    """For a given confusion matrix, calculate how many predictions where in-group.
    Ignores all labels not within any group."""
    #First, clean up. Check all labels if they exist in cm
    groups = copy.copy(groups)
    for g in groups:
        for i_l in range(len(g)-1,-1,-1):
            if not g[i_l] in cm.labels:
                print "Ignoring label", g[i_l]
                del g[i_l]
    #Now calculate
    ingrp = 0
    for g in groups:
        for l1 in g: 
            for l2 in g:
                try:
                    ingrp += cm.matrix[cm.labels.index(l1),cm.labels.index(l2)]
                except ValueError, ve:
                    print "ValueError: label %i or %i not in ConfusionMatrix" %(l1,l2)
    all_labels = np.concatenate(groups)
    all_predictions = cm.matrix[[cm.labels.index(l) for l in all_labels]][:,[cm.labels.index(l) for l in all_labels]].sum()
    return float(ingrp)/all_predictions, ingrp
    
