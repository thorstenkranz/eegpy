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

"""My classifier-warehouse"""
__all__ = ["some_svms","svms_for_CombinedClassifier"]

def some_svms():
    """Returns a couple of FeatureSelectionClassifiers
    based on SVMs with different numbers of features and/or
    sensitivity measure"""
    clfr1 = FeatureSelectionClassifier(
          SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
          SensitivityBasedFeatureSelection(
                OneWayAnova(),
                FixedNElementTailSelector(500, mode='select', tail='upper')
            ),
              descr="LinSVM on 500 (ANOVA)"
          )
    clfr2 = FeatureSelectionClassifier(
              SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
              SensitivityBasedFeatureSelection(
                SVM().getSensitivityAnalyzer(transformer=Absolute),
                FixedNElementTailSelector(500, mode='select', tail='upper')
            ),
              descr="LinSVM on 500 (SVM)"
          )
    clfr3 = SVM()
    clfr4 = FeatureSelectionClassifier(
              SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
              SensitivityBasedFeatureSelection(
                SVM().getSensitivityAnalyzer(transformer=Absolute),
                FractionTailSelector(0.05, mode='select', tail='upper'),
            ),
              descr="LinSVM on 5 % (SVM)"
          )
    return [clfr1,clfr2,clfr3,clfr3]


def svms_for_CombinedClassifier():
    """For my iEEG study, I use a CombinedClassifier. The components are defined here"""
    clfrs = []
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    #SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    OneWayAnova(),
                    FixedNElementTailSelector(500, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 500 (Anova)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    #SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    OneWayAnova(),
                    FixedNElementTailSelector(300, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 300 (Anova)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    #SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    OneWayAnova(),
                    FixedNElementTailSelector(200, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 200 (Anova)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    #SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    OneWayAnova(),
                    FixedNElementTailSelector(500, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 100 (Anova)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    #OneWayAnova(),
                    FixedNElementTailSelector(500, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 500 (SVM)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    #OneWayAnova(),
                    FixedNElementTailSelector(300, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 300 (SVM)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    #OneWayAnova(),
                    FixedNElementTailSelector(200, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 200 (SVM)"
                 )
                )
    clfrs.append(FeatureSelectionClassifier(
                  SVM(descr = "libsvm.LinSVM(C=def)", probability = 1),
                   SensitivityBasedFeatureSelection(
                    SVM(descr = "libsvm.LinSVM(C=def)", probability = 1).getSensitivityAnalyzer(transformer=mvpa.misc.transformers.Absolute),
                    #OneWayAnova(),
                    FixedNElementTailSelector(500, mode='select', tail='upper')
                   ),
                 descr="LinSVM on 100 (SVM)"
                 )
                )
    return clfrs


