import copy

import mvpa
from mvpa.suite import CrossValidatedTransferError, TransferError, ConfusionMatrix, TreeClassifier, Classifier, SplitClassifier, SVM
from mvpa.misc.state import ClassWithCollections
import numpy as np
import hcluster
from eegpy.analysis.mvpa.splitters import OnePerLabelSplitter

class LinkageClassifier(Classifier):
    def __init__(self,clf,splitter=None, splitter_clf=None,*args,**kwargs):
        Classifier.__init__(self,*args,**kwargs)
        self._clf = clf
        self._confusion = None
        self._dataset = None
        self.ulabels = None
        if splitter == None:
            self._splitter = OnePerLabelSplitter()
        else:
            self._splitter = splitter
        if splitter_clf==None:
            self._splitter_clf=self._splitter
        else:
            self._splitter_clf==splitter_clf
        self.linkage = None
        self._tree = None
        self._tree_clf = None

    def _train(self, trainset):
        self._dataset = trainset
        self.ulabels=trainset.uniquelabels
        # Do cross-validation for normal classifier
        self.cvterr = CrossValidatedTransferError(TransferError(self._clf),self._splitter,enable_states=["confusion"])
        self.cvterr(self._dataset)
        # From the confusion matrix, calculate linkage and tree-structure
        # First prepare distance matrix from confusion matrix
        dist = self.cvterr.confusion.matrix
        dist = dist.max()-dist # Kind of inversion. High values in confusion -> similar -> small distance
        dist = (dist+dist.T)/2 # Distance must be symmetric (property of a norm)
        dist -= np.diag(np.diag(dist)) # Distance to self must be zero -> make diagonal elements zero
        # Calculate linkage matrix
        self.linkage = hcluster.linkage(hcluster.squareform(dist))
        # Build tree and according TreeClassifier
        self.tree = hcluster.to_tree(self.linkage)
        self._tree_clf = self.build_tree_classifier_from_linkage_tree(self.tree)[0]
        self._tree_clf.train(trainset)
        #print "Trained on", self.ulabels

    def build_tree_classifier_from_linkage_tree(self, tree):
        # Here SVM is only used as it can do single-class classification (pseudo)
        if tree.left.is_leaf() and tree.right.is_leaf():
            return (SplitClassifier(self._clf.clone(),self._splitter_clf),[tree.left.id, tree.right.id])
            #return (self._clf.clone(),[tree.left.id, tree.right.id])
        elif tree.left.is_leaf():
            clf1, ids1 = (SVM(), [tree.left.id])
            clf2, ids2 = self.build_tree_classifier_from_linkage_tree(tree.right)
            return (TreeClassifier(SplitClassifier(self._clf.clone(),self._splitter_clf),
            #return (TreeClassifier(self._clf.clone(),
                                    {"c%02i"%tree.left.id:([self.ulabels[i] for i in ids1],clf1), "c%02i"%tree.right.id:([self.ulabels[i] for i in ids2],clf2)}), 
                    ids1+ids2)
        elif tree.right.is_leaf():
            clf1, ids1 = self.build_tree_classifier_from_linkage_tree(tree.left)
            clf2, ids2 = (SVM(), [tree.right.id])
            return (TreeClassifier(SplitClassifier(self._clf.clone(),self._splitter_clf),
            #return (TreeClassifier(self._clf.clone(),
                                    {"c%02i"%tree.left.id:([self.ulabels[i] for i in ids1],clf1), "c%02i"%tree.right.id:([self.ulabels[i] for i in ids2],clf2)}), 
                    ids1+ids2)
        else:
            clf1, ids1 = self.build_tree_classifier_from_linkage_tree(tree.left)
            clf2, ids2 =  self.build_tree_classifier_from_linkage_tree(tree.right)
            return (TreeClassifier(SplitClassifier(self._clf.clone(),self._splitter_clf),
            #return (TreeClassifier(self._clf.clone(),
                                    {"c%02i"%tree.left.id:([self.ulabels[i] for i in ids1],clf1), "c%02i"%tree.right.id:([self.ulabels[i] for i in ids2],clf2)}), 
                    ids1+ids2)
    
    def untrain(self):
        self._tree_clf = None
        super(LinkageClassifier,self).untrain()

    def _predict(self,samples):
        if self._tree_clf != None:
            return self._tree_clf.predict(samples)
        else:
            raise ValueError("Classifier wasn't yet trained, so cannot predict.")

    def dendrogram(self):
        #import pylab as p
        if not self.linkage == None:
            hcluster.dendrogram(self.linkage,labels=np.unique(self._dataset.labels))

