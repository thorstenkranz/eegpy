#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import mvpa
from mvpa.datasets.splitters import Splitter

import numpy as N

"""Define some custom Splitter"""

class OnePerLabelSplitter(Splitter):
    """For each label, put one sample for validation in dataset 2, leave equal amounts of the others in dataset 1.
    I.e. for 5 labels with 4,7,5,6,9 samples, create 4 splits."""
    
    def __init__(self, strategy="random", **kwargs):
        """Initialize the OnePerLabelSplitter
   
          :Parameters:
            kwargs
              additional parameters are passed to the 'Splitter' base class
        """
        Splitter.__init__(self,**kwargs)
        self._strategy = strategy

    #Override method from Splitter
    def splitcfg(self, dataset):
        """Return splitcfg for a given dataset"""
        self._ds_labels = dataset.labels
        #return self._getSplitConfig(eval('dataset.unique' + self.__splitattr))
        return self._getSplitConfig(dataset.uniquechunks)

    def _getSplitConfig(self,uniqueattrs):
        """Huka chaka, wuka waka! (TAK: Whatever that means...)
        """
        uniquelabels = N.unique(self._ds_labels)
        labelpositions = [uniqueattrs[self._ds_labels==l] for l in uniquelabels]
        min_l = min([len(lp) for lp in labelpositions])
        if not min_l > 0:
            print "Each label must be represented at least once"
            raise ValueError("Each label must be represented at least once")
        if self._strategy == "random":
            # Shuffle positions
            for lp in labelpositions:
                N.random.shuffle(lp)
        else:
            pass
        # Only pick the min_l first labelpositions
        labelpositions = [lp[:min_l] for lp in labelpositions]
        for lp in labelpositions:
            N.sort(lp)
        #print "OPLS:", [len(lp) for lp in labelpositions], labelpositions
        all_positions = N.concatenate(labelpositions,axis=0)
        #print "OPLS:", all_positions
        split_list = []
        for i in range(min_l):
            #split1 = copy.deepcopy(labelpositions)
            #split2 = [lp[i] for lp in split1]
            #for j in range(len(split1)):
            #    #print j, split1[j].shape,
            #    split1[j] = N.delete(split1[j],i)
            #    #print split1[j].shape
            #split_list.append((N.concatenate(split1,axis=0),split2))
            split2 = [lp[i] for lp in labelpositions]
            split1 = list(set(all_positions).difference(set(split2)))
            #print "OPLS:", len(all_positions),len(split1), len(split2)
            split_list.append((split1,split2))
        return split_list

    def __str__(self):
        """String summary over the object
        """
        return "OnePerLabelSplitter "
