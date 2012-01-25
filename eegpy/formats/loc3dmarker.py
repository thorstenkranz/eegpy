# -*- coding: utf-8 -*-
"""Read location marker files as exported by Loc3djr and PyLocator"""

import os.path
import sys

import numpy as np

class Loc3dMarkers(object):
    """Handle information about markers from Loc3djr."""
    def __init__(self,fn=None):
        self._fn = fn
        if fn!=None:
            self.read_from_file(fn)
        else:
            self._labels = None
            self._xyz = None
            self._sizes = None
            self._colors = None

    def __str__(self):
        return "Loc3dMarkers (%i entries, file: %s)"%(self.count,str(self._fn))

    def __repr__(self):
        return str(self)

    def __add__(self,other):
        rv = Loc3dMarkers()
        rv._labels = np.concatenate([self._labels,other._labels],axis=0)
        rv._xyz = np.concatenate([self._xyz,other._xyz],axis=0)
        rv._sizes = np.concatenate([self._sizes,other._sizes],axis=0)
        rv._colors = np.concatenate([self._colors,other._colors],axis=0)
        return rv

    def read_from_file(self,fn):
        """Reads an ASCII-file as exported from Loc3dJR"""
        fh = open(fn,"r")
        labels = []
        xyz = []
        sizes = []
        colors = []
        for line in fh.readlines():
            try:
                if not line.startswith("#"):
                    label,x,y,z,size,r,g,b = line.split(",")
                    labels.append(label)
                    xyz.append([x,y,z])
                    sizes.append(size)
                    colors.append((float(r),float(g),float(b)))
            except IOError, ioe:
                print "IOError:", ioe
        self._labels = np.array(labels)
        self._xyz = np.array(xyz).astype("f")
        self._sizes = np.array(sizes).astype("f")
        self._colors = np.array(colors)

    @property
    def labels(self):
        return self._labels

    @property
    def xyz(self):
        return self._xyz

    @property
    def sizes(self):
        return self._sizes

    @property
    def colors(self):
        return self._colors

    @property
    def count(self):
        return self.labels.shape[0]


if __name__=="__main__":
    lm = Loc3dMarkers("/media/story/SchlafGed/iEEG/Implantationsdaten/Can7/canseven_electroden.txt")
    print lm.labels
        


        
