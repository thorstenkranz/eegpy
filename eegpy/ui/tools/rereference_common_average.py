#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command-line program for rereferencing eeg to common-average"""

import sys
import getopt
import shutil

import numpy as np

import eegpy

debug = False
input = None
output = None
exclude_channels = []


def main(argv):
    parse_args(sys.argv[1:])
    if debug:
        print input, output, exclude_channels
    #copy input-file to ouput-file
    try:
        shutil.copyfile(input,output)
    except Exception,e:
        print "Cannot copy file %s to %s,"%(input,output), e
        usage()
        sys.exit(2)
    #open eeg-file and rereference
    try:
        eeg = eegpy.open_eeg(output)
        bsl = np.ones((eeg.num_channels),np.bool)
        for ec in exclude_channels:
            if ec<eeg.num_channels:
                bsl[ec]=True
        start = 0
        while start<eeg.num_datapoints:
            if debug:
                print start
            length=10000
            if start+length>=eeg.num_datapoints:
                length=eeg.num_datapoints-start
            eeg[start:start+length,:] -= np.repeat(eeg[start:start+length,bsl].mean(axis=1).reshape(length,1),eeg.num_channels,axis=1)
            start+=length
        #for i in range(eeg.num_datapoints):
        #    eeg[i,:] -= eeg[i,bsl].mean()
        #    if debug:
        #        if i%1000==999:
        #            print i 
    except ValueError,e: #Evtl. anpassen
        pass

def parse_args(argv):
    global input, output, exclude_channels
    try:                                
        opts, args = getopt.getopt(argv, "hdi:o:e:", ["help", "input="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    #print opts
    #print args
    for opt,arg in opts:
        if opt in ["-h","--help"]:
            usage()
            exit(0)
        elif opt in ["-d"]:
            global debug
            debug = True
        elif opt in ["-i","--in"]:
            input = arg
        elif opt in ["-o","--out"]:
            output = arg
        elif opt in ["-e", "--exclude"]:
            exclude_channels = eval(arg)
    try:
        if not all(type(ch) == int for ch in exclude_channels):
                raise ValueError("exclude-channels must be all integers")
    except Exception, e:
        print "An error occured:", e
        usage()
        exit(2)
            
    
def usage():
    print """Usage details:"""
    print "=" * 30
    print """
-h, --help:
    Show usage information
-d:
    Print debugging information    
-i, --in:
    Specify input filename
-o, --out:
    Specify output filename
-e, --exclude_channels:
    Which channels to exclude from averaging
"""

if __name__== "__main__":
    main(sys.argv[1:])