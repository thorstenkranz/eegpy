#!/usr/bin/env python
# -*- coding: utf-8 -*-

standard_bands = {"delta":[0.5,3.],
                  "theta":[3.,7.5],
                  "alpha":[7.5,12.],
                  "beta":[12.,30.],
                  "gamma":[30.,100.],
                  "low-gamma":[30.,60.],
                  "high-gamma":[60.,90.]
                  }

standard_bands_centers = {}
for k in standard_bands.keys():
    standard_bands_centers[k] = (standard_bands[k][1]+standard_bands[k][0])/2


def get_for_names(names = []):
    assert type(names) == type([]), "names-parameter must be a list of band-names"
    rv = []
    for name in names:
        try:
            rv.append(standard_bands[name])
        except KeyError, e:
            raise KeyError("Given band unknown. Must be one of %s" %str(standard_bands.keys()))
    return rv
        

def get_centers_for_names(names = []):
    rv = []
    for name in names:
        try:
            rv.append(standard_bands_centers[name])
        except KeyError, e:
            raise KeyError("Given band unknown. Must be one of %s" %str(standard_bands.keys()))
    return rv
