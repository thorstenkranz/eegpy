#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR


def filterRecursively(x, func, *args, **kwargs):
    if len(x.shape)==1:
        x=func(x, *args, **kwargs)
    else:
        for i in range(x.shape[-1]):
            filterRecursively(x[...,i],func,*args,**kwargs)
            