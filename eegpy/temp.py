# -*- coding: utf-8 -*-
"""
Simple extension of mktemp to enable usage in with-statement
"""
import tempfile
import os

class UnlinkingTempfile(object):
    """Create a temporary filename. Unlink this file when leaving.
    
    Usage::
        with UnlinkingTempfile() as filename:
            with open(filename,"w") as fh:
                fh.write("Test")
    
    The created tempfile will automatically be removed again.
    """
    def __init__(self, suffix=None, prefix=tempfile.template, dir=None):
        self._suffix = "" if suffix is None else suffix
        self._prefix = prefix
        self._dir = dir
        
        self._fn = None
        
    def __enter__(self):
        self._fn = tempfile.mktemp(self._suffix, self._prefix, self._dir)
        return self._fn

    def __exit__(self, type_, value, traceback):
        try:
            if os.path.exists(self._fn):
                os.unlink(self._fn)
        finally:
            self._fn = None
    

