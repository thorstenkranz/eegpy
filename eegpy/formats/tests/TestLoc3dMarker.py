import os
import numpy as np
import time
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises, raises
from tempfile import mktemp
from eegpy.formats.loc3dmarker import Loc3dMarkers

test_data = """A,-23.58,-18.00,-20.74,5.0,0.0,0.0,1.0
B,-25.50,-18.00,25.94,5.0,0.0,0.0,1.0
C,-57.27,-18.00,-1.49,5.0,0.0,0.0,1.0
D,-9.14,-18.00,59.63,5.0,0.0,0.0,1.0"""

TMP_FN = None

def setup():
    global TMP_FN 
    TMP_FN = mktemp()
    with open(TMP_FN, "w") as fh:
	fh.write(test_data)

def teardown():
    if TMP_FN is not None and os.path.exists(TMP_FN):
    	os.unlink(TMP_FN)


def test_load_marker():
    markers = Loc3dMarkers(TMP_FN)
    assert_equal(4, markers.count)
    assert_equal("A", markers.labels[0])
    assert_equal(5.0, markers.sizes[0])

#@raises(ValueError)
#def test_nextpow2_negative_x():
#    nextpow2(-1)

