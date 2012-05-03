import numpy as np
import time
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises
from eegpy.helper import fconv

def test_fconv_for_random_arrays():
    x = np.random.random((128))
    y = np.random.random((128))
    _assert_numpy_convolve_gives_same_results_as_fconv_for(x,y)

def test_fconv_for_smaller_random_arrays():
    x = np.random.random((60))
    y = np.random.random((60))
    _assert_numpy_convolve_gives_same_results_as_fconv_for(x,y)

def _assert_numpy_convolve_gives_same_results_as_fconv_for(x,y):
    conv_np = np.convolve(x,y)
    conv_fconv = fconv(x,y)
    assert_array_almost_equal(conv_np, conv_fconv)
        
def test_fconv_is_faster_than_convolve():
    x = np.random.random((6000))
    y = np.random.random((6000))
    t1 = time.time()
    np.convolve(x,y)
    t2 = time.time()
    fconv(x,y)
    t3 = time.time()
    assert_true((t3-t2)<(t2-t1))

def test_fconv_returns_complex_when_input_is_complex():
    x = np.random.random((10)).astype(np.complex)
    y = np.random.random((10)).astype(np.complex)
    res = fconv(x,y)
    assert_true(res.dtype==np.complex)

def test_fconv_returns_double_when_input_is_double():
    x = np.random.random((10)).astype(np.double)
    y = np.random.random((10)).astype(np.double)
    res = fconv(x,y)
    assert_true(res.dtype==np.double)


