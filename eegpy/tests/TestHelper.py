import numpy as np
import time
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises, raises
from eegpy.helper import fconv, nextpow2, is_power_of_2

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

def test_nextpow2():
    for x, np2 in [(0,0),(2,2),(3,4),(5,8),(33,64),(170,256)]:
        yield _check_nextpow2, x, np2

def _check_nextpow2(x,np2):
    return nextpow2(x)==np2

def test_nextpow2_float():
    for x, np2 in [(1.8,2),(19.1,32),(63.0,64),(78.4,128)]:
        yield _check_nextpow2, x, np2

@raises(ValueError)
def test_nextpow2_negative_x():
    nextpow2(-1)

def test_is_power_of_2():
    for x, is_pow_2 in [
            (0,True),
            (1,True),
            (2,True),
            (3,False),
            (4,True),
            (256, True),
            (8., True),
            (-2, True),
            (-3,False),
            ("2",False)]:
        yield lambda x,y: is_power_of_2(x) == y, x, is_pow_2

@raises(ValueError)
def test_is_power_of_2_wrong_argument_type():
    is_power_of_2("abc")
