import numpy as np
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises, raises

from eegpy.stats.cluster import ClusterSearch1d
from eegpy.filter.smoothing import smooth

        
#Helper-methods
def _make_data_with_delta(n):
    ar1 = np.zeros((n))
    ar1[n/2] = 1
    return ar1

def _check_data_decay_from_center(data):
    def _assert_range_decreases(data,indices):
        last_val=10**30
        for i in indices:
            print i, data[i], last_val
            assert data[i]<=last_val 
            last_val=data[i]
            if abs(last_val)<10**-10:
                return True
    n=len(data)
    _assert_range_decreases(data,range(n/2,n,1))
    _assert_range_decreases(data,range(n/2,-1,-1))
    return True


def test_DeltaN1000W100():
    #arange
    data = _make_data_with_delta(1000)
    #act
    data_smooth = smooth(data,100)
    #assert
    assert data.shape == data_smooth.shape
    assert _check_data_decay_from_center(data_smooth)

def test_DeltaN100W10():
    #arange
    data = _make_data_with_delta(100)
    #act
    data_smooth = smooth(data,10)
    #assert
    assert data.shape == data_smooth.shape
    assert _check_data_decay_from_center(data_smooth)

def test_DeltaN10W100():
    #arange
    data = _make_data_with_delta(100)
    #act
    data_smooth = smooth(data,10)
    #assert
    assert data.shape == data_smooth.shape
    assert _check_data_decay_from_center(data_smooth)

def test_DeltaN100W10FlatWindow():
    n=100
    w=10
    #arange
    data = _make_data_with_delta(n)
    #act
    data_smooth = smooth(data,w,"flat")
    #assert
    assert data.shape == data_smooth.shape
    assert_array_almost_equal(data_smooth[n/2-w/2+1:n/2-w/2+1+w],np.ones((w))/w)


