import os
import numpy as np
import time
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_false
from tempfile import mktemp
from eegpy.formats.f32 import F32, F32filtered

from eegpy.temp import UnlinkingTempfile


def test_write():
    for n_ch in [1,10,20]:
        for n_dp in [1,10,100,1000,10000]:
            data = np.random.random((n_dp, n_ch)).astype("d")
            yield _internal_write, data

def _internal_write(data):
    with UnlinkingTempfile(".f32") as tmp_fn:
        eeg = F32(tmp_fn, "w+", shape=data.shape)
        eeg[:,:] = data
        eeg.close()
        assert_true(os.path.getsize(tmp_fn) > 1024)
        
def test_write_and_read_back():
    for n_ch in [1,10,20]:
        for n_dp in [1,10,100,1000,10000]:
            data = np.random.random((n_dp, n_ch)).astype("d")
            yield _internal_write_and_read_back, data

def _internal_write_and_read_back(data):
    with UnlinkingTempfile(".f32") as tmp_fn:
        eeg = F32(tmp_fn, "w+", shape=data.shape)
        eeg[:,:] = data
        eeg.close()
        
        eeg2 = F32(tmp_fn)
        read_data = np.array(eeg2[:])
        assert_array_almost_equal(data, read_data)

def test_channel_names_are_automatically_created():
    with UnlinkingTempfile(".f32") as tmp_fn:
        eeg = F32(tmp_fn, "w+", shape=(1,10))
        assert_equal(10, len(eeg.channel_names))
        for ch in eeg.channel_names:
            assert_false(ch=="")
            
def test_Fs_is_persisted():
    Fs = 1234.0
    with UnlinkingTempfile(".f32") as tmp_fn:
        eeg = F32(tmp_fn, "w+", shape=(1,1), Fs=Fs)
        eeg.close()
        eeg2 = F32(tmp_fn)
        assert_equal(Fs, eeg2.Fs)

def test_filtered_reading():
    with UnlinkingTempfile(".f32") as tmp_fn:
	data = np.random.random((100,1))
        eeg = F32(tmp_fn, "w+", shape=data.shape)
        eeg[:,:] = data
        eeg.close()
        
	eeg_filtered = F32filtered(tmp_fn, lambda x: 3*x)
	read_data = eeg_filtered[:]

        assert_array_almost_equal(data, read_data/3)


