import numpy
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal
from eegpy import EventTable

class TestEvents:

    testValues = {"a":[10,123,123,199,2**5],
                  "b":[1238,12380,129387,8,8,9,8]}

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_CreateBlankEventTable (self):
        evt = EventTable()
        assert_true( evt!=0  ) 

    def test_CreateEvtFillInSomeValuesAndGetThemBack(self):
        #arange
        evt = EventTable()
        #act
        for k in self.testValues.keys():
            evt[k] = self.testValues[k]
        #assert
        for k in self.testValues.keys():
            assert len(evt[k]) == len(self.testValues[k])
            sortedValues = self.testValues[k]
            sortedValues.sort()
            assert_array_equal(sortedValues,evt[k])
        
