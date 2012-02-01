import numpy
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises
from eegpy import EventTable

class TestEvents:

    testValues = {"a":[10,123,123,199,2**5],
                  "b":[1238,12380,129387,8,8,9,8],
                  "c":[238042395724,2342349,243678387,876547,4568456,34696,4863],
                 }

    def setUp(self):
        self.evt = EventTable()

    def tearDown(self):
        pass

    def test_CreateBlankEventTable (self):
        assert_true( self.evt!=None ) 

    def test_CreateEvtFillInSomeValuesAndGetThemBack(self):
        #arange
        #act
        for k in self.testValues.keys():
            self.evt[k] = self.testValues[k]
        #assert
        for k in self.testValues.keys():
            assert len(self.evt[k]) == len(self.testValues[k])
            sortedValues = self.testValues[k]
            sortedValues.sort()
            assert_array_equal(sortedValues,self.evt[k])
        
    def test_AddTrigger(self):
        self.evt.add_trigger(self.testValues.keys()[0],self.testValues.values()[0][0])
        
    def test_AddTriggerRaisesExceptionForListAsTime(self):
        with assert_raises(ValueError) as cm:
            self.evt.add_trigger(self.testValues.keys()[0],self.testValues.values()[0])

    def test_AddTriggerType(self):
        self.evt.add_trigger_type(self.testValues.keys()[0],self.testValues.values()[0])

    def test_AddTriggerTypeRaisesValueErrorForIntAsList(self):
        with assert_raises(TypeError) as cm:
            self.evt.add_trigger_type(self.testValues.keys()[0],self.testValues.values()[0][0])

    def test_Keys(self):
        for k in self.testValues.keys():
            self.evt[k] = self.testValues[k]
        tv_keys = self.testValues.keys()
        evt_keys = self.evt.keys()
        tv_keys.sort()
        evt_keys.sort()
        assert_array_equal(tv_keys,evt_keys)
        
