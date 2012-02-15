import numpy as np
from numpy.testing import (assert_array_almost_equal,
                                   assert_array_equal)
from nose.tools import assert_true, assert_equal, assert_raises

from eegpy.stats.cluster import ClusterSearch1d
from eegpy.filter.smoothing import smooth

        
#Helper-methods
def make_data_without_cluster():
    ar1 = np.random.random((100,15))
    ar2 = np.random.random((100,15))
    return [ar1,ar2]

def make_data_with_one_cluster():
    ar1,ar2 = make_data_without_cluster()
    ar1[40:60,:] += np.hanning(20).reshape(-1,1).repeat(15,axis=1)
    return [ar1,ar2]

#Create Testdata
data_without_cluster = make_data_without_cluster()
data_with_one_cluster = make_data_with_one_cluster()

#Actual test fixture
class TestEvents:

    def setUp(self):
        self.data_with_one_cluster = [ar.copy() for ar in data_with_one_cluster] 
        self.data_without_cluster = [ar.copy() for ar in data_without_cluster] 

    def tearDown(self):
        pass

    #def test_CreateBlankEventTable (self):
    #    assert_true( self.evt!=None ) 

    def test_FindOneCluster(self):
        #arange
        cl1d = ClusterSearch1d(self.data_with_one_cluster,num_surrogates=100)
        #act
        results = cl1d.search()        
        cluster = results[1][0]
        #assert
        assert len(results[1]) == 1
        assert cluster[0]>40
        assert cluster[1]<60

    def test_FindNoCluster(self):
        #arange
        cl1d = ClusterSearch1d(self.data_without_cluster,num_surrogates=100)
        #act
        results = cl1d.search()        
        #assert
        assert len(results[1]) == 0
        assert results[4].min()>0.05
        
