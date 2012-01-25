#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mvpa
from mvpa.mappers.base import Mapper, ChainMapper

import numpy as np
from scipy.signal import detrend, hilbert
import pylab as p

from eegpy.helper import is_power_of_2
from eegpy.analysis.wavelet import wavepower_lin, wt_analyze
from eegpy.stats.cluster import ClusterSearch1d, ClusterSearchDiscreteWavelet
from eegpy.filter.freqfilt import filtfilt, butter, _filtfilt
from eegpy.filter import pcafilt

"""Define some custom Mappers"""

class Cluster1dPermutationMapper(Mapper):
    """Mapper that enables the search for significant clusters 
    according to Marvis, Oostenveld (2007), 
    "Nonparametric statistical testing of EEG and MEG data".
    """
    def __init__(self,num_surrogates=10,blocksize=None,*args,**kwargs):
        Mapper.__init__(self,*args,**kwargs)
        self._num_surrogates = num_surrogates
        self._clusters = []  # Will become list of tuples (channel, band, (start,end), probability)
        self._is_trained = False
        #Shape for each sample of the data this mapper was trained on
        self._sample_shape = None
        self._blocksize = blocksize

    def _make_samples_4d(self, d_tmp):
        """We expect 4d data, but sometimes might also hand in less dimensions"""
        dts = d_tmp.shape
        if len(dts)==2:
            d_tmp.shape = (dts[0],dts[1],1,1)
        elif len(dts)==3:
            d_tmp.shape = (dts[0],dts[1],dts[2],1)
        elif len(dts)==4:
            pass
        else:
            raise ValueError("Samples must be one of 1d, 2d or 3d.")
        return d_tmp

    def train(self,dataset):
        """Do the actual search for significant cluster. 
        Enables mapper to do the forward mapping.

        """
        data = []
        self._clusters=[]
        self._sample_shape=dataset.samples_original.shape[1:]
        for i_l, lab in enumerate(dataset.uniquelabels):
            d_tmp = dataset.samples_original[dataset.labels==lab]
            self._make_samples_4d(d_tmp)
            data.append(d_tmp)
        self.search_clusters(data)
        self._is_trained = True

    def forward(self,data):
        if not data.shape[1:] == self._sample_shape:
            print data.shape[1:], self._sample_shape
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        self._make_samples_4d(data)
        #Shorthands
        bs = self._blocksize
        #Make outsamples
        if bs==None:
            num_features = len(self.clusters)
        else:
            num_features = np.sum([np.ceil(float(cl[2][1]-cl[2][0])/bs) for cl in self.clusters])
        rv = np.zeros((data.shape[0],num_features))
        for i_d in range(data.shape[0]):
            i_f = 0 #Index of feature
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (channel,band,(start,end),prob) = cluster
                if bs==None:
                    rv[i_d,i_c] = data[i_d,start:end,channel,band].mean()
                else:
                    for avgstart in np.arange(start,end,bs):
                        #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                        rv[i_d,i_f] = data[i_d,avgstart:min(avgstart+bs,end),channel,band].mean()
                        i_f+=1
        return rv

    def untrain(self):
        self._is_trained=False
        self._clusters=[]
        self._sample_shape = None

    @property
    def is_trained(self):
        return self._is_trained

    @property
    def clusters(self):
        return self._clusters

    def search_clusters(self, data):
        """This function performs the search for significant clusters. 
           Uses a class from eegpy.stats.cluster
        """
        plotid=0
        print data[0].shape
        for i_ch in range(data[0].shape[2]):
            for i_b in range(data[0].shape[3]):
                #print i_ch,i_b
                cs_data = [d[:,:,i_ch,i_b].T for d in data]
                #print len(cs_data), [csd.shape for csd in cs_data]
                fs,sct,scp,ct,cp = ClusterSearch1d(cs_data,num_surrogates=self._num_surrogates).search()
                #print "CPM: fs.shape=",fs.shape
                #print ct, cp
                for i_c in range(len(sct)):
                    self._clusters.append( (i_ch,i_b,sct[i_c],scp[i_c]) ) 
                #Do plot if significant cluster found
                #if len(sct)>0:
                #    #print "CPM: fs.shape=",fs.shape, fs.dtype
                #    #print "CPM: fs[499:510]", fs[498:510]
                #    p.plot(np.array(fs))
                #    p.title("Channel %i, band %i"%(i_ch,i_b))
                #    for cluster in sct:
                #        p.axvspan(cluster[0],cluster[1],color="y",alpha=0.2)
                #    p.savefig("/tmp/fbpm%03d.png"%plotid)
                #    plotid+=1
                #    p.clf()

class FrequencyBandPowerMapper(Mapper):
    """For all frequeny bands that are given to the constructor, filter the data,
    calculate hiblert transformation and calculate the power. Maybe do a baseline correction."""
    
    def __init__(self, bands, border=2, baseline_interval=None, data_interval=None, Fs=1000.0, *args, **kwargs):
        """Initialize the FrequencyBandPowerMapper
   
          :Parameters:
            bands: list of 2-tuples
              list of (fl,fh), edge-frequencies for filtering
            border: int
              order of the filter to be used
            baseline_interval: 2-tuple of ints
              indices of the interval to be used for baseline correction
            data_interval: 2-tuple of ints
              indices of the interval to be used for cluster search
            args, kwargs: 
              additional arguments passed to parent class
        """
        Mapper.__init__(self,*args,**kwargs)
        self._bands = bands
        self._border = border
        self._baseline_interval = baseline_interval
        self._data_interval = data_interval
        self._Fs = Fs
        #TODO: Perform some checks, for valid intervals etc.
        #precalculate filter parameters (butter)
        self._bs = list(np.zeros((len(bands))))
        self._as = list(np.zeros((len(bands))))
        for i_b in range(len(bands)):
            [self._bs[i_b],self._as[i_b]]=butter(border,[bands[i_b][0]/(Fs/2),bands[i_b][1]/(Fs/2)], btype="band")
        

    def train(self,dataset):
        """Do the actual search for significant cluster. 
        Enables mapper to do the forward mapping.

        """
        #data = []
        #self._clusters=[]
        self._sample_shape=dataset.samples_original.shape[1:]
        #TODO: assert intervals in sample shape

        ##print "# ", dataset.samples.shape
        #for i_l, lab in enumerate(dataset.uniquelabels):
        #    d_tmp = dataset.samples_original[dataset.labels==lab]
        #    data.append(self.calculate_power(d_tmp))
        self._is_trained = True

    def calculate_power(self, data):
        rv = np.zeros((data.shape[0],data.shape[1],data.shape[2],len(self._bands)))
        #Shorthands
        bli = self._baseline_interval
        di = self._data_interval
        for i_d in range(data.shape[0]): #Loop over samples
            for i_b, band in enumerate(self._bands): #Loop over frequency bands
                #print "FBPM.c_p: data[i_d,:,:].shape", data[i_d,:,:].shape
                for  i_ch in range(data.shape[2]):
                    rv[i_d,:,i_ch,i_b] = abs(hilbert(_filtfilt(self._bs[i_b],self._as[i_b],data[i_d,:,i_ch])))
                    #rv[i_d,:,i_ch,i_b] = abs(data[i_d,:,i_ch])
                #p.plot(np.arange(rv.shape[1]),rv[i_d,:,:,i_b])
                #assert 1==0
                #Baseline correction
                if bli!=None:
                    rv[i_d,:,:,i_b] /= rv[i_d,bli[0]:bli[1],:,i_b].mean(axis=0).reshape(1,-1).repeat(rv.shape[1],axis=0)
        if di == None:
            return rv
        else:
            #Return only the selected data interval
            return rv[:,di[0]:di[1],...]

    def forward(self,data):
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        rv = self.calculate_power(data)
        return rv

    def __repr__(self):
        """String summary over the object
        """
        return "FrequencyBandPowerMapper"


class PCAClusterMapper(Cluster1dPermutationMapper):
    """Calculates global mean over trials, then performs PCA on this.
    Projection is then performed for each trial. 
    Then searches for clusters in the PCs, only in a certain fraction 
    of the components."""
    def __init__(self, baseline_interval=None, data_interval=None, use_fraction=0.2, *args, **kwargs):
        """Initialize the PCAClusterMapper
   
          :Parameters:
            baseline_interval: 2-tuple of ints
              indices of the interval to be used for baseline correction
            data_interval: 2-tuple of ints
              indices of the interval to be used for cluster search
            use_fraction: float
              only the first int(round(n_components*use_fraction)) components are used,
              0<use_fraction<=1
            args, kwargs: 
              additional arguments passed to parent class
        """
        Cluster1dPermutationMapper.__init__(self,*args,**kwargs)
        self._baseline_interval = baseline_interval
        self._data_interval = data_interval
        self._use_fraction = use_fraction
        self._pcanalyzer = None
        #TODO: Perform some checks, for valid intervals etc.

    def train(self,dataset):
        """Prepare mapper to behave as desired. Do PCA, then search clusters.

        :Parameters:
          dataset: MVPA Dataset
            The dataset the mapper is trained for
        """
        data = []
        self._clusters=[] 
        self._sample_shape=dataset.samples_original.shape[1:]
        #Global mean over trials
        avg_data = dataset.samples_original.mean(axis=0)
        bli = self._baseline_interval #Baseline-correction
        avg_data -= avg_data[bli[0]:bli[1],:].mean(axis=0).reshape(1,-1).repeat(avg_data.shape[0],axis=0) 
        #PCA on this
        self._pcanalyzer = pcafilt.PCAnalyzer(avg_data) #Does training of node, but no execution
        #print "# ", dataset.samples.shape
        for i_l, lab in enumerate(dataset.uniquelabels):
            d_tmp = dataset.samples_original[dataset.labels==lab]
            data.append(self.unmix(d_tmp))
        self.search_clusters(data)
        self._is_trained = True

    @property
    def num_components(self):
        return int(round(self._use_fraction*self._sample_shape[1]))

    def unmix(self,data):
        """Perform PCA transform."""
        if self._pcanalyzer==None:
            raise RuntimeError("PCAnalyzer not trained yet.")
        #Baseline-correction
        bli = self._baseline_interval
        data_bl = np.zeros_like(data)
        for i_d in range(data.shape[0]):
            for i_ch in range(data.shape[2]):
                data_bl[i_d,:,i_ch] = data[i_d,:,i_ch]/data[i_d,bli[0]:bli[1],i_ch].mean()
        #PCA
        rv = np.zeros_like(data)
        for i_d in range(data.shape[0]):
            rv[i_d,:,:] = self._pcanalyzer.unmix(data[i_d,:,:]) #Do projection as learned in train()
        #return only data interval
        di = self._data_interval
        return rv[:,di[0]:di[1],:]

    def search_clusters(self, data):
        """This function performs the search for significant clusters. 
           Uses a class from eegpy.stats.cluster
        """
        plotid=0
        for i_cp in range(self.num_components):
            print "Searching in component %i of %i"%(i_cp+1,self.num_components)
            cs_data = [d[:,:,i_cp].T for d in data]
            #print len(cs_data), [csd.shape for csd in cs_data]
            fs,sct,scp,ct,cp = ClusterSearch1d(cs_data,num_surrogates=self._num_surrogates).search()
            #print "CPM: fs.shape=",fs.shape
            #print ct, cp
            for i_c in range(len(sct)):
                self._clusters.append( (i_cp,sct[i_c],scp[i_c]) ) 
            #Do plot if significant cluster found
            if len(sct)>0:
                #print "CPM: fs.shape=",fs.shape, fs.dtype
                #print "CPM: fs[499:510]", fs[498:510]
                p.plot(np.array(fs))
                p.title("Component %i"%(i_cp))
                for cluster in sct:
                    p.axvspan(cluster[0],cluster[1],color="y",alpha=0.2)
                p.savefig("/tmp/pcacm%03d.png"%plotid)
                plotid+=1
                p.clf()

    def forward(self,data):
        """Do the forward mapping"""
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        #Make outsamples, shape is nsamples x nclusters
        rv = np.zeros((data.shape[0],len(self.clusters)))
        components = self.unmix(data)
        for i_d in range(data.shape[0]):
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (comp,(start,end),prob) = cluster
                #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                rv[i_d,i_c] = components[i_d,start:end,comp].mean()
        return rv

    def __repr__(self):
        return "PCAClusterMapper"

class PCAClusterMapperBlockAverage(PCAClusterMapper):
    """Same behaviour as parent class, just does block-averaging within clusters
    instead of averaging over the whole clusters.
    """
    def __init__(self,blocksize=10,*args,**kwargs):
        PCAClusterMapper.__init__(self,*args,**kwargs)
        self._blocksize = blocksize

    def forward(self,data):
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        #Shorthands
        bs = self._blocksize
        #Make outsamples
        num_features = np.sum([np.ceil(float(cl[1][1]-cl[1][0])/bs) for cl in self.clusters])
        rv = np.zeros((data.shape[0],num_features))
        components = self.unmix(data)
        for i_d in range(data.shape[0]):
            i_f = 0 #Index of feature
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (comp,(start,end),prob) = cluster
                for avgstart in np.arange(start,end,bs):
                    #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                    rv[i_d,i_f] = components[i_d,avgstart:min(avgstart+bs,end),comp].mean()
                    i_f+=1
        return rv

    def __repr__(self):
        return "PCAClusterMapperBlockAverage"

class WaveletClusterMapper(Cluster1dPermutationMapper):
    """Calculates WaveletTransformation and power for each trial (discrete). 
    Then searches for clusters in the Wavelet powers.

    Normalization is performed as follows: data are expected to be 
    1.5 x as long as the desired wavelet-data (e.g. 1536 for 1024 wavelet-coefficients),
    reaching from -0.5 to 1.0 of time analyzed. We then do two Wavelet-Transformations, 
    one from -0.5 to 0.5 and one from 0 to 1. The power from 0 to 1 is then divided
    by the power immediately before t=0, i.e. in the middle of the -0.5 to 0.5-interval.
    """
    def __init__(self, wavelet="db4", normalize=True, num_wts_to_take=512, *args, **kwargs):
        """Initialize the WaveletClusterMapper
   
          :Parameters:
            normalize: bool
              Flag whether normalization should be applied. 
            args, kwargs: 
              additional arguments passed to parent class
        """
        Cluster1dPermutationMapper.__init__(self,*args,**kwargs)
        self._wavelet = wavelet
        self._normalize = normalize
        self._num_wts_to_take = num_wts_to_take
        #TODO: Perform some checks, for valid intervals etc.

    def train(self,dataset):
        """Prepare mapper to behave as desired. Calculate wavelet power, then search clusters.

        :Parameters:
          dataset: MVPA Dataset
            The dataset the mapper is trained for. Data must range from -0.5 
            to 1.0 (in units of time-range to be analyzed, 
            e.g. from -512 to 1024 for 1024 datapoints to be used.)
        """

        data = []
        self._sample_shape=dataset.samples_original.shape[1:]
        self._clusters=[] 
        self.search_clusters(data)
        self._is_trained = True

    def calculate_wt_power(self,data):
        """For given EEG-data, calculate baseline-corrected Wavelet-power"""
        #Checks
        data_len = dataset.samples_original.shape[1]
        if not is_power_of_2(data_len*2/3):
            raise ValueError("length of data must be 1.5 x a power of 2")
        # Shorthands 
        nwtt = self._num_wts_to_take
        # Calculate
        rv = np.zeros((data.shape),"d")
        for i_s in range(data.shape[0]):
            for i_ch in range(data.shape[2]):
                #TODO: Noch nicht implementiert
                rv[i_s,:,i_ch] = wt_power()
    
    def search_clusters(self, data):
        """This function performs the search for significant clusters. 
           Uses a class from eegpy.stats.cluster
        """
        plotid=0
        for i_ch in range(data[0].shape[2]):
            for i_b in range(data[0].shape[3]):
                #print i_ch,i_b
                cs_data = [d[:,:,i_ch,i_b].T for d in data]
                #print len(cs_data), [csd.shape for csd in cs_data]
                fs,sct,scp,ct,cp = ClusterSearch1d(cs_data,num_surrogates=self._num_surrogates).search()
                #print "CPM: fs.shape=",fs.shape
                #print ct, cp
                for i_c in range(len(sct)):
                    self._clusters.append( (i_ch,i_b,sct[i_c],scp[i_c]) ) 
                #Do plot if significant cluster found
                if len(sct)>0:
                    #print "CPM: fs.shape=",fs.shape, fs.dtype
                    #print "CPM: fs[499:510]", fs[498:510]
                    p.plot(np.array(fs))
                    p.title("Channel %i, band %i"%(i_ch,i_b))
                    for cluster in sct:
                        p.axvspan(cluster[0],cluster[1],color="y",alpha=0.2)
                    p.savefig("/tmp/fbpm%03d.png"%plotid)
                    plotid+=1
                    p.clf()

    def transform(self, data):
        """Perform transformation and normalization (if desired)"""
        data_len = data.shape[1]
        rv = np.zeros((data.shape[0],data.shape[1]*2/3,data.shape[2]),"d")
        for i_tr in range(rv.shape[0]):
            for i_ch in range(rv.shape[2]):
                rv[i_tr,:,i_ch] = wavepower_lin(data[i_tr,data_len/3:,i_ch],self._wavelet,normalise=False)
        if self._normalize:
            normdata = np.zeros((data.shape[0],data.shape[1]*2/3,data.shape[2]),"d")
            for i_tr in range(normdata.shape[0]):
                for i_ch in range(normdata.shape[2]):
                    normdata[i_tr,:,i_ch] = wavepower_lin(data[i_tr,:data_len*2/3,i_ch],self._wavelet,normalise=False)
            #Do normalization

    def forward(self,data):
        """Do the forward mapping"""
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        #Make outsamples, shape is nsamples x nclusters
        rv = np.zeros((data.shape[0],len(self.clusters)))
        components = self.unmix(data)
        for i_d in range(data.shape[0]):
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (comp,(start,end),prob) = cluster
                #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                rv[i_d,i_c] = components[i_d,start:end,comp].mean()
        return rv

    def __repr__(self):
        return "PCAClusterMapper"

class MorletPowerMapper(Cluster1dPermutationMapper):
    """For all frequeny bands that are given to the constructor, filter the data,
    calculate hiblert transformation and calculate the power. Maybe do a baseline correction."""
    
    def __init__(self, freqs=None, baseline_interval=None, data_interval=None, Fs=1000.0, *args, **kwargs):
        """Initialize the FrequencyBandPowerMapper
   
          :Parameters:
            freqs: list of 2-tuples
              list of frequencies for Morlet-Wavelets
            baseline_interval: 2-tuple of ints
              indices of the interval to be used for baseline correction
            data_interval: 2-tuple of ints
              indices of the interval to be used for cluster search
            Fs: float
              sampling frequency in Hz
            args, kwargs: 
              additional arguments passed to parent class
        """
        Cluster1dPermutationMapper.__init__(self,*args,**kwargs)
        if freqs == None:
            self._freqs = np.logspace(0,np.log2(100),10,True,2)
        else:
            self._freqs = np.array(freqs)
        self._baseline_interval = baseline_interval
        self._data_interval = data_interval
        self._Fs = Fs
        #TODO: Perform some checks, for valid intervals etc.
        

    def train(self,dataset):
        """Do the actual search for significant cluster. 
        Enables mapper to do the forward mapping.

        """
        data = []
        self._clusters=[]
        self._sample_shape=dataset.samples_original.shape[1:]
        #print "# ", dataset.samples.shape
        for i_l, lab in enumerate(dataset.uniquelabels):
            d_tmp = dataset.samples_original[dataset.labels==lab]
            data.append(self.calculate_power(d_tmp))
        print "Searching cluster...", 
        self.search_clusters(data)
        print "finished"
        self._is_trained = True

    def calculate_power(self, data):
        rv = np.zeros((data.shape[0],data.shape[1],data.shape[2],len(self._freqs)))
        #Shorthands
        bli = self._baseline_interval
        di = self._data_interval
        #iterate over channels and trials
        for i_d in range(data.shape[0]): #Loop over samples (e.g. trials, not timepoints!)
            #print "MPM.c_p: i_d, data[i_d,:,:].shape", i_d, data[i_d,:,:].shape
            for i_ch in range(data.shape[2]):
                rv[i_d,:,i_ch,:] = abs(wt_analyze(data[i_d,:,i_ch],self._freqs))**2
            #p.plot(np.arange(rv.shape[1]),rv[i_d,:,:,i_b])
            #assert 1==0
            #Baseline correction
            for i_f in range(len(self._freqs)):
                if bli!=None:
                    rv[i_d,:,:,i_f] /= rv[i_d,bli[0]:bli[1],:,i_f].mean(axis=0).reshape(1,-1).repeat(rv.shape[1],axis=0)
        if di == None:
            return rv
        else:
            #Return only the selected data interval
            return rv[:,di[0]:di[1],...]

    def forward(self,data):
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        #Make outsamples, shape is nsamples x nclusters
        rv = np.zeros((data.shape[0],len(self.clusters)))
        power_ar = self.calculate_power(data)
        for i_d in range(data.shape[0]):
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (channel,band,(start,end),prob) = cluster
                #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                rv[i_d,i_c] = power_ar[i_d,start:end,channel,band].mean()
        return rv

    def __repr__(self):
        """String summary over the object
        """
        return "FrequencyBandPowerMapper"

class MorletPowerMapperBlockAverage(MorletPowerMapper):
    """Same behaviour as parent class, just does block-averaging within clusters
    instead of averaging over the whole clusters.
    """
    def __init__(self,blocksize=10,*args,**kwargs):
        MorletPowerMapper.__init__(self,*args,**kwargs)
        self._blocksize = blocksize

    def forward(self,data):
        if not data.shape[1:] == self._sample_shape:
            raise ValueError("Forward mapping failed. Each sample of the data given must have shape %s"%str(self._sample_shape))
        #Shorthands
        bs = self._blocksize
        #Make outsamples
        num_features = np.sum([np.ceil(float(cl[2][1]-cl[2][0])/bs) for cl in self.clusters])
        rv = np.zeros((data.shape[0],num_features))
        power_ar = self.calculate_power(data)
        for i_d in range(data.shape[0]):
            i_f = 0 #Index of feature
            for i_c in range(len(self.clusters)):
                cluster = self.clusters[i_c]
                (channel,band,(start,end),prob) = cluster
                for avgstart in np.arange(start,end,bs):
                    #print "channel, band, start, end, power_ar.shape",channel, band, start, end, power_ar.shape
                    rv[i_d,i_f] = power_ar[i_d,avgstart:min(avgstart+bs,end),channel,band].mean()
                    i_f+=1
        return rv

    def __repr__(self):
        """String summary over the object
        """
        return "FrequencyBandPowerMapperBlockAverage"

if __name__ == "__main__":
    from mvpa.suite import MaskedDataset
    noiselevel = 20 
    condition1 = np.random.random((50,500,10))*noiselevel
    normfactor = np.hanning(20).sum()
    for i in range(50):
        for j in range(10):
            condition1[i,:,j] = np.convolve(condition1[i,:,j],np.hanning(20),mode="same")/normfactor
    condition2 = np.random.random((43,500,10))*noiselevel
    for i in range(43):
        for j in range(10):
            condition2[i,:,j] = np.convolve(condition2[i,:,j],np.hanning(20),mode="same")/normfactor
    pseudoekp = np.hanning(200).reshape(1,-1).repeat(50,axis=0)
    pseudoekp = pseudoekp.reshape(50,200,1).repeat(10,axis=-1)
    condition1[:,200:400]+=pseudoekp[:]
    condition2[:,200:400]-=pseudoekp[:43]
    dataset = MaskedDataset(samples = np.concatenate([condition1[:,20:-20,:],condition2[:,20:-20,:]],axis=0), labels=[0]*50+[1]*43)
    cm = Cluster1dPermutationMapper(num_surrogates=20,blocksize=None)
    cm.train(dataset)
    print cm.clusters
    cl_data = cm.forward(dataset.samples_original)
    print cl_data.shape
    cm = Cluster1dPermutationMapper(num_surrogates=20,blocksize=10)
    cm.train(dataset)
    print cm.clusters
    cl_data = cm.forward(dataset.samples_original)
    print cl_data.shape
    cm = Cluster1dPermutationMapper(num_surrogates=20,blocksize=4)
    cm.train(dataset)
    print cm.clusters
    cl_data = cm.forward(dataset.samples_original)
    print cl_data.shape
    cm = Cluster1dPermutationMapper(num_surrogates=20,blocksize=1)
    cm.train(dataset)
    print cm.clusters
    cl_data = cm.forward(dataset.samples_original)
    print cl_data.shape
    #mpm = MorletPowerMapper(np.arange(1,20,10),baseline_interval=[0,200])
    #mpm.train(dataset)
    #import subprocess
    #i=0
    #while i<10000000000:
    #    if i % 100 == 0:
    #        print i
    #        subprocess.call(["free","-m"])
    #    tmp = np.random.random((50,460,10))
    #    tmp2 = mpm.forward(tmp)
    #    i+=1
    #assert 1==0
    #fs, signif_cluster_times, signif_cluster_probs, cluster_times, cluster_probs = cs.search()

    ##Plotting for a better understanding
    import pylab as p
    p.subplot(211)
    p.plot(condition1.mean(axis=1),label="Condition 1")
    p.plot(condition2.mean(axis=1),label="Condition 2")
    p.ylabel("signal [a.u.]")
    p.subplot(212)
    #for i_c in range(len(cluster_times)):
    #    start, end = cluster_times[i_c]
    #    p.axvspan(start,end,color="b",alpha=0.3)
    #    p.text(start,2,"%.2f"%cluster_probs[i_c],fontsize=8)
    signif_cluster_times = [cl[2] for cl in cm.clusters]
    for i_c in range(len(signif_cluster_times)):
        start, end = signif_cluster_times[i_c]
        p.axvspan(start,end,0.95,1.0,color="b",alpha=1.0)
    #p.plot(fs)
    p.xlabel("timepoints")
    p.ylabel("f-values")
    p.show()
