#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import sys
import time
from Queue import Queue

import mvpa
from mvpa.datasets.splitters import Splitter

import numpy as N
np = N
from scipy.stats import f_oneway, percentileofscore

from eegpy.helper import fprobi

"""Clustering algorithm like used in field trip, see Maris/Oostenveld"""

class ClusterSearch1d:
    """Search for clusters in 1d-data"""

    array_list_error_message = "array_list must be a list of 2d-arrays: dim1 is time, dim2 repetitions"
    num_surrogates_error_message = "num_surrogates must be a positive integer"

    def __init__(self,array_list,stat_fun=None,threshold=1.67,num_surrogates=1000):
        """Initialization

           :Parameters:
             array_list: list 
               List of 2d-arrays containing the data, dim 1: timepoints, dim 2: elements of groups
             stat_fun : function
               function called to calculate statistics, must accept 1d-arrays as arguments (default: scipy.stats.f_oneway)
             threshold: float
               Threshold to use for deciding on statistical significance (uncorrected)
             num_surrogates: int
               number of permutation-surrogates to use for testing
        """
        #Some checks
        try:
            for ar in array_list:
                if not len(ar.shape) == 2:
                    raise ValueError(self.array_list_error_message)
        except Exception, e:
            raise ValueError(self.array_list_error_message)
        if not num_surrogates>0 or type(num_surrogates)!=int:
            raise ValueError(self.num_surrogates_error_message)


        self._al = array_list
        if stat_fun==None:
            self._sf = self.f_oneway
        else:
            self._sf = stat_fun
        self._threshold = threshold
        self._ns = num_surrogates
        #self._not_anova = not_anova

    #make read-only properties
    @property
    def array_list(self):
        return self._al
    
    @property
    def stat_fun(self):
        return self._sf
    
    @property
    def threshold(self):
        return self._threshold
    
    @property
    def num_surrogates(self):
        return self._ns

    #More interesting methods...
    def f_oneway(self,*args):
        """Call scipy.stats.f_oneway, but return only f-value"""
        return f_oneway(*args)[0]

    #THE method...
    def search(self):
        """For a list of 2d-arrays of data, e.g. power values, calculate some
        statistics for each timepoint (dim 1) over groups.  Do a cluster
        analysis with permutation test like in Maris, Oostenveld (2007)
        for calculating corrected p-values
        """
        #Create Shorthands
        al = self._al
        sf = self._sf

        #print len(al), [ar.shape for ar in al]
        ns_trs = [ar.shape[1] for ar in al] # Number of trials for each group
        #if not self._not_anova:
        #    crit_f = fprobi(len(al)-1,N.sum(ns_trs)-1,self._threshold) #Critical F-value
        #else:
        crit_f = self._threshold
        l=al[0].shape[0]
        #print "CS1d: l=",l
        #Calculate Anova (or other stat_fun)
        fs = N.zeros((l),"d")
        for i in range(l):
            anova_ars = [ar[i,:] for ar in al]
            fs[i] = sf(*anova_ars)
        clusters = self.find_clusters(fs,crit_f,"greater")
        if len(clusters)>0:
            cluster_stats = [self.calc_cluster_stats(c,fs) for c in clusters]
            cluster_ps = N.ones((len(clusters)),"d")
            cluster_stats_hist = N.zeros((self._ns)) #For making histogram (not visually) and finding percentile
            ar_shuffle = N.concatenate(al,axis=1)
            for i_s in range(self._ns):
                #Mache Liste mit Indices fuer alle Trials, permutiere, zerlege in Einzellisten der Laengen ns_trs
                indices_lists = N.split(N.random.permutation(sum(ns_trs)),N.cumsum(ns_trs)[:-1])
                #print ar_shuffle.shape, ar_shuffle
                ar_shuffle_list = [ar_shuffle[:,indices] for indices in indices_lists]
                #print "ar_shuffle_list shapes", [ar.shape for ar in ar_shuffle_list]
                fs_surr = N.zeros((l))
                for i in range(l):
                    anova_ars_perm = [ar[i,:] for ar in ar_shuffle_list]
                    fs_surr[i] = sf(*anova_ars_perm)
                clusters_perm = self.find_clusters(fs_surr,crit_f,"greater")
                #print "clusters_perm", clusters_perm
                if len(clusters_perm)>0:
                    cluster_stats_perm = [self.calc_cluster_stats(c,fs_surr) for c in clusters_perm]
                    cluster_stats_hist[i_s] = max(cluster_stats_perm)
                else:
                    cluster_stats_hist[i_s] = 0
            cluster_ps[:] = [percentileofscore(cluster_stats_hist,cluster_stats[i_cl]) for i_cl in range(len(clusters))]
            cluster_ps[:] = (100.0 - cluster_ps[:]) / 100.0 # From percent to fraction
            #print "CS1d: fs.shape[0]=",fs.shape[0]
            #Set cluster-ps in clusters
            self._set_cluster_ps(clusters,cluster_ps)
            return fs, N.array(clusters)[cluster_ps<=0.05], cluster_ps[cluster_ps<=0.05], N.array(clusters), cluster_ps
        else:
            return fs,N.array([]),N.array([]),N.array([]),N.array([])
    
    def find_clusters(self,ar,thres,cmp_type="greater"):
        """For a given 1d-array (test statistic), find all clusters which
        are above/below a certain threshold. Returns a list of 2-tuples.
        """
        #clusters =  []
        if not cmp_type in ["lower","greater","abs_greater"]:
            raise ValueError("cmp_type must be in [\"lower\",\"greater\",\"abs_greater\"]")
        ar = N.concatenate([N.array([thres]),ar,N.array([thres])])
        if cmp_type=="lower":
            ar_in = (ar<thres).astype(N.int)
        elif cmp_type=="greater":
            ar_in = (ar>thres).astype(N.int)
        else: #cmp_type=="abs_greater":
            ar_in = (abs(ar)>thres).astype(N.int)
        ar_switch = N.diff(ar_in)
        inpoints = N.arange(ar.shape[0])[ar_switch>0]
        outpoints = N.arange(ar.shape[0])[ar_switch<0]
        #print inpoints, outpoints
        in_out = N.concatenate([inpoints.reshape(-1,1),outpoints.reshape(-1,1)],axis=1)
        clusters = [(c[0],c[1]) for c in in_out]
        return clusters

    def calc_cluster_stats(self,cluster,fs):
        return N.sum(fs[cluster[0]:cluster[1]])

    def _set_cluster_ps(self,cs,ps):
        pass


class ClusterBandPower:
    """A 1d cluster in band-based frequency analysis"""
    def __init__(self,ch,band,se,prob,descr = None):
        """Initialization

           :Parameters:
             ch : int 
               channel-number
             band : tuple
               2-tuple containing the edge-frequencies of the frequency-band
             se : tuple
               start/end of cluster, in samples, relative to trigger
             prob : float
               Propability for accepting the null-hypothesis, as calculated via ???cluster???
        """
        if not ch > 0 and not type(ch)==int:
            raise ValueError("ch must be a positive int")
        self._ch = ch
        #TODO: Further checks
        self._band = band
        self._se = se
        self._prob = prob
        self._descr = descr

    #make read-only properties
    @property
    def ch(self):
        return self._ch

    @property
    def band(self):
        return self._band

    @property
    def se(self):
        return self._se

    @property
    def prob(self):
        return self._prob

    @property
    def descr(self):
        return self._descr

    @property
    def size(self):
        return self._se[1]-self._se[0]

    def __repr__(self):
        """String representation of the cluster-object
        """
        return "ClusterBandPower(ch=%i, band=%s, se=%s, prob=%.3f, descr=%s)" % (self._ch, str(self._band), str(self._se), self.prob,self.descr)

    def __str__(self):
        """String summary over the object
        """
        return "ClusterBandPower, channel %i, band %s, start/end %s, probability %.3f" % (self._ch, str(self._band), str(self._se), self.prob)

class ClusterNd:
    """A Nd-cluster is represented by a boolean numpy-array, shape like the base data, and """
    def __init__(self,ar,p):
        pass
    #TODO: implement...

class ClusterDiscreteWavelets(object):
    """Representation of a cluster in wavelet-power,
    when linearized (stored in a 1d-array)"""
    def __init__(self,mask,prob=1.0):
        """Initialize the cluster. 
        :Parameters:
          mask: array
            array of booleans, True for all coefficients that are part of the cluster
        """
        self._mask = mask
        self._prob=prob

    @property
    def mask(self):
        return self._mask

    @property
    def mask_as_list(self):
        sc_lengths = []
        l = self.mask.shape[0]
        while l>1:
            l/=2
            sc_lengths.insert(0,l)
        sc_lengths.insert(0,l)
        assert sum(sc_lengths) == self.mask.shape[0], "The values automatically calculated for sc_lengths don't make sense!"
        rv = []
        for j in range(len(sc_lengths)):
            offset = sum(sc_lengths[:j])
            rv.append(self.mask[offset:offset+sc_lengths[j]])
        return rv


    def get_prob(self):
        return self._prob
    
    def set_prob(self,p):
        if not 0<=p<=1:
            raise ValueError("probabiliy must be in [0;1]")
        self._prob=float(p)

    def calc_cluster_stats(self,fs):
        """For some array of f-values, sum over all fs contained in the mask"""
        return fs[self._mask].sum()

    def is_connected_to(self,other):
        """Checks if cluster is connected to other cluster"""
        m1 = self.mask_as_list
        m2 = other.mask_as_list
        for i in range(len(m1)):
            if m1[i].sum()>0:
                for i1 in range(len(m1[i])):
                    if m1[i][i1]: #found an element of the cluster, now look in other
                        if m2[i][max(0,i1-1):min(i1+2,len(m2[i]))].sum()>0: #Same row
                            return True
                        try: #row below
                            if m2[i-1][i1/2]: #In row below, only one coeff is connected
                            #if m2[i-1][max(0,i1-1)/2:min(i1+2,len(m2[i-1]))/2].sum()>0: 
                                return True
                        except IndexError,ve:
                            pass #Just ignore Error
                        except ValueError,ve:
                            pass #Just ignore Error
                        try: #row above
                            if m2[i+1][i1*2] or m2[i+1][i1*2+1]: #In row above, two neighbors
                            #if m2[i+1][max(0,i1-1)*2:min(i1+2,len(m2[i+1]))*2].sum()>0: 
                                return True
                        except ValueError,ve:
                            pass #Just ignore Error
                        except IndexError,ve:
                            pass #Just ignore Error
        # If in all those lines no connection was found, return False
        return False
                            
    def union_with(self,other,do_check=True):
        """Checks if cluster is connected to other cluster and then returns union"""
        if do_check:
            if not self.is_connected_to(other):
                raise ValueError("Clusters without connection cannot be united")
        return ClusterDiscreteWavelets(self.mask+other.mask)

    probability = property(get_prob,set_prob)
        

class ClusterSearchDiscreteWavelet(ClusterSearch1d):
    """Search for clusters in data from discrete wavelet-transforms"""

    def __init__(self,array_list,stat_fun=None,threshold=1.67,num_surrogates=1000):
        """Initialization

           :Parameters:
             array_list: list 
               List of 2d-arrays containing the data, dim 1: linearized wavelet-powers, dim 2: elements of groups
             stat_fun : function
               function called to calculate statistics, must accept 1d-arrays as arguments (default: scipy.stats.f_oneway)
             threshold: float
               Threshold to use for deciding on statistical significance (uncorrected)
             num_surrogates: int
               number of permutation-surrogates to use for testing
        """
        #TODO: Do some checks
        ClusterSearch1d.__init__(self,array_list,stat_fun,threshold,num_surrogates)
    
    def find_clusters(self,ar,thres,cmp_type="greater"):
        """For a given 1d-array (test statistic), find all clusters which
        are above/below a certain threshold. Returns a list of 2-tuples.
        """
        #clusters =  []
        if not cmp_type in ["lower","greater","abs_greater"]:
            raise ValueError("cmp_type must be in [\"lower\",\"greater\",\"abs_greater\"]")
        ar = N.concatenate([N.array([thres]),ar,N.array([thres])])
        if cmp_type=="lower":
            ar_in = (ar<thres).astype(N.int)
        elif cmp_type=="greater":
            ar_in = (ar>thres).astype(N.int)
        else: #cmp_type=="abs_greater":
            ar_in = (abs(ar)>thres).astype(N.int)
        ar_switch = N.diff(ar_in)
        inpoints = list(N.arange(ar.shape[0])[ar_switch>0])
        outpoints = list(N.arange(ar.shape[0])[ar_switch<0])
        #Do not allow cluster for motherwavelet power
        if inpoints[0] == 0: 
            inpoints.pop(0)
            outpoints.pop(0)
        #print inpoints, outpoints
        #print inpoints, outpoints
        #in_out = N.concatenate([inpoints.reshape(-1,1),outpoints.reshape(-1,1)],axis=1)
        #clusters = [(c[0],c[1]) for c in in_out]

        #insert out/inpoints at cluster that span across scales (go over 512/256/128/...)
        points_to_insert = []
        for i in range(len(inpoints)):
            exp_in = int(np.log2(inpoints[i]))
            exp_out = int(np.log2(outpoints[i]))
            if exp_in!=exp_out:
                #print exp_in, exp_out
                for i_e,e in enumerate(range(exp_in+1,exp_out+1)):
                    #print e
                    points_to_insert.append((i+1+i_e,int(2**e),int(2**e-1)))
                    #print inpoints
                    #print outpoints
        #Now insert the new points, starting from the end
        #print inpoints[:10]
        #print outpoints[:10]
        points_to_insert.reverse()
        for p in points_to_insert:
            inpoints.insert(p[0],p[1])
            outpoints.insert(p[0]-1,p[2])
        #print "later"
        #print inpoints[:10]
        #print outpoints[:10]

        #Now make list of cluster
        clusters = []
        for i in range(len(inpoints)):
            mask = N.zeros((ar.shape[0]-2),np.bool) #-2: because of pre/appending of thres
            #print "mask.shape",mask.shape
            mask[inpoints[i]:outpoints[i]] = True
            clusters.append(ClusterDiscreteWavelets(mask))
        #Unite clusters if necessary
        goon_uniting = True
        while goon_uniting:
            for i in range(len(clusters)):
                restart = False
                for j in range(i+1,len(clusters)):
                    try:
                        clusters[i] = clusters[i].union_with(clusters[j])
                        clusters.pop(j)
                        restart = True
                        break
                    except ValueError:
                        restart=False
                if restart:
                    break
            if not restart:
                goon_uniting=False

        return clusters

    def calc_cluster_stats(self,cluster,fs):
        return cluster.calc_cluster_stats(fs)

    def _set_cluster_ps(self,cs,ps):
        for i in range(len(cs)):
            cs[i].probability = ps[i]

class Cluster2d(object):
    """Representation of a cluster in wavelet-power,
    when linearized (stored in a 1d-array)"""
    def __init__(self,mask,prob=1.0):
        """Initialize the cluster. 
        :Parameters:
          mask: array
            2d-array of booleans, True for all coefficients that are part of the cluster
        """
        self._mask = mask
        self._prob=prob

    @property
    def mask(self):
        return self._mask

    def get_prob(self):
        return self._prob
    
    def set_prob(self,p):
        if not 0<=p<=1:
            raise ValueError("probability must be in [0;1]")
        self._prob=float(p)

    def calc_cluster_stats(self,fs):
        """For some array of f-values, sum over all fs contained in the mask"""
        return fs[self._mask].sum()

    def is_connected_to(self,other):
        """Checks if cluster is connected to other cluster"""
        m1 = self.mask
        m2 = other.mask
        #for i1 in range(m1.shape[0]):
        #    for i2 in range(m1.shape[1]):
        #        if m1[i1,i2]:
        #            try:
        #                if m2[i1,i2-1]:
        #                    return True
        #            except Exception,e:
        #                pass
        #            try:
        #                if m2[i1,i2+1]:
        #                    return True
        #            except Exception,e:
        #                pass
        #            try:
        #                if m2[i1-1,i2]:
        #                    return True
        #            except Exception,e:
        #                pass
        #            try:
        #                if m2[i1+1,i2]:
        #                    return True
        #            except Exception,e:
        #                pass
        for i in range(m1.shape[1]): #iterate over frequencies
            if m1[:,i].sum()>0 and m2[:,max(0,i-1):min(m2.shape[1],i+2)].sum()>0:
                for i1 in range(m1.shape[0]):
                    if m1[i1,i]: #found an element of the cluster, now look in other
                        if m2[i1,max(0,i-1):min(m2.shape[1],i+2)].sum()>0: #Same row
                            return True
        # If in all those lines no connection was found, return False
        return False
                            
    def union_with(self,other,do_check=True):
        """Checks if cluster is connected to other cluster and then returns union"""
        if do_check:
            if not self.is_connected_to(other):
                raise ValueError("Clusters without connection cannot be united")
        return Cluster2d(self.mask+other.mask)

    def is_neighbor(self,i,j):
        """Check if a point is a neighbor of a cluster"""
        if self.mask[max(i-1,0):min(i+2,self.mask.shape[0]),max(j-1,0):min(j+2,self.mask.shape[1])].sum():
            return True
        else:
            return False

    def add_point(self,i,j):
        """Add a point at coordinates i,j, to cluster-mask"""
        try:
            self.mask[i,j] = True
        except ValueError:
            raise ValueError("coordinates (%i,%i) to be added to the cluster-mask are invalid"%(i,j))

    probability = property(get_prob,set_prob)
        

class ClusterSearch2d():
    """Search for clusters in data from Morlet wavelet analysis"""

    def __init__(self,array_list,stat_fun=None,threshold=1.67,num_surrogates=1000,strategy="permute",cmp_type="greater"):
        """Initialization

           :Parameters:
             array_list: list 
               List of 3d-arrays containing the data, dim 1: timepoints, dim 2: frequencies, dim 3: elements of groups
             stat_fun : function
               function called to calculate statistics, must accept 1d-arrays as arguments (default: scipy.stats.f_oneway)
             threshold: float
               Threshold to use for deciding on statistical significance (uncorrected)
             num_surrogates: int
               number of permutation-surrogates to use for testing
             strategy: "permute" or "pairwise"
               How to generate the surrogates. For unpaired tests, just permute. For paired tests on a per-subject level, 
               just exchange within subject.
        """
        #TODO: Do some checks
        self._al = array_list
        if stat_fun==None:
            self._sf = self.f_oneway
        else:
            self._sf = stat_fun
        self._threshold = threshold
        self._ns = num_surrogates
        self._strategy = strategy
        self._cmp_type = cmp_type
        #self._not_anova = not_anova

    #make read-only properties
    @property
    def array_list(self):
        return self._al
    
    @property
    def stat_fun(self):
        return self._sf
    
    @property
    def threshold(self):
        return self._threshold
    
    @property
    def num_surrogates(self):
        return self._ns

    #More interesting methods...
    def f_oneway(self,*args):
        """Call scipy.stats.f_oneway, but return only f-value"""
        return f_oneway(*args)[0]

    #THE method...
    def search(self):
        """For a list of 3d-arrays of data, e.g. power values, calculate some
        statistics for each timepoint (dim 1) and frequency (dim 2) over groups.  Do a cluster
        analysis with permutation test like in Maris, Oostenveld (2007)
        for calculating corrected p-values
        """
        #Create Shorthands
        al = self._al
        sf = self._sf

        #print len(al), [ar.shape for ar in al]
        ns_trs = [ar.shape[2] for ar in al] # Number of trials for each group
        print "ns_trs:", ns_trs
        if self._threshold == "f_auto":
            crit_f = fprobi(len(al)-1,N.sum(ns_trs)-1,0.05) #Critical F-value
        else:
            crit_f = self._threshold
        l1=al[0].shape[0]
        l2=al[0].shape[1]
        #print "l1,l2", l1,l2
        #print "CS1d: l=",l
        #Calculate Anova (or other stat_fun)
        fs = N.zeros((l1,l2),"d")
        for i1 in range(l1):
            for i2 in range(l2):   
                anova_ars = [ar[i1, i2,:] for ar in al]
                #print "anova_ars shapes", [ar.shape for ar in anova_ars]
                fs[i1,i2] = sf(*anova_ars)
        clusters = self.find_clusters(fs,crit_f,self._cmp_type)
        if len(clusters)>0:
            cluster_stats = [self.calc_cluster_stats(c,fs) for c in clusters] 
            cluster_ps = N.ones((len(clusters)),"d")
            cluster_stats_hist = N.zeros((self._ns)) #For making histogram (not visually) and finding percentile
            #ar_shuffle = N.concatenate(al,axis=2) #Vermeiden, Speicherprobleme!
            for i_s in range(self._ns):
                print "Surrogat Nr. %i"%i_s
                sys.stdout.flush()
                #t_start = time.time()
                ar_shuffle_list = self.make_surrogate_array_list(al,self._strategy)
                #print "Listen sortieren: %.2f"%(time.time()-t_start)
                #print "ar_shuffle_list shapes", [ar.shape for ar in ar_shuffle_list]
                #t_start = time.time()
                fs_surr = N.zeros((l1,l2),"d")
                for i1 in range(l1):
                    for i2 in range(l2):
                        anova_ars_perm = [ar[i1,i2,:] for ar in ar_shuffle_list]
                        fs_surr[i1,i2 ] = sf(*anova_ars_perm)
                #print "ANOVAs rechnen: %.2f"%(time.time()-t_start)
                #t_start = time.time()
                clusters_perm = self.find_clusters(fs_surr,crit_f,self._cmp_type)
                #print "Cluster finden: %.2f"%(time.time()-t_start)
                #print "clusters_perm", clusters_perm
                if len(clusters_perm)>0:
                    cluster_stats_perm  = [self.calc_cluster_stats(c,fs_surr) for c in clusters_perm]
                    cluster_stats_hist[i_s] = max(cluster_stats_perm)
                else:
                    cluster_stats_hist[i_s] = 0
            cluster_ps[:] = [percentileofscore(cluster_stats_hist,cluster_stats[i_cl]) for i_cl in range(len(clusters))]
            cluster_ps[:] = (100.0 - cluster_ps[:]) / 100.0 # From percent to fraction
            #print "CS1d: fs.shape[0]=",fs.shape[0]
            #Set cluster-ps in clusters
            self._set_cluster_ps(clusters,cluster_ps)
            return fs, N.array(clusters)[cluster_ps<=0.05], cluster_ps[cluster_ps<=0.05], N.array(clusters), cluster_ps
        else:
            return fs,N.array([]),N.array([]),N.array([]),N.array([])
    
    def find_clusters(self,ar,thres,cmp_type="greater"):
        """For a given 2d-array (test statistic), find all clusters which
        are above/below a certain threshold. 
        """
        #clusters =  []
        if not cmp_type in ["lower","greater","abs_greater"]:
            raise ValueError("cmp_type must be in [\"lower\",\"greater\",\"abs_greater\"]")
        #ar = N.concatenate([N.ones((1,ar.shape[1]))*thres,ar,N.ones((1,ar.shape[1]))*thres],axis=0)
        #p.figure(2)
        #p.imshow(ar)
        #time.sleep(10)
        clusters = []

        if cmp_type=="lower":
            ar_in = (ar<thres).astype(N.bool)
        elif cmp_type=="greater":
            ar_in = (ar>thres).astype(N.bool)
        else: #cmp_type=="abs_greater":
            ar_in = (abs(ar)>thres).astype(N.bool)
        
        #iteriere ueber Zeiten und Frequenzen, füge zu exisierenden Clustern hinzu
        already_visited = np.zeros(ar_in.shape,np.bool)
        for i_s in range(ar_in.shape[0]): #i_s wie i_sample
            for i_f in range(ar_in.shape[1]):
                if not already_visited[i_s,i_f]:
                    if ar_in[i_s,i_f]:
                        #print "Anzahl cluster:", len(clusters) 
                        mask = np.zeros(ar_in.shape,np.bool)
                        check_queue = Queue()
                        check_queue.put((i_s,i_f))
                        while not check_queue.empty():
                            pos_x,pos_y = check_queue.get()
                            if not already_visited[pos_x,pos_y]:
                                #print pos_x,pos_y
                                already_visited[pos_x,pos_y] = True
                                if ar_in[pos_x,pos_y]:
                                    mask[pos_x,pos_y] = True
                                    for coords in [(pos_x-1,pos_y),(pos_x+1,pos_y),(pos_x,pos_y-1),(pos_x,pos_y+1)]: #Direct Neighbors
                                        if 0<=coords[0]<ar_in.shape[0] and 0<=coords[1]<ar_in.shape[1]:
                                            check_queue.put(coords)
                        clusters.append(Cluster2d(mask))
        return clusters

    def make_surrogate_array_list(self,al,strategy):
        """Make new list of arrays for surrogate test"""
        ns_trs = [ar.shape[2] for ar in al] # Number of trials for each group
        if strategy=="permute":
            #Mache Liste mit Indices fuer alle Trials, permutiere, zerlege in Einzellisten der Laengen ns_trs
            indices_lists = N.split(N.random.permutation(sum(ns_trs)),N.cumsum(ns_trs)[:-1])
            ar_shuffle_list = []
            for il in indices_lists:
                #print "il", il
                ar = np.zeros((al[0].shape[0],al[0].shape[1],len(il)))
                for i_idx, idx in enumerate(il):
                    al_idx = 0
                    #print idx, al[al_idx].shape[2]
                    #print type(idx), type(al[al_idx].shape[2])
                    #print (idx>=(al[al_idx].shape[2]))
                    while idx>=(al[al_idx].shape[2]):
                        idx-=al[al_idx].shape[2]
                        al_idx+=1
                    ar[:,:,i_idx] = al[al_idx][:,:,idx]
                ar_shuffle_list.append(ar)
            return ar_shuffle_list
        elif strategy=="pairwise":
            ar_shuffle_list = []
            assert len(np.unique(ns_trs))==1, "For pairwise testing, all groups must have same length"
            for i in range(len(al)):
                ar_shuffle_list.append(np.zeros_like(al[i]))
                switch_list = np.permutation(range(len(ns_trs))) #In welche Reihenfolge umsortieren
                for j in range(al[i].shape[-1]):
                    ar_shuffle_list[-1][...,j] = al[i][...,switch_list[j]]
            return ar_shuffle_list
        else:
            raise ValueError("strategy must be one of 'permute' and 'pairwise'")

    def calc_cluster_stats(self,cluster,fs):
        return cluster.calc_cluster_stats(fs)

    def _set_cluster_ps(self,cs,ps):
        for i in range(len(cs)):
            cs[i].probability = ps[i]

if __name__ == "__main__":
    #from eegpy.analysis.wavelet import wavepower_lin
    #from eegpy.plot.wavelet import plot_wavedec_lin
    #import pylab as p
    #m1 = np.zeros(512,np.bool)
    #m1[89:110]=True
    #c1 = ClusterDiscreteWavelets(m1)
    #m2 = np.zeros(512,np.bool)
    #m2[140:180]=True
    #c2 = ClusterDiscreteWavelets(m2)
    #m3 = np.zeros(512,np.bool)
    #m3[256:270]=True
    #c3 = ClusterDiscreteWavelets(m3)
    #print c1.is_connected_to(c2), c2.is_connected_to(c1)
    #p.figure(1)
    #p.subplot(511)
    #plot_wavedec_lin(c1.mask,vmin=0,vmax=1,aspect="auto")
    #p.subplot(512)
    #plot_wavedec_lin(c2.mask,vmin=0,vmax=1,aspect="auto")
    #p.subplot(513)
    #plot_wavedec_lin(c3.mask,vmin=0,vmax=1,aspect="auto")
    #p.subplot(514)
    #c12 = c1.union_with(c2)
    #plot_wavedec_lin(c12.mask,vmin=0,vmax=1,aspect="auto")
    ##p.subplot(515)
    ##c123 = c12.union_with(c3)
    ##plot_wavedec_lin(c123.mask,vmin=0,vmax=1,aspect="auto")

    #p.figure(2)
    #p.subplot(311)
    #wts1 = np.zeros((1024,10))
    #for i in range(10):
    #    ar1 = np.random.random((1024))
    #    wts1[:,i] = wavepower_lin(ar1,"db4",normalise=False)
    #plot_wavedec_lin(wts1.mean(axis=1))
    #p.colorbar()
    #p.subplot(312)
    #wts2 = np.zeros((1024,10))
    #for i in range(10):
    #    ar2 = np.random.random((1024))
    #    #ar2 += np.sin(2*np.pi*60.*np.arange(0,1.024,0.001))
    #    ar2[0:300] += np.sin(2*np.pi*60.*np.arange(0,0.300,0.001))
    #    ar2[600:800] += np.sin(2*np.pi*30.*np.arange(0,0.200,0.001))
    #    wts2[:,i] = wavepower_lin(ar2,"db4",normalise=False)
    #plot_wavedec_lin(wts2.mean(axis=1))
    #p.colorbar()

    #cs = ClusterSearchDiscreteWavelet([wts1,wts2],num_surrogates=5)
    #fs,scs,scps,cs,cps = cs.search()
    ##clusters = cs.find_clusters((c12.mask+c3.mask).astype("d"),0.5)
    #p.subplot(313)
    #for i,c in enumerate(scs):
    #    print c,c.probability
    #    plot_wavedec_lin(np.ma.MaskedArray(np.ones(c.mask.shape)*(i+1),~c.mask),vmin=0,vmax=len(scs)+1,aspect="auto")

    #p.colorbar()
    #p.show()
    from scipy.stats import scoreatpercentile
    from eegpy.analysis.wavelet import wt_analyze
    import pylab as p
    import time
    freq1 = 10
    ampl1 = 5
    freq2 = 45
    ampl2 = 4
    xs = np.arange(0,2.5,0.001)
    ar1 = np.random.random((2500,20))*3
    ar2 = np.random.random((2500,20))*3
    mask1 = np.zeros((2500))
    mask1[1300:1800] += np.hanning(500)
    mask2 = np.zeros((2500))
    mask2[1100:1500] += np.hanning(400)
    for i_t in range(ar1.shape[1]):
        phase = np.random.random()*np.pi*2
        ar1[:,i_t] += np.sin(2*np.pi*freq1*xs+phase)*mask1*ampl1
        ar2[:,i_t] += np.sin(2*np.pi*freq2*xs+phase)*mask2*ampl2

    p.subplot(211)
    p.imshow(ar1.T,aspect="auto", interpolation="nearest")
    p.subplot(212)
    p.imshow(ar2.T,aspect="auto", interpolation="nearest")

    wt_freqs = np.linspace(1,50,50,True)
    wts1 = np.zeros((ar1.shape[0],len(wt_freqs),ar1.shape[1]),"d")
    wts2 = np.zeros((ar2.shape[0],len(wt_freqs),ar2.shape[1]),"d")
    for i_t in range(ar1.shape[1]):
        print i_t
        wts1[:,:,i_t] = abs(wt_analyze(ar1[:,i_t],wt_freqs,Fs=1000.0))**2
        wts1[:,:,i_t] /= wts1[500:800,:,i_t].mean(axis=0).reshape(1,-1).repeat(wts1.shape[0],axis=0)
        wts1[:,:,i_t] = 10*np.log10(wts1[:,:,i_t]) #Calculate dB?
        wts2[:,:,i_t] = abs(wt_analyze(ar2[:,i_t],wt_freqs,Fs=1000.0))**2
        wts2[:,:,i_t] /= wts2[500:800,:,i_t].mean(axis=0).reshape(1,-1).repeat(wts2.shape[0],axis=0)
        wts2[:,:,i_t] = 10*np.log10(wts2[:,:,i_t]) #Calculate dB?

    p.figure(2)
    p.subplot(2,1,1)
    p.imshow((wts1.mean(axis=2)).T,interpolation="nearest",aspect="auto",origin="lower",cmap=p.cm.jet)
    p.colorbar()
    p.subplot(2,1,2)
    p.imshow((wts2.mean(axis=2)).T,interpolation="nearest",aspect="auto",origin="lower",cmap=p.cm.jet)
    p.colorbar()
    #print "Pause"
    #time.sleep(5)
    #assert 1==0

    cs2d = ClusterSearch2d([wts1,wts2],num_surrogates=2)
    fs,cl_sig,cl_sig_ps,cl,cl_ps = cs2d.search()
    p.figure(3)
    vmax = scoreatpercentile(fs.flatten(),95)
    for i_c in range(len(cl_sig)):
        p.subplot(len(cl_sig),1,i_c+1)
        p.imshow(fs.T,interpolation="nearest",aspect="auto",origin="lower",cmap=p.cm.gray,vmin=0,vmax=4)
        p.contour(fs.T,[1.7],colors="r")
        p.imshow(np.ma.MaskedArray(fs,~cl_sig[i_c].mask).T,interpolation="nearest",aspect="auto",origin="lower",cmap=p.cm.jet,vmin=0,vmax=vmax)



