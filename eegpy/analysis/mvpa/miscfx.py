# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PyMVPA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
"""Misc function performing operations on datasets.

All the functions defined in this module must accept dataset as the
first argument since they are bound to Dataset class in the trailer.
"""

__docformat__ = 'restructuredtext'

import numpy as N

from mvpa.datasets.base import Dataset, datasetmethod
from mvpa.base.dochelpers import table2string
from mvpa.misc.support import getBreakPoints
from mvpa.base import externals, warning

if __debug__:
    from mvpa.base import debug

if externals.exists('scipy'):
    from mvpa.datasets.miscfx_sp import detrend


@datasetmethod
def zscore(dataset, mean=None, std=None,
           perchunk=True, baselinelabels=None,
           pervoxel=True, targetdtype='float64'):
    """Z-Score the samples of a `Dataset` (in-place).

    ADAPTED FROM MVPA, BUT WITH RETURN-VALUES

    `mean` and `std` can be used to pass custom values to the z-scoring.
    Both may be scalars or arrays.

    All computations are done *in place*. Data upcasting is done
    automatically if necessary into `targetdtype`

    If `baselinelabels` provided, and `mean` or `std` aren't provided, it would
    compute the corresponding measure based only on labels in `baselinelabels`

    If `perchunk` is True samples within the same chunk are z-scored independent
    of samples from other chunks, e.i. mean and standard deviation are
    calculated individually.

    !!! Returns mean, std !!!
    """

    if __debug__ and perchunk \
      and N.array(dataset.samplesperchunk.values()).min() <= 2:
        warning("Z-scoring chunk-wise and one chunk with less than three "
                "samples will set features in these samples to either zero "
                "(with 1 sample in a chunk) "
                "or -1/+1 (with 2 samples in a chunk).")

    #perchunk not supported by my version
    if perchunk:
        raise ValueError("Zscoring perchunk not supported by this version")

    # cast to floating point datatype if necessary
    if str(dataset.samples.dtype).startswith('uint') \
       or str(dataset.samples.dtype).startswith('int'):
        dataset.setSamplesDType(targetdtype)

    def doit(samples, mean, std, statsamples=None):
        """Internal method."""

        if statsamples is None:
            # if nothing provided  -- mean/std on all samples
            statsamples = samples

        if pervoxel:
            axisarg = {'axis':0}
        else:
            axisarg = {}

        # calculate mean if necessary
        if mean is None:
            mean = statsamples.mean(**axisarg)

        # de-mean
        samples -= mean

        # calculate std-deviation if necessary
        # XXX YOH: would that be actually what we want?
        #          may be we want actually estimate of deviation from the mean,
        #          which per se might be not statsamples.mean (see above)?
        #          if logic to be changed -- adjust ZScoreMapper as well
        if std is None:
            std = statsamples.std(**axisarg)

        # do the z-scoring
        if pervoxel:
            samples[:, std != 0] /= std[std != 0]
        else:
            samples /= std

        return samples, std, mean

    if baselinelabels is None:
        statids = None
    else:
        statids = set(dataset.idsbylabels(baselinelabels))

    # for the sake of speed yoh didn't simply create a list
    # [True]*dataset.nsamples to provide easy selection of everything
    if perchunk:
        pass
        #for c in dataset.uniquechunks:
        #    slicer = N.where(dataset.chunks == c)[0]
        #    if not statids is None:
        #        statslicer = list(statids.intersection(set(slicer)))
        #        dataset.samples[slicer] = doit(dataset.samples[slicer],
        #                                       mean, std,
        #                                       dataset.samples[statslicer])
        #    else:
        #        slicedsamples = dataset.samples[slicer]
        #        dataset.samples[slicer] = doit(slicedsamples,
        #                                       mean, std,
        #                                       slicedsamples)
    elif statids is None:
        ns, std, mean = doit(dataset.samples, mean, std, dataset.samples)
    else:
        ns, std, mean = doit(dataset.samples, mean, std, dataset.samples[list(statids)])

    #Return mean, std for reuse
    #print "Mean:", mean, "Std:", std
    return mean, std

