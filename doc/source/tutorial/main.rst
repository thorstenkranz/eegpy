=========================
Tutorial - Starting Point
=========================

The intentions of this tutorial are to

* serve as a starting point to learn how to use eegpy
* show some of the key features of eegpy.

The basis for our tutrial is a little data set, taken from a simple
EEG study, featuring two different conditions.

We will start off by just clicking through the data (with *eegview*,
included in *eegpy*). Then, we will programmatically load these using
the built-in methods of *eegpy*. Afterwards, we're going to apply
some frequency filtering to the data. Then we'll do some segmentation,
extract the signals for each *trial* and plot the average signal for
each condition, the *event related potentials* (ERPs).

.. admonition:: Schedule

   .. toctree::
      :maxdepth: 2

      viewing
      loading
      filtering
      segmentation
      plotting
      wavelet
      cluster1d
      cluster2d
           
.. note::
        .. plot::
            :include-source:

            import matplotlib.pyplot as p
            import numpy as np

            p.plot(np.arange(0,10,0.1),
                   np.sin(np.arange(0,10,0.1)),
                   "rD",label="sine wave")

