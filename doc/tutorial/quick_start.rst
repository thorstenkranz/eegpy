.. _tutorial_quick_start:

Quick start
-----------
In this section we'll give you a "hands-on" example on how to use the different
parts of eegpy to process your eeg-data. To do all steps on your own, just grep
our example dataset from http://www.eegpy.org/data/eegpy_sample1.tar.gz,
extract the archive, go to the created folder and then open an interactive
python session. We recommend using `IPython <http://www.ipython.org>`_ which
offers significantly enhanced functionality compared to the standard
Python-Interpreter [#f2]_.

In this folder, you'll find the following files:
* dotlr.f32
* dotlr.vmrk

The first one is the data-file. It's in the *F32-format*, a format that is
based on header informations and float32-values. It is the standard-format of
eegpy, but eegpy supports many other formats. Additionally, converters are
included in eegpy.

The other one is a Marker-file as the Brain Products "Vision Analyzer" software
produces when you export data.

The data come from a simple EEG-study consisting of visual stimuli and ...

What we want to do now is a simple ERP (Event-Related-Potentials) analysis.
First, we want to apply a frequency-filter to our data, afterwards calculate
the ERPs for each stimulus type. In the end, we want to plot the results, i.e.
the ERPs for each channel and each stimuls type and additionally we want to
have a look at the topography of the P300-component which we will encounter in
our analysis.

Let's start...
^^^^^^^^^^^^^^
First, we do the necessary imports::

    import eegpy

Accessing the data
^^^^^^^^^^^^^^^^^^
Next, we'll open the data file::

    eeg = eegpy.F32("dotlr.f32","r")

This will open it read-only, so that we do not change our original data.
Accessing your data now is very intuitive, basically like accessing a NumPy
array. You can do slicing on it: the first index is for *datapoints*, the
second one for *channels*. It is good convention (and also used throughout
eegpy) that the first index is used to address the datapoint /
timepoint.

Say, you want to get all data of channel 3 [#f1]_, you just type::

    data = eeg[:,3]

This will put all data of channel three into the variable ``data``. To check
this, do a simple plot command using the `Matplotlib
<http://matplotlib.sourceforge.net>`_::
    
    plot(data[10000,12000])

which should show a plot of the timepoints 10000 through 12000 of channel 3.
Something like::

    plot(eeg[10000,12000,:])

will plot all channels in the same time-window.

Filtering the data
^^^^^^^^^^^^^^^^^^
Applying some frequency-filter to our data is very easy. To do this, we should
first make a copy of our data::

    eeg.copy_to_file("dotlr_f.f32")
    out = eegpy.F32("dotlr_f.f32","r+")

This will copy the contents of our file to the new file named ``"dotlr_f.f32"``
which we open as "read-write" afterwards.

The actual process of filtering is very easy now::

    from eegpy.filter import filtfilt_low
    out[...] = filtfilt_low(70.0,out[...])

We just import the function ``filtfilt_low`` [#f3]_ and use to to filter the contents of our EEG-file with a lowpass-filter with a corner-frequency of 70.0 Hz. That's it.

Parsing the Triggers
^^^^^^^^^^^^^^^^^^^^
Next we want to read the information about the stimulus timing from the
marker-file. Managing such events is very important. eegpy offers a special
class to do this, ``eegpy.EventTable``. It can be constructed from multiple
data-formats, from dictionaries and more.

Reading the *.vmrk*-file is as easy as::

    events = eegpy.EventTable("dotlr.vmrk")

This can now be used in a dictionary-like fashion::

    In [6]: print event.keys()
    ['', 'S 41', 'S 42', 'S101', 'S 31', 'S  1', 'S 32', 
    'S  2', 'S 12', 'S 11', 'S 22', 'S 21']

    In [7]: print event_table["S101"]
    [43066, 48369, 56051, ... , 1195140, 1199566]

``EventTable``-objects offer a number of convenience-functions, e.g. for conversion if you extract parts of an EEG. For more information, see ???.

Calculating ERPs
^^^^^^^^^^^^^^^^
After these preparational steps, we can go and calculate ERPs from our data. 



.. rubric:: Footnotes
.. [#f2] From now on, we assume you're using IPython, started with something
   like `ipython -pylab`. If you don't have IPython available, you need to type
   ``from pylab import *`` and ``ion()`` at the beginning of your session to
   enable the plotting abilities of pylab.
.. [#f1] The counting is zero based, just as it is the case for all Python lists and arrays.
.. [#f3] This method uses the Butterworth-algorithm (scipy.signal.butter) to calculate
   filter coefficients. It then filters the data forwards and backwards in time
   in order to preserve phases (comp. ???)

