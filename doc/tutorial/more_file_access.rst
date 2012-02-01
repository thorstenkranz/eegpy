.. _tutorial_more_file_access:

More on file access
-------------------

As mentioned before, the standard datafile format of eegpy is the f32 Format.
It has the advantage that it contains the data in a continous way that makes it
possible to address them via a numpy.memmap. This makes working with the eeg
easier and faster.


