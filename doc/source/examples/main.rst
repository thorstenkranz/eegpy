.. _examples_main:

==========
Examples
==========

In this section we'll show some examples of using eegpy. Usually, they 
are self-contained, i.e. necessary data are being generated artificially.

If not, the source of the data is mentioned.


.. toctree::
    :maxdepth: 2
    :hidden:
	
    using_F32
    erp
    wavelet_power

Reading and writing the F32 data format
===========================================================

F32 is a simple data format developed at the Univerity of Bonn, Department for
Epileptology. It's definition is freely available. eegpy offers a handy class,
:py:class:`eegpy.F32`, to easily read from and write to such files. 

:ref:`Look at this example for a quick overview<examples_f32>`

Plotting time-locked data, e.g., Event-Related-Potentials
=========================================================

In EEG analysis, especially cognitive studies, one often has repetitive events
and wants to analyse the event-related potentials. In :ref:`this section
<examples_erp>` we show how to do this in eegpy.

Analysing temporal and spectral fluctuations of power in a signal
===================================================================

Morlet wavelets are commonly used to get some insight on the temporal and
spectral properties of a signal. Also for this application we 
:ref:`have an example<examples_wavelet>`.


