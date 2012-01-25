#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""eegpy - processing and analysis of EEG-data
===============================================

Python-library for the processing and analysis of EEG/MEG-data and similar data.

.. packagetree::
   :style: UML

Dependencies
--------------

For the sake of rapid developement and to prevent me from reinventing the wheel, I have a handful of
external dependencies. These are:

  - numpy: Very important (oblig.)
  - scipy: Filter, Hilbert, Fourier usw. (oblig.)
  - matplotlib: Plotten von Daten, Topomaps etc.
  - pywt: wavelet-transfomration and decomposition
  - mdp: Contains some useful algorithms, especially PCA / ICA.
  - (biosig4python:) Reading EDF, BDF, GDF... Maybe include in distribution?
  - for special features some further dependencies might be made

Only numpy and scipy are obligatory. If other packages are missing, some parts of eegpy won't work.

Features
--------
  * Lesen von einigen wichtigen **Datenformaten**. Auf jeden Fall EDF und F32, andere können auch später erst dazu kommen. Schnittstelle muss einheitlich sein!
  * Datentyp für Zeitpunkte / Events
  * Export von Trial-Zeitreihen Arrays zur schnellen Wiederverwendung
  * Erkennung und Entfernen von Artefakten
  * Viele **Analysemethoden**, z.B. 
    * EKPs
    * Wavelet-Analyse
    * Verschiedene Spektrale Maße
    * Maße der nichtlinearen Zeitreihenanalyse
    * ...
    * Dabei sinnvolles Ablegen der Ergebnisse auf der Festplatte, um anschließende **Gruppenstatistik** zu vereinfachen.
  * **Visualisierung** von Ergebnissen
  * **Quellenlokalisation** (sLoreta u.a.)



Things to know:
  * If applicable, array-indices are always ordered like (idxDatapoints, idxChannels, furtherIndices...). furtherIndices include but are not restricted to the Index of Trials.
  * All filter functions designed for 1d-Arrays can be applied to arrays of arbitrary dimension. The filter is then applied to all 1d Arrays contained. See above for Index-Order.
"""

__docformat__ = "restructuredtext en"


from eegpy.formats.open_eeg import open_eeg
from eegpy.helper import load
from eegpy.formats.f32 import F32
from eegpy.events import EventTable

#from pylab import imshow
