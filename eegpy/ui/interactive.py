#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################
# Module-Import #
#################

#eegpy-modules
try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.ui.FreqFiltWin import FreqFiltWin
    from eegpy.ui.PCAFiltWin import PCAFiltWin
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

#############################################
# Functions for use in interactive Sessions #
#############################################

def ui_freqFilt(x_in=None, x_out=None):
    """Loads the GUI for frequency-filtering"""
    ffw = FreqFiltWin(x_in=x_in, x_out=x_out)
    gtk.main()
    return ffw.varChOut.get_var()

def ui_pcaFilt(x_in=None, x_out=None):
    """Loads the GUI for pca-filtering"""
    pfw = PCAFiltWin(x_in=x_in, x_out=x_out)
    gtk.main()
    return pfw.varChOut.get_var()

if __name__ == "__main__":
    try:
        import eegpy.formats.edf
        r = eegpy.formats.edf.EdfReader("/home/thorsten/Documents/Downloads/Newtest17-256.bdf")
        ar2 = r.getData(0,2000)
        #ar3 = r.getData(0,1000)
    except Exception:
        pass
    #ar2_2 = None
    ar2_2 = ui_freqFilt(ar2)
    print ar2_2
    