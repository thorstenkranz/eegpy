#!/usr/bin/env python
# -*- coding: utf-8 -*-


#eegpy-modules
try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.filter.freqfilt import filtfilt, butter
    from eegpy.ui.widgets.iowidgets import VarChooser
    from eegpy.ui.icon import eegpy_logo
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

import sys

class EegpyInfoDialog(gtk.Dialog):
    """A dialog to show some short infos"""
    def __init__(self, title="eegpy - information", text=""):
        gtk.Dialog.__init__(self, title, None, gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT, (gtk.STOCK_OK, gtk.RESPONSE_OK))    
        self.hbox = gtk.HBox()
        self.vbox.pack_start(self.hbox, padding=20)
        logo = gtk.Image()
        logo.set_from_pixbuf(eegpy_logo("large").scale_simple(147,134,gtk.gdk.INTERP_BILINEAR))
        self.hbox.pack_start(logo, padding=20)
        self.hbox.pack_start(gtk.Label(text), padding=20)
        self.vbox.show_all()
        #response = dialog_label.run()
        #dialog_label.destroy()

def show_info_dialog(title="eegpy - information", text=""):
    dialog = EegpyInfoDialog(title, text)
    response=dialog.run()
    dialog.destroy()

