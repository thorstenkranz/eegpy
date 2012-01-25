#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

import os.path

#print __path__

icon_names = {"add_trigger": "add_trigger.png", "add_trigger_type":"add_trigger_type.png"}

def image_from_eegpy_stock(icon_name):
    rv = gtk.Image()
    rv.set_from_file(os.path.join(__path__[0],icon_names[icon_name]))
    return rv

def eegpy_logo(size="small"):
    """Returns the eegpy-logo as Pixbuf"""
    rv = None
    if size=="large":
        rv = gtk.gdk.pixbuf_new_from_file(os.path.join(__path__[0],"eegpy.png"))
    else:
        rv = gtk.gdk.pixbuf_new_from_file(os.path.join(__path__[0],"eegpy_small.png"))
    return rv