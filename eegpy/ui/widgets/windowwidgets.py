#!/usr/bin/env python
# -*- coding: utf-8 -*-


#eegpy-modules
try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.filter.freqfilt import filtfilt, butter
    from eegpy.ui.widgets.iowidgets import VarChooser
    from eegpy.ui.icon import image_from_eegpy_stock, eegpy_logo
    from eegpy.ui.widgets.dialogwidgets import show_info_dialog
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

try:
    from matplotlib.axes import Subplot
    # uncomment to select /GTK/GTKAgg/GTKCairo
    from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
    from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
    import matplotlib
    #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, NavigationToolbar
    from matplotlib.figure import Figure, SubplotParams
    from matplotlib.axis import Axis
    import matplotlib.cm
except ImportError:
    raise FATALERROR('Error while importing matplotib. Please visit http://matplotlib.sf.net for more information.')

import sys


class EegpyBaseWin(gtk.Window):
    programName = "eegpy: Window-widget"
    panePosDef = 0
    x=None
    y=None
    
    # Konstruktor
    def __init__(self):
        gobject.threads_init()
        gtk.Window.__init__(self,gtk.WINDOW_TOPLEVEL)
        self.setupGUI()
        
    def cb_quit(self, b):
        self.window.hide()
        gtk.main_quit()
    
    # This callback quits the program
    def delete_event(self, widget, event, data=None):
        self.hide()
        gtk.main_quit()
        return False
            
    def setupGUI(self):
        self.set_default_size(700,500)
        self.set_title(self.programName)

        # Set a handler for delete_event that immediately
        # exits GTK.
        self.connect("delete_event", self.delete_event)

        # Sets the border width of the window.
        self.set_border_width(0)
        self.hpane = gtk.HPaned()
        self.add(self.hpane)
        #self.setupMenu()
        self.inner_pane = gtk.VPaned()
        self.hpane.add1(self.inner_pane)

        self.upper_hbox = gtk.HBox()
        self.lower_vbox = gtk.VBox()
        self.inner_pane.add1(self.upper_hbox)
        self.inner_pane.add2(self.lower_vbox)

        self.setupCanvas()        
        #self.setupLog()
        
        self.set_icon_list(eegpy_logo("small"), eegpy_logo("large"))
        self.show_all()
        self.panePosDef = self.hpane.get_position()
    
    def setupLog(self):
        self.frameLog = gtk.Frame("Log")
        self.logbuffer = gtk.TextBuffer(None)
        self.lvSw = gtk.ScrolledWindow()
        self.logview = gtk.TextView(self.logbuffer)
        self.logview.set_editable(False)
        #self.logview.set_buffer(self.logbuffer)
        self.logview.set_border_width(1)
        self.lvSw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC);
        self.lvSw.add(self.logview)
        self.frameLog.add(self.lvSw)
        self.upper_hbox.pack_end(self.frameLog)
        pass
    
    def log(self,message):
        """logs a message to the log window"""
        message=message+"\n"
        self.logbuffer.insert(self.logbuffer.get_start_iter(),message)#_at_cursor(message, len(message))
        #self.lvSw.set_vadjustment(self.lvSw.get_vadjustment().lower)
            
    def setupCanvas(self):
        self.f = Figure(figsize=(5,4), dpi=100, subplotpars=SubplotParams(left=0.06, top=0.95, right=0.97, bottom=0.1))
        self.setupSubplots()
        self.canvas = FigureCanvas(self.f)
        #self.canvas.show()
        self.lower_vbox.pack_start(self.canvas)
        #self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.toolbar = NavigationToolbar( self.canvas, self.window )
        self.lower_vbox.pack_start(self.toolbar, False, False)
    
    def setupSubplots(self):
        self.f.clear()
        self.a = self.f.add_subplot(111)
        self.subplAxes = self.f.get_axes()

class EegpyWin(EegpyBaseWin):
    def __init__(self, x_in=123456, x_out=123456):
        gobject.threads_init()
        self.x_in = x_in
        self.x_out = x_out
        EegpyBaseWin.__init__(self)
    
    def setupGUI(self):
        EegpyBaseWin.setupGUI(self)
        self.varChIn = VarChooser("Input", io="in", var_= self.x_in, onActivate=self.cb_Buttons)
        self.upper_hbox.pack_start(self.varChIn,expand=False,fill=False,padding=1)
        self.varChOut = VarChooser("Output", io="out", var_= self.x_out, onActivate=self.cb_Buttons)
        self.upper_hbox.pack_start(self.varChOut,expand=False,fill=False,padding=1)
        
        self.setupLog()
        
        self.show_all()
        self.panePosDef = self.hpane.get_position()