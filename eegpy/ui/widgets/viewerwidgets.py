# -*- coding: utf-8 -*-
import os.path
import sys
import time
from threading import Thread, Lock
import pickle

#Eigene Imports
import eegpy
from eegpy.formats import f32
from eegpy.misc import FATALERROR
from eegpy.ui.icon import image_from_eegpy_stock, eegpy_logo
from eegpy.ui.widgets.dialogwidgets import show_info_dialog
from markerwidgets import Marker, MarkerWithAverage
from trigmanwidgets import TriggerManager
from eegpy.analysis.wavelet import wavedec_lin


try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

from numpy import *

try:
    import matplotlib
    from matplotlib.axes import Subplot
    from matplotlib.mlab import specgram
    # uncomment to select /GTK/GTKAgg/GTKCairo
    from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
    # or NavigationToolbar for classic
    from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
    #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, NavigationToolbar
    from matplotlib.figure import Figure, SubplotParams
    from matplotlib.axis import Axis
    import matplotlib.cm as cm
except ImportError:
    raise FATALERROR('Error while importing matplotib. Please visit http://matplotlib.sf.net for more information.')

try:
    import pywt
except ImportError:
    raise FATALERROR('PyWavelet not found.')

class EEGPlot:
    programName = "eegpy Viewer 0.0.3"
    offset=100
    f32Loaded = False
    chList = []
    data=[]
    ovData = []
    ovTs = []
    ts=[]
    fn = "/media/Extern/public/Experimente/AudioStroop/2008-01-21_Leon_AS/eeg/leon_leerlauf_f.f32" #None
    reader = None
    markerfilename = None
    markerparser = None
    markMaker=None
    trigManager=None
    timesOffset = 0
    tsFactor = 1 #Faktor für Umrechnung der Zeiten, wird wichtig beim zeichnen der Markers
    boolPlotMarkers = True
    OvThread = None
    f32Lock = Lock()
    panePosDef = 0
    showAnalysis = False # Wenn True, dann wird eine Analyse des ersten ausgewählten Kanals gezeigt.
    whichAnal = "Specgram"
    
    # Konstruktor
    def __init__(self,fn=None):
        gobject.threads_init()
        self.setupGUI()
        if fn==None:
            self.set_filename("/media/Extern/public/Experimente/AudioStroop/2008-01-21_Leon_AS/eeg/leon_leerlauf_f.f32")
        else:
            self.set_filename(fn)
    # Our callback.
    # The data passed to this method is printed to stdout
    def callback(self, widget, data=None):
        print "Hello again - %s was pressed" % data
    
    def cb_plot(self,widget,data=None):
        self.plot()
        
    def cb_plotOverview(self,widget,data=None):
        self.plot_overview()
        
    def cb_canvas(self, event):    
        if event.inaxes == self.a2:
            #print "In Übersicht"            
            self.sbStartpoint.set_value(event.xdata)
            self.plot()
        elif event.inaxes == self.a:
            #print event.xdata, event.ydata
            if self.markMaker != None:
                self.markMaker.add(event.xdata)
            self.plot_data()
            #pass
            #mrkr = Marker()
    
    def cb_col_toggled( self, cell, path, user_data):
        model, column = user_data
        model[path][column] = not model[path][column]
        #print model, path, model[path][column]
        if column == 2:
            self.plot()
        else:
            self.plot_data()
        return     
    
    def cb_MvAccels(self, action):
        #print "Hallo", 
        #print action.get_name()
        if not self.f32Loaded:
            return False   
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        
        if action.get_name() == "Llleft":
            start=start-int(5*length*stride)
        elif action.get_name() == "Lleft":
            start=start-int(1*length*stride)
        elif action.get_name() == "Left":
            start=start-int(0.2*length*stride)
        elif action.get_name() == "Right":
            start=start+int(0.2*length*stride)
        elif action.get_name() == "Rright":
            start=start+int(1*length*stride)
        elif action.get_name() == "Rrright":
            start=start+int(5*length*stride)
        if start < 0:
            start=0
        if start+stride*length>self.reader.numDatapoints:
            start=self.reader.numDatapoints-stride*length
        self.sbStartpoint.set_value(start)
        self.plot()
         
    def cb_MvButtons(self, widget):
        if not self.f32Loaded:
            return False   
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        
        if widget == self.btLlleft:
            start=start-int(5*length*stride)
        elif widget == self.btLleft:
            start=start-int(1*length*stride)
        elif widget == self.btLeft:
            start=start-int(0.2*length*stride)
        elif widget == self.btRight:
            start=start+int(0.2*length*stride)
        elif widget == self.btRright:
            start=start+int(1*length*stride)
        elif widget == self.btRrright:
            start=start+int(5*length*stride)
        if start < 0:
            start=0
        if start+stride*length>self.reader.numDatapoints:
            start=self.reader.numDatapoints-stride*length
        self.sbStartpoint.set_value(start)
        #self.plot()
        
    def cb_SelAccels(self, action):
        def select(model, path, iter,data):
            if data=="all":
                model.set_value(iter,2,True)
            elif data=="none":
                model.set_value(iter,2,False)
            return False     # keep the foreach going
        
        #print "Hallo", 
        #print action.get_name()
        if not self.f32Loaded:
            return False   
        
        if action.get_name() == "PlotAll":
            self.tree.foreach(select,"all")
        elif action.get_name() == "PlotNone":
            self.tree.foreach(select,"none")
        self.setup_subplots()
        self.plot()
    
    def cb_activate_radio_action(self, action, current):
        self.whichAnal = current.get_name()
        self.plot_data()
        #print 'Radio action "%s" selected, %s'% (current.get_name(),action.get_name())

    def cb_SpinBs(self, widget):
        if not self.f32Loaded:
            return False  
        
        if widget in [self.sbStartpoint, self.sbNDataP, self.sbStride]:
            self.setup_subplots()
            self.plot()
        elif widget == self.sbScale:
            self.offset = self.sbScale.get_value()
            #print "Offset: %f" % self.offset
            self.plot_data()
    
    def cb_showAnalysis(self, widget):
        if widget == self.btAnal:
            self.showAnalysis = not self.showAnalysis
            self.setup_subplots()
            self.plot()
    
    
    def cb_open(self, b):
        dialog = gtk.FileChooserDialog("Open f32-File..", None, gtk.FILE_CHOOSER_ACTION_OPEN, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("EEG f32")
        filter.add_pattern("*.f32")
        dialog.add_filter(filter)
        
        filter = gtk.FileFilter()
        filter.set_name("All files")
        filter.add_pattern("*")
        dialog.add_filter(filter)
        
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            self.set_filename(dialog.get_filename())
            #print dialog.get_filename(), 'selected'
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()
        pass
    
    def cb_startMarking(self, b):
        if self.markMaker == None:
            self.markMaker = MarkerWithAverage(self)
        else:
            try:
                self.markMaker.window.show_all()
            except AttributeError,e:
                "markMaker ist irgendwie falsch..."
                
    def cb_manageTriggers(self, b):
        if self.trigManager == None:
            self.trigManager = TriggerManager(self)
        else:
            try:
                self.trigManager.window.show()
                #self.trigManager.mainBox.show_all()
            except AttributeError,e:
                "trigManager ist irgendwie falsch..."
    
    def cb_about(self, b):
        text = """This is the %s. 

It is the EEG-Viewer supplied with eegpy, an open-source project
for the analysis of eeg-data.

More informations are available at our website, http://eegpy.sf.net."""%(self.programName)
        show_info_dialog("About",text)
        
    def cb_quit(self, b):
        self.window.hide()
        gtk.main_quit()
    
    def addLine(self, widget, data=None):
        #print "Hello again - %s was pressed" % data
        iter = self.tree.insert_before(None, None)
        self.tree.set_value(iter, 0, 1)
        self.tree.set_value(iter, 1, data)

    # This callback quits the program
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False
    
                
    def set_filename(self, fn):
        """Set the filename to use, check if file exists and start f32reader """
        self.fn = fn
        self.chList = []
        self.data=[]
        self.ovData = []
        self.ovTs = []    
        self.f32Loaded = False
        if os.path.exists(self.fn):
            try:
                self.reader = f32.F32(self.fn,"r")
                self.f32Loaded = True
                self.window.set_title(self.programName + " - " + self.fn)
                #Kanalliste füllen
                self.tree.clear()
                #base = self.tree.append(None)
                #self.tree.set(base ,1,"Channels:")
                i=0
                for c in self.reader.channel_names:                                
                    iter = self.tree.append(None)
                    self.tree.set(iter, 0,i , 1,c, 2,True, 3,False)
                    #print i, c
                    #self.tree.set_value(iter, 0, c)
                    i += 1
                self.pane.set_position(self.panePosDef)
                #self.scale_startpoint["to"] = self.reader.numDatapoints-self.scale_nDataP["from"]-1 
                self.sbStartpoint.set_adjustment(gtk.Adjustment(0,0,self.reader.numDatapoints-100-1,100,1000))
                self.a.clear()
                #self.plot()
                #Übersicht plotten
                #gtk.gdk.threads_enter()    
                #self.OvThread = Thread(target=self.getOvData)
                #self.OvThread.start()
                self.get_overview_data()
                evt_exts = [".evt",".vmrk"]
                for e in evt_exts:
                    if os.path.exists(os.path.splitext(self.fn)[0]+e):
                        self.cb_manageTriggers(None)
                        self.trigManager.add(os.path.splitext(self.fn)[0]+e)
                        break
                self.plot()
                #gtk.gdk.threads_leave()
            except Exception, e:
                self.fn= None
                #tkMessageBox.showwarning("Falsches Format","Dies ist keine gültige f32-Datei!")
                print e
        
    def setupGUI(self):
        # Create a new window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_default_size(700,500)
        self.window.set_title(self.programName)

        # Set a handler for delete_event that immediately
        # exits GTK.
        self.window.connect("delete_event", self.delete_event)

        # Sets the border width of the window.
        self.window.set_border_width(0)
        self.mainBox = gtk.VBox()
        self.window.add(self.mainBox)
        self.setup_menu()
        self.hpane = gtk.HPaned()
        self.mainBox.pack_start(self.hpane)
        self.pane = gtk.VPaned()
        self.hpane.add1(self.pane)
        #self.window.add(self.pane)
        #self.upperTable = gtk.Table(1, 3)
        self.bTable = gtk.Table(3, 6, True)
        self.vbox = gtk.VBox()
        self.hbox = gtk.HBox()
        self.hbox2 = gtk.HBox()
        self.pane.add1(self.hbox)
        self.pane.add2(self.hbox2)
        self.hbox2.pack_start(self.vbox)
        
        self.setup_SpinBoxes()
        self.setup_buttons()
        self.setup_canvas()
        self.vbox.pack_start(self.toolbar, False, False)
        self.setup_TreeView()
        self.window.set_icon_list(eegpy_logo("small"))#, eegpy_logo("large"))
        self.window.show_all()
        self.panePosDef = self.pane.get_position()
    
    
    def setup_menu(self):
        self.ui = '''<ui>
    <menubar name="MenuBar">
      <menu action="File">
        <menuitem action="Open"/>
        <menuitem action="Marks"/>
        <menuitem action="Triggers"/>
        <separator/>
        <menuitem action="Quit"/>
      </menu>
      <menu action="Movement">
        <menuitem action="Llleft"/>
        <menuitem action="Lleft"/>
        <menuitem action="Left"/>
        <menuitem action="Right"/>
        <menuitem action="Rright"/>
        <menuitem action="Rrright"/>
      </menu>
      <menu action="ChSel">
        <menuitem action="PlotAll"/>
        <menuitem action="PlotNone"/>
        <menuitem action="PlotMarked"/>
        <menuitem action="NotPlotMarked"/>
      </menu>
      <menu action="Analysis">
        <menuitem action="Specgram"/>
        <menuitem action="XCorr"/>
        <menuitem action="Avg"/>
        <menuitem action="Wavelet"/>
      </menu>
      <menu action="Help">
        <menuitem action="About"/>
      </menu>
    </menubar>
    </ui>'''
        # Create a UIManager instance
        self.uimanager = gtk.UIManager()
        
        # Add the accelerator group to the toplevel window
        self.accelgroup = self.uimanager.get_accel_group()
        self.window.add_accel_group(self.accelgroup)
        
        # Create an ActionGroup
        actiongroup = gtk.ActionGroup('UIManagerExample')
        self.actiongroup = actiongroup

        # Create a ToggleAction, etc.
        #actiongroup.add_toggle_actions([('Mute', None, '_Mute', '<Control>m',
         #                                'Mute the volume', self.mute_cb)])

        # Create actions
        actiongroup.add_actions([('Open', gtk.STOCK_OPEN, '_Open f32-file', None, 'Open an f32-file from the file-system', self.cb_open),                                  
                                 ('Marks', None, '_Mark time-points', "<Control>m", 'Mark timepoints in file for later use', self.cb_startMarking),
                                 ('Triggers', None, 'Manage _trigger-files', "<Control>t", 'Manage the external trigger-files used for displaying certain events in the eeg.', self.cb_manageTriggers),
                                 ('Quit', gtk.STOCK_QUIT, '_Quit me!', None,'Quit the Program', self.cb_quit),
                                 ('File', None, '_File'),
                                 ('Movement', None, '_Movement'),
                                 ('Analysis', None, '_Analysis'),
                                 ('Mute', None, '_Mute'),
                                 ('Llleft', None, "5x Left" , "<Control><Alt><Shift>Left", "Llleft", self.cb_MvAccels),
                                 ('Lleft', None, "1x Left" , "<Control><Alt>Left", "Lleft", self.cb_MvAccels),
                                 ('Left', None, "0.2x Left" , "<Control>Left", "Left", self.cb_MvAccels),
                                 ('Right', None, "0.2x Right" , "<Control>Right", "Right", self.cb_MvAccels),
                                 ('Rright', None, "1x Right" , "<Control><Alt>Right", "Rright", self.cb_MvAccels),
                                 ('Rrright', None, "5x Right" , "<Control><Alt><Shift>Right", "Rrright", self.cb_MvAccels)])
        
        actiongroup.add_actions([('ChSel', None, '_Channel-Selection'),                                  
                                 ('PlotAll', None, 'Mark all channels for plotting', None, 'Mark all channels for plotting', self.cb_SelAccels),                                  
                                 ('PlotNone', None, "Unmark all channels" , None, None, self.cb_SelAccels),
                                 ('PlotMarked', None, "Mark selected channels" , None, None, self.cb_SelAccels),
                                 ('NotPlotMarked', None, "Unmark selected" , None, None, self.cb_SelAccels)])
        
        actiongroup.add_actions([('About', None, "About" , None, None, self.cb_about),
                                 ('Help', None, '_Help')])
        
        (ANAL_SPEC,ANAL_XCORR,ANAL_AVG) = range(3)
        
        anal_entries = (
          ( "Specgram", None,                               # name, stock id
            "_Spectogram", None,                      # label, accelerator
            "Show a spectogram in the Analysis-window", ANAL_SPEC ),                      # tooltip, value
          ( "XCorr", None,                             # name, stock id
            "_Cross-correlation", None,                    # label, accelerator
            "Compute the cross-correlation", ANAL_XCORR ),                    # tooltip, value
          ( "Avg", None,                             # name, stock id
            "_Average over channels", None,                    # label, accelerator
            "Compute the average over channels", ANAL_AVG ),                    # tooltip, value
          ( "Wavelet", None,                             # name, stock id
            "_Wavelet", None,                    # label, accelerator
            "Wavelet-decomposition of the signal", ANAL_AVG )
        )

        actiongroup.add_radio_actions(anal_entries, ANAL_SPEC, self.cb_activate_radio_action)

        #actiongroup.get_action('Quit').set_property('short-label', '_Quit')
        
        # Add the actiongroup to the uimanager
        self.uimanager.insert_action_group(actiongroup, 0)
        
        # Add a UI description
        self.uimanager.add_ui_from_string(self.ui)
        
        # Create a MenuBar
        self.menuBar = self.uimanager.get_widget('/MenuBar')
        
        self.frMenuBar = gtk.Frame()
        self.frMenuBar.set_property("shadow-type", gtk.SHADOW_OUT)
        #self.menuBar = gtk.MenuBar()
        self.frMenuBar.add(self.menuBar)
        #self.frMenuBar.show_all()
        self.mainBox.pack_start(self.frMenuBar,expand=False)
               
    def setup_SpinBoxes(self):
        self.frameSBs = gtk.Frame("Adjust View")
        self.hbox.pack_start(self.frameSBs,expand=False,padding=50)
        self.tableSBs = gtk.Table(4,2,False)
        self.frameSBs.add(self.tableSBs)
        self.tableSBs.attach(gtk.Label("Start"),0,1,0,1)
        self.sbStartpoint = gtk.SpinButton(gtk.Adjustment(0,0,10000,100,1000))
        self.sbStartpoint.connect("value_changed",self.cb_SpinBs)
        self.tableSBs.attach(self.sbStartpoint,1,2,0,1)
        self.tableSBs.attach(gtk.Label("Datapoints"),0,1,1,2)
        self.sbNDataP = gtk.SpinButton(gtk.Adjustment(1100,100,20000,100,1000))
        self.sbNDataP.connect("value_changed",self.cb_SpinBs)
        self.tableSBs.attach(self.sbNDataP,1,2,1,2)
        self.tableSBs.attach(gtk.Label("Stepsize"),0,1,2,3)
        self.sbStride = gtk.SpinButton(gtk.Adjustment(5,1,400,1,100))
        self.sbStride.connect("value_changed",self.cb_SpinBs)
        self.tableSBs.attach(self.sbStride,1,2,2,3)
        self.tableSBs.attach(gtk.Label("Scale"),0,1,3,4)
        self.sbScale = gtk.SpinButton(gtk.Adjustment(100,1,10000,10,100))
        self.sbScale.connect("value_changed",self.cb_SpinBs)
        self.tableSBs.attach(self.sbScale,1,2,3,4)
        for i in [self.sbStartpoint, self.sbNDataP, self.sbStride, self.sbScale]:
            i.set_numeric(True)
        
    def setup_buttons(self):
        # Create move-buttons
        self.btLlleft = gtk.Button("<<<")        
        self.btLlleft.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btLlleft, 0, 1, 0, 1)
        #self.btLlleft.show()
        
        self.btLleft = gtk.Button("<<")
        self.btLleft.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btLleft, 1, 2, 0, 1)
        #self.btLleft.show()
        
        self.btLeft = gtk.Button("<")        
        self.btLeft.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btLeft, 2, 3, 0, 1)
        #self.btLeft.show()
        
        self.btRight = gtk.Button(">")
        self.btRight.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btRight, 3, 4, 0, 1)
        #self.btRight.show()
        
        self.btRright = gtk.Button(">>")
        self.btRright.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btRright, 4, 5, 0, 1)
        #self.btRright.show()
        
        self.btRrright = gtk.Button(">>>")
        self.btRrright.connect("clicked", self.cb_MvButtons)
        self.bTable.attach(self.btRrright, 5, 6, 0, 1)
        #self.btRrright.show()
        
        # Create "Plot" buttons
        self.btPlot = gtk.Button("Plot")
        self.btPlot.connect("clicked", self.cb_plot)
        self.bTable.attach(self.btPlot, 0, 3, 1, 2)
        #self.btPlot.show()     
        
        self.btPlotOv = gtk.Button("Plot Overview")
        self.btPlotOv.connect("clicked", self.cb_plotOverview)
        self.bTable.attach(self.btPlotOv, 3, 6, 1, 2)
        #self.btPlotOv.show()   
        
        #Radio-Buttons
        self.lbTime = gtk.Label("Time:")
        self.bTable.attach(self.lbTime,0,1,2,3)
        self.rbSamples = gtk.RadioButton(None,"#")
        self.rbSamples.connect("clicked", self.cb_plot)
        self.bTable.attach(self.rbSamples,1,2,2,3)
        self.rbSecs = gtk.RadioButton(self.rbSamples,"s")
        self.rbSecs.connect("clicked", self.cb_plot)
        self.bTable.attach(self.rbSecs,2,3,2,3)
        self.rbMins = gtk.RadioButton(self.rbSamples,"m")
        self.rbMins.connect("clicked", self.cb_plot)
        self.bTable.attach(self.rbMins,3,4,2,3)
        self.rbHours = gtk.RadioButton(self.rbSamples,"h")
        self.rbHours.connect("clicked", self.cb_plot)
        self.bTable.attach(self.rbHours,4,5,2,3)
        self.rbDays = gtk.RadioButton(self.rbSamples,"d")
        self.rbDays.connect("clicked", self.cb_plot)
        self.bTable.attach(self.rbDays,5,6,2,3)
        
        self.hbox.pack_start(self.bTable, False, padding=50)
        
        self.bTable2 = gtk.Table(1, 1, True)        
        self.btAnal = gtk.Button("Show/Hide Analysis")
        self.btAnal.connect("clicked", self.cb_showAnalysis)
        self.bTable2.attach(self.btAnal,0,1,0,1)
        self.hbox.pack_start(self.bTable2, False, padding=50)
        
    def setup_canvas(self):
        self.f = Figure(figsize=(5,4), dpi=100, subplotpars=SubplotParams(left=0.06, top=0.95, right=0.97, bottom=0.1,hspace=0))
        self.setup_subplots()
        self.canvas = FigureCanvas(self.f)
        #self.canvas.show()
        self.vbox.pack_start(self.canvas)
        #self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.toolbar = NavigationToolbar( self.canvas, self.window )
        self.canvas.mpl_connect('button_press_event', self.cb_canvas)
    
    def setup_subplots(self):
        self.f.clear()
        if self.showAnalysis:
            self.a = self.f.add_subplot(212)
            self.a2 = self.f.add_subplot(414)
            self.a3 = self.f.add_subplot(211)
            self.subplAxes = self.f.get_axes()
            self.subplAxes[0].set_position([0.06,0.05,0.91,0.50])
            self.subplAxes[2].set_position([0.06,0.55,0.91,0.35])
            self.subplAxes[1].set_position([0.06,0.96,0.91,0.02])
            self.subplAxes[1].set_xticks([])
            self.subplAxes[1].set_yticks([])
        else:
            self.a = self.f.add_subplot(111)
            self.a2 = self.f.add_subplot(414)
            self.a3 = None
            self.subplAxes = self.f.get_axes()
            self.subplAxes[1].set_position([0.06,0.96,0.91,0.02])
            self.subplAxes[1].set_xticks([])
            self.subplAxes[1].set_yticks([])
        
        #self.a.autoscale_view(tight=True)
        #self.subplAxes[0].set_position(0,0.1,1,0.9)
        #self.subplAxes[1].set_aspect(0.01)
        #self.subplAxes[1].set_yscale('log')        
    
    def setup_TreeView(self):
        self.tvScrolledWin = gtk.ScrolledWindow()
        self.tvScrolledWin.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.tree = gtk.TreeStore(gobject.TYPE_INT, gobject.TYPE_STRING, gobject.TYPE_BOOLEAN, gobject.TYPE_BOOLEAN)
        self.treeV = gtk.TreeView(self.tree)
        self.treeV.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        renderer = gtk.CellRendererText()
        self.col1 = gtk.TreeViewColumn("Number", renderer,text=0)
        self.col1.set_resizable(True)
        self.col1.set_min_width(20)
        #self.col1.set_sizing(gtk.TREE_VIEW_COLUMN_AUTOSIZE)
        self.treeV.append_column(self.col1)
        self.col2 = gtk.TreeViewColumn("Channel Name", renderer, text=1)
        self.col2.set_resizable(True)
        self.col2.set_min_width(20)
        self.treeV.append_column(self.col2)
        #Dritte Spalte für Highlight
        renderer1 = gtk.CellRendererToggle()
        renderer1.set_property('activatable', True)
        renderer1.connect("toggled", self.cb_col_toggled, (self.tree, 2))
        self.col3 = gtk.TreeViewColumn("Plot?", renderer1, active=2)
        self.col3.set_resizable(True)
        self.col3.set_min_width(20)
        self.treeV.append_column(self.col3)
        renderer2 = gtk.CellRendererToggle()
        renderer2.set_property('activatable', True)
        renderer2.connect('toggled', self.cb_col_toggled, (self.tree, 3))
        self.col4 = gtk.TreeViewColumn("-", renderer2, active=3)
        self.col4.set_resizable(True)
        self.col4.set_min_width(20)
        self.treeV.append_column(self.col4)
        #self.treeV.show()
        self.tvScrolledWin.add(self.treeV)
        #self.tvScrolledWin.show_all()
        #self.hbox.pack_start(self.tvScrolledWin)
        self.hpane.add2(self.tvScrolledWin)
        
           
    def plot(self):
        if not self.f32Loaded:
            return False
        #Zuerst Start und Stride zurecht schieben
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        
        if length*stride>self.reader.numDatapoints:
            stride = self.reader.numDatapoints/length
            self.sbStride.set_value(stride)
        
        if start+stride*length>self.reader.numDatapoints:
            start=self.reader.numDatapoints-stride*length
            self.sbStartpoint.set_value(start)
        #Dann DAten holen
        #tStart = time.time()
        self.get_data()
        #print "getData: Dauer: ", time.time()-tStart, "s"
        if self.chList == []:
            return False
        #Dann Daten plotten
        #tStart = time.time()
        self.plot_data() 
        #print "plotData. Dauer: ", time.time()-tStart, "s"
        try:
            self.plot_overview_data()   
            pass
        except Exception:
            pass
        
    def plot_overview(self):
        if not self.f32Loaded:
            return False
        #Zuerst Start und Stride zurecht schieben
        self.sbStartpoint.set_value(0)
        if (self.sbStride.get_range()[1]*self.sbNDataP.get_value())<self.reader.numDatapoints:
            self.sbNDataP.set_value(self.reader.numDatapoints/self.sbStride.get_range()[1])
            self.sbStride.set_value(self.sbStride.get_range()[1])
        else:
            self.sbStride.set_value(self.reader.numDatapoints/self.sbNDataP.get_value())
        self.plot()      
        
    def get_data(self):
        #Zeit stoppen
        #tStart = time.time()
        #Test des lw
        self.chList = []
        #(model,paths) = self.treeV.get_selection().get_selected_rows()
        
        for item in self.tree:
            #print "---",item
            if item[2]:
                #print "   ---",item[0]
                self.chList.append(item[0])
            #print item, model[item], model[item][0], model[item][1]
            #print self.chList
        #print "Welche Kanäle sind gewählt? Dauer: ", time.time()-tStart, "s"
        #tStart = time.time()
        if self.chList == []:
            return False
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        #if self.scale_strideMulti.get()>0:
        #    stride = stride*self.scale_strideMulti.get()
        #print "Start/Laenge/Schritt? Dauer: ", time.time()-tStart, "s"
        tStart=time.time()
        self.f32Lock.acquire()    
        self.data = self.reader.getData(start,length,stride,self.chList)
        self.f32Lock.release()
        #print "Daten holen. Dauer: ", time.time()-tStart, "s"
        #tStart = time.time()
        self.get_times()
        #print "getTimes: Dauer: ", time.time()-tStart, "s"
        #self.ts = [(start+n*stride) for n in range(length)]
        pass
        
    def plot_data(self):
        
        #print self.data
        #print self.data.shape
        self.canvas.hide()
        if self.chList == []:
            print "Keine Kanaele"
            return False
        
        self.a.clear()
        self.subplAxes[0].grid(color='#AAAAAA', linestyle='--', linewidth=0.5)
        
        #Daten umschreiben mit Offset
        #print "Starte umschreiben"
        plotdata = zeros(self.data.shape,"d")
        i=0
        #tStart = time.time()
        for d in range(self.data.shape[1]):
            #plotdata.append([(v-i*self.offset) for v in d])
            plotdata[:,d] = self.data[:,d]-i*self.offset
            i+=1
            #dataObj = []
        #print "Umschreiben fertig"
        #print "plotdata machen. Dauer: ", time.time()-tStart, "s"
        #print plotdata, plotdata.shape
        
        #tStart = time.time()
        #Plotte die Zeitreihe. Überprüfe, ob der Kanal highlighted werden soll
        #(model,paths) = self.treeV.get_selection().get_selected_rows()    
        sgData = []#None  
        plotColor="k"  
        self.a.plot(self.ts,plotdata,plotColor,linewidth=0.5)
        for i in range(plotdata.shape[1]):
            #print len(self.ts), len(plotdata[:,i])
            if self.tree[self.chList[i]][3]:
                self.a.plot(self.ts,plotdata[:,i],"r",linewidth=0.5)
            #    plotColor = "r"
                #if sgData == None:
                #    sgData = plotdata[:,i]
                sgData.append(plotdata[:,i])
            #tPlotStart=time.time()
            #self.a.plot(self.ts,plotdata[:,i],plotColor,linewidth=0.5)
            #print "Einen Kanal plotten. Dauer: ", time.time()-tPlotStart, "s"
            #print i
        #print "Zeitreihen plotten. Dauer: ", time.time()-tStart, "s"
        
        #ytics verändern
        if self.offset > 0:
            ytics = [-(y*self.offset) for y in range(len(self.chList))]
            yticsNames = [self.reader.channel_names[i] for i in self.chList]
            self.subplAxes[0].set_yticks(ytics)
            self.subplAxes[0].set_yticklabels(yticsNames)
            
        else:
            self.subplAxes[0].set_yticks([])
        
        #self.a.set_xlim
        #self.a.autoscale_view(tight=True, scaley=False)
        #self.a.autoscale_view(tight=False,scalex=False, scaley=True)
        #Plotte Marker
        if(self.boolPlotMarkers):
            self.plot_markers()
        
        #xlim/ylim anpassen 
        minx = min(self.ts)
        maxx = max(self.ts)
        #print minx, maxx
        self.subplAxes[0].set_xlim((minx, maxx))
        if self.offset > 0: 
            self.subplAxes[0].set_ylim(0-(len(self.chList)-1)*self.offset-self.offset,self.offset)
        else:
            self.subplAxes[0].set_ylim(-self.offset,self.offset)
        #Show some analysis?
        if self.showAnalysis:
            try:
                self.a3.clear()
            except Exception:
                print "Fehler beim Clearen von a3"
                pass
            if self.whichAnal == "Specgram":
                #print "sg"
                if len(sgData) > 0:
                    #print "sg2"
                    self.a3.specgram(sgData[0],Fs=self.reader.samplRate/self.sbStride.get_value())
                    #Pxx = None
                    #for s in sgData:
                        #if Pxx == None:
                        #    res = specgram(s,Fs=self.reader.samplRate/self.sbStride.get_value())
                        #    Pxx=res[0]
                        #    freqs = res[1]
                        #else:
                        #    Pxx = Pxx + specgram(s,Fs=self.reader.samplRate/self.sbStride.get_value())[0]
                        #print Pxx
                        #print Pxx.shape
                    #im = self.a3.imshow(Pxx,extent=(0,1,freqs[0],freqs[-1]),interpolation="nearest", origin="lower",aspect="auto",cmap=cm.hot)
                    #self.f.colorbar(im)
                    self.a3.set_xticks([])
            elif self.whichAnal == "XCorr":
                #print "xc"
                if len(sgData) >= 2:
                    #print "xc"
                    tmpData = correlate(sgData[0],sgData[1],mode="same")
                    #print tmpData.shape, sgData[0].shape, sgData[1].shape
                    self.a3.plot(tmpData)
            elif self.whichAnal == "Avg":
                #print "xc"
                if len(sgData) >= 2:
                    #print "xc"
                    tmpData = sgData[0]
                    for i in range(1,len(sgData)):
                        tmpData+=sgData[i]
                    tmpData /= len(sgData)
                    self.a3.plot(tmpData)
            elif self.whichAnal == "Wavelet":
                #print "xc"
                if len(sgData) >= 1:
                    tmpData = sgData[0]
                    #wts = wavedec_lin(tmpData,"db4")
                    #self.a3.plot(wts)
                    mld = pywt.wavedec(tmpData, "db4")
                    #print mld
                    level=0
                    self.a3.plot(linspace(0,len(tmpData),mld[-1].shape[0]),mld[-1],label=str(level))
                    self.a3.twinx()
                    for sc in mld[:-1]:
                        self.a3.plot(linspace(0,len(tmpData),sc.shape[0]),sc,label=str(level))
                        level+=1
                    self.a3.legend()
        sys.stdout.flush()
        #self.hpane.set_position(self.hpane.get_position()-1
        self.canvas.show()
        #print "Fertig show"
    

    def get_overview_data(self):
        #self.ovData = self.reader.get
        #print "Starte getOvData"
        #self.f32Lock.acquire()
        #tmpData = self.reader.getOverviewData(100)
        #self.f32Lock.release()
        #self.ovData = tmpData.mean(1)
        #self.ovTs = range(0,self.reader.numDatapoints,self.reader.numDatapoints/tmpData.shape[0])[:100] 
        #self.plotOvData()  
        self.f32Lock.acquire()
        tmpData = self.reader.getOverviewData(1000,range(0,len(self.reader.channel_names),2))
        self.f32Lock.release()
        self.ovData = tmpData.mean(1) 
        self.ovTs = range(0,self.reader.numDatapoints,self.reader.numDatapoints/tmpData.shape[0])[:1000]
        self.setup_subplots()    
        self.plot_overview_data()  
        
    
    def plot_overview_data(self):
        #self.canvas.hide()
        self.a2.clear()
        self.a2.plot(self.ovTs, self.ovData)
        self.subplAxes[1].set_xticks([])
        self.subplAxes[1].set_yticks([])
        self.a2.autoscale_view(tight=True)
        #read scales
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        #set xmin and xmax for rectangle
        xmin = start
        xmax = start+stride*length
        #plot rectangle
        self.a2.axvspan(xmin,xmax,fc="g",alpha=0.4)
        self.canvas.show()
    
    def plot_markers(self):
        """Plottet alle Marker im betrachteten x-Bereich"""
        markersColors = [ '#9E7448' , '#658049' , '#C03DC6' , '#945563' , '#194EB3' , '#718EB6' , '#B03C42' , '#C38591' , '#A2357B' , '#7569A2' , '#7E8DA2' , '#282760' , '#92C08A' , '#372F88' , '#7E3789' , '#768966' , '#A97C21' , '#ACB3C2' , '#1B7266' , '#40569E' , '#BCA4C4' , '#A8996A' , '#584539' , '#9B3F4B' , '#595F58' , '#53798E' , '#6C7183' , '#72AD71' , '#1EA041' , '#7DABBA' , '#24371B' , '#6F7122' , '#9D20B5' , '#593851' , '#9C1A2D' , '#AC57C5' , '#378B48' , '#4223B2' , '#75B32F' , '#B4664E' , '#2E6C85' , '#903EA9' , '#3ABD19' , '#9F6CC5' , '#C01619' , '#6FAD8D' , '#40C753' , '#7F6199' , '#23853D' , '#258816' , '#34518F' , '#4B168D' , '#8CBB95' , '#1B8036' , '#76959F' , '#A6759E' , '#A4C059' , '#498283' , '#A26F5C' , '#50BB35' , '#228796' , '#905A45' , '#AB3694' , '#1B8EBF' , '#92B96F' , '#4D6E4E' , '#701431' , '#5BC547' , '#37495C' , '#486E1B' , '#86449B' , '#B4B68E' , '#958115' , '#20415F' , '#214D4F' , '#2AA42C' , '#5E52B9' , '#178E39' , '#A5166D' , '#A49F79' , '#49C06A' , '#6D881B' , '#701A67' , '#9090AC' , '#377832' , '#51C53B' , '#B53936' , '#334C1D' , '#A147C5' , '#C2C628' , '#363265' , '#633AC5' , '#75AB9C' , '#333079' , '#4B2945' , '#773A57' , '#51B3B2' , '#84A598' , '#274B15' , '#BD4268']
        xrange = self.a.get_xlim()
        #print xrange
        if self.markMaker != None:
            #print "rangeToTake:", [float(x)*self.tsFactor for x in xrange]
            gtkMarkerMarks = self.markMaker.getMarks(rangeToTake=[float(x)*self.tsFactor for x in xrange])
        if self.trigManager != None:
            #print "rangeToTake:", [float(x)*self.tsFactor for x in xrange]
            trigs = self.trigManager.getMarks(rangeToTake=[float(x)*self.tsFactor for x in xrange])
        
        try: 
            textSwitch=0
            for i in trigs.keys():
                for j in trigs[i]:
                    #print (j-self.timesOffset)/self.tsFactor
                    #if (j-self.timesOffset)/self.tsFactor>xrange[0] and (j-self.timesOffset)/self.tsFactor<xrange[1]:
                        #print (j-self.timesOffset)/self.tsFactor
                    try:
                        colorIdx=int(i[-3:])%97
                    except Exception,e:
                        colorIdx=0
                    self.a.axvline((j-self.timesOffset)/self.tsFactor, lw=2, alpha=0.5, color=markersColors[colorIdx])
                    self.a.text((j-self.timesOffset)/self.tsFactor,0+textSwitch*self.offset/2,i,fontsize=6)
                    textSwitch=textSwitch-1
        except Exception, e:
            #print "Keine Markers geladen..."
            pass
        
        try:
            for tp in gtkMarkerMarks:
                self.a.axvline((tp-self.timesOffset)/self.tsFactor, lw=2, alpha=0.5, color="#00937B")
        except Exception, e:
            #print "Keine Markers geladen..."
            pass
#        try:
#            for tp in trigs:
#                self.a.axvline((tp-self.timesOffset)/self.tsFactor, lw=1.5, alpha=0.5, color="#000000")
#        except Exception, e:
#            #print "Keine Markers geladen..."
#            pass
    
    def get_times(self):
        #Zeiten umschreiben
        start=int(self.sbStartpoint.get_value())
        length=int(self.sbNDataP.get_value())
        stride=int(self.sbStride.get_value())
        self.ts = [(start-self.timesOffset+n*stride) for n in range(length)]
        samplR = self.reader.samplRate
        self.tsFactor = 1
        
        if self.rbSecs.get_active(): # Sekunden
            self.ts = [float(t)/samplR for t in self.ts]
            self.tsFactor*=samplR
        elif self.rbMins.get_active(): # Minuten
            self.ts = [float(t)/(samplR*60.) for t in self.ts]
            self.tsFactor*=samplR*60
        elif self.rbHours.get_active(): # Stunden
            self.ts = [float(t)/(samplR*3600.) for t in self.ts]
            self.tsFactor*=samplR*3600
        elif self.rbDays.get_active(): #Tage
            self.ts = [float(t)/(samplR*86400.) for t in self.ts]
            self.tsFactor*=samplR*86400

                               
def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    if len(sys.argv)>1:
        #print sys.argv[1], "wurde als Argument übergeben."
        if os.path.exists(sys.argv[1]):
            p = EEGPlot(sys.argv[1])
        else:
            p = EEGPlot()
    else:
        p = EEGPlot()
    
    main()