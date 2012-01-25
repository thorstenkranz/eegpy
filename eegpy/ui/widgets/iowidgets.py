#!/usr/bin/env python
# -*- coding: utf-8 -*-

import eegpy
from eegpy.misc import FATALERROR

#Third-party
try:
    import numpy
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

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
import time

#se = gtk.stock_lookup(gtk.STOCK_EXECUTE)
#se = list(x for x in se)
#print se
#se[1] = " "
#se = tuple(se)
#print se
#gtk.stock_add([tuple(se)])

class VarChooser(gtk.Frame):
    """Widget for choosing an input/output-variable."""
    onActivate = None
    table = None
    entry = None
    var_ = None
    
    def __init__(self, label="", io="in", var_=123456, onActivate = None):
        assert type(label) == type(""), "varChooser(): Label is not a string"
        assert io in ["in","out"], "varChooser(): The io-parameter has an illegal value: "+io+". Choose either \"in\" or \"out\"."
        
        #print label
        #sys.stdout.flush()
        self.var_ = var_
        
        gtk.Frame.__init__(self,label)
        self.hbox = gtk.HBox(False,0)
        self.add(self.hbox)
        
        if io=="in":
            if type(self.var_) == type(123456):
                if self.var_ == 123456:
                    self.label = gtk.Label("Loading from files not implemented yet.")
                    self.hbox.pack_start(self.label, expand=False, fill=True)
            elif type(self.var_) == type(numpy.zeros((1,1))):     
                self.label = gtk.Label("Array of shape\n"+str(self.var_.shape)+"\nloaded")
                self.hbox.pack_start(self.label, expand=False, fill=True) 
                
                self.button = gtk.Button(stock=gtk.STOCK_EXECUTE)
                #self.button.set_label("")
                #print self.button.get_label()
                self.button.connect("clicked",self.cb_onActivate)      
                self.hbox.pack_start(self.button, expand=False, fill=True)
            else:
                raise TypeError, "varChooser(): variable-type is not adequate."
        else: #io=="out"
            if type(self.var_) == type(123456):
                if self.var_ == 123456:
                    self.label = gtk.Label("Writing to files not implemented yet.")
                    self.hbox.pack_start(self.label, expand=False, fill=True)
            else:     
                self.label = gtk.Label("Called interactively.\nResults will be returned.")
                self.hbox.pack_start(self.label, expand=False, fill=True) 
                
                self.button = gtk.Button(stock=gtk.STOCK_SAVE) 
                self.button.connect("clicked",self.cb_onActivate)      
                self.hbox.pack_start(self.button, expand=False, fill=True)

        self.onActivate = onActivate
    
    def cb_onActivate(self,widget):
        if not self.onActivate == None:
            try:
                self.onActivate(self)
            except TypeError, e:
                print """varChooser: Calling the defined callback failed.
It must be a method of the signature callback(widget), where widget will be set to this varChooser-instance.
The error was: %s""", e
        
    def get_text(self):
        return self.entry.get_text()
    
    def get_var(self):
        if self.var_ == None:
            "Loading from files not implemented yet."
        elif type(self.var_) == type(numpy.zeros((1,1))):
            return self.var_
        return None
    
    def set_var(self, x):
        self.var_ = x
        #print self.var_
        #self.label.set_text("")


class VarNameCompleterEntry(gtk.Entry):
    """An enhanced Entry-Field for auto-completion of variable-names."""
    KEY_TAB = 65289
    def __init__(self, completion_getter=None):
        gtk.Entry.__init__(self)

        if not completion_getter == None:
            self.completion_getter = completion_getter
        else:
            print "Standard"
            self.completion_getter = self.get_completion

        self.completion = gtk.EntryCompletion()
        self.completion.set_model(None) # EntryCompletion is overzealous.
        # Only give it a model when we've explicitly asked for completion.
        self.completion.set_inline_selection(True)
        self.set_completion(self.completion)
        self.completion.set_minimum_key_length(1)
        self.completion.set_text_column(0)

        self.connect('key-press-event', self.entry_keypress_cb)
        self.completion.connect('match-selected', self.match_cb)

    def entry_keypress_cb(self, widget, event):
        if event.keyval == self.KEY_TAB:
            liststore = self.completion_getter(self.get_text())

            if len(liststore) == 1:
                self.set_text(liststore[0][0])
                self.set_position(-1)
            else:
                self.completion.set_model(liststore)
                self.completion.complete()
                gtk.main_do_event(gtk.gdk.Event(gtk.gdk.KEY_PRESS))
            return True
        else:
            self.completion.set_model(None)
            return False

    def match_cb(self, completion, model, iter):
        print model[iter][0], 'was selected'
        completion.set_model(None)
        return

    def get_completion(self, prefix):
        """This method selects those global variables which start with the string in the Entry-field"""
        liststore = gtk.ListStore(str)
        for k,v in globals().items():
            print k
            if k.startswith(prefix):
                liststore.append([k])
        for k,v in locals().items():
            print k
            if k.startswith(prefix):
                liststore.append([k])
                #print s
        return liststore
    
class ComponentChooser(gtk.Window):
    """Used for selecting components to use for remixing"""
    
    components = None # Contains the components to choose from
    data = None # Contains the original data
    choose_mask = []
    
    def __init__(self, _components = None, _data=None):
        self.components = _components
        self.data = _data
        
        gtk.Window.__init__(self,gtk.WINDOW_TOPLEVEL)
        self.set_default_size(700,500)
        self.set_title("Choose components to use...")
        self.connect("delete_event", self.delete_event)
        self.connect("destroy", self.delete_event)
        
        self.vbox = gtk.VBox()
        self.add(self.vbox)
        self.f = Figure(figsize=(5,4), dpi=100, subplotpars=SubplotParams(left=0.06, top=0.95, right=0.97, bottom=0.1))
        self.setupSubplots()
        self.plotAll()
        self.canvas = FigureCanvas(self.f)
        self.canvas.mpl_connect('button_press_event', self.cb_canvas)
        self.vbox.pack_start(self.canvas)
        self.toolbar = NavigationToolbar( self.canvas, self.window )
        self.vbox.pack_start(self.toolbar, False, False)
        self.show_all()
        
    def delete_event(self, widget, event, data=None):
        self.hide()
        return True
    
    def cb_canvas(self, event):
        """Called when the figure is clicked.""" 
        #print event.button
        if event.button == 3:   
            for i,a in enumerate(self.axs):
                if event.inaxes == a:
                    #print "Subplot %i geclickt" % i
                    self.choose_mask[i] = bool(True-self.choose_mask[i])
                    #print self.choose_mask[i]
                    #a.grid()
                    a.cla()
                    #TODO: plotte...
                    a.plot(self.components[:,i])
                    #
                    if not self.choose_mask[i]:
                        a.text(numpy.array(a.get_xlim()).mean(),numpy.array(a.get_ylim()).mean(),"X",horizontalalignment='center',verticalalignment='center',fontsize=20,color="red")
                    self.canvas.draw()
                    
        
    def setupSubplots(self):
        self.f.clear()
        try:
            ncols = 4
            nrows = self.components.shape[1]/ncols +1
            self.axs = []
            self.choose_mask = []
            for i in range(self.components.shape[1]):
                self.axs.append(self.f.add_subplot(nrows,ncols,i+1))
                self.choose_mask.append(True)
        except Exception, e:
            print "Error while setting up subplots", e
    
    def plotAll(self):
        for i,a in enumerate(self.axs):
            a.plot(self.components[:,i]) 
    
    def set_components(self):
        pass
    
    def set_data(self):
        """"""
        pass
    