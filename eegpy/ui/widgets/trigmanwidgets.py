# -*- coding: utf-8 -*-

import os.path
import sys, time
import pickle

import pygtk
pygtk.require('2.0')
import gobject
import gtk

import numpy as n

from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.axes import Subplot
from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
from matplotlib.figure import Figure, SubplotParams
from matplotlib.axis import Axis

from eegpy.formats import f32
from Eeg.Formats import bvFormat
from eegpy.events import EventTable
#import pylab




class TriggerManager:
    """Managing the display of Triggers from files."""
       
    def __init__(self, plot=None):
        self.etDict = {}
        self.plot = plot
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_default_size(200,300)
        self.window.set_title("Trigger-Manager")

        # Set a handler for delete_event that immediately
        # exits GTK.
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.delete_event)
        
        self.mainBox = gtk.VBox()
        self.window.add(self.mainBox)
        
        self.tvScrolledWin = gtk.ScrolledWindow()
        self.tvScrolledWin.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.tree = gtk.TreeStore(gobject.TYPE_STRING,gobject.TYPE_INT)
        #self.treeS = gtk.TreeModelSort(self.tree)
        self.treeV = gtk.TreeView(self.tree)
        self.treeV.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        renderer = gtk.CellRendererText()
        self.col1 = gtk.TreeViewColumn("File ...", renderer,text=0)
        self.treeV.append_column(self.col1)
        self.col2 = gtk.TreeViewColumn("Offset", renderer,text=1)
        self.treeV.append_column(self.col2)
        self.treeV.show()
        self.tvScrolledWin.add(self.treeV)
        self.tvScrolledWin.show_all()
        #self.hbox.pack_start(self.tvScrolledWin)
        self.mainBox.pack_start(self.tvScrolledWin)
        
        #self.cbx_pauseMarking = gtk.CheckButton("pause marking")
        #self.mainBox.pack_start(self.cbx_pauseMarking,expand=False)
        
        self.btAdd = gtk.Button("Add Trigger-file")
        self.btAdd.connect("clicked", self.cb_add)
        self.mainBox.pack_start(self.btAdd,expand=False)
        self.btRemove = gtk.Button("Remove marked")
        self.btRemove.connect("clicked", self.cb_remove)
        self.mainBox.pack_start(self.btRemove,expand=False)
        #self.btSave = gtk.Button("Save to pickle-file")
        #self.btRemove.set_size_request(-1, 20)
        #self.btSave.connect("clicked", self.cb_save)
        #self.mainBox.pack_start(self.btSave,expand=False)
        
        self.window.show_all()
    
    def delete_event(self, widget, event=None, data=None):
        self.window.hide()
        return True
    
    def cb_add(self, widget):
        dialog = gtk.FileChooserDialog("Open trigger-File..", None, gtk.FILE_CHOOSER_ACTION_OPEN, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("EEG-events (eegpy,BVA)")
        filter.add_pattern("*.vmrk")
        filter.add_pattern("*.evt")
        dialog.add_filter(filter)
        
        filter = gtk.FileFilter()
        filter.set_name("All files")
        filter.add_pattern("*")
        dialog.add_filter(filter)
        
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            self.add(dialog.get_filename())
            print dialog.get_filename(), 'selected'
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()
        pass
    
    def cb_remove(self, widget):
        def remove(model, path, iter):
            model.remove(iter)
            return False     # keep the foreach going
        
        pathlist = self.treeV.get_selection().get_selected_rows()[1]
        iterlist = [self.tree.get_iter(row) for row in pathlist]
        for row in iterlist:
            try:
                print self.tree.get(row,0)[0]
                del self.etDict[self.tree.get(row,0)[0]]
            except KeyError, e:
                print "KeyError", e, self.etDict
            self.tree.remove(row)
        self.plot.plot()
            #print row
    
    def cb_save(self, widget):
        marklist = []
        def append(model, path, iter, user_data):
            marklist.append(self.tree.get(iter,0)[0])
            
        dialog = gtk.FileChooserDialog("Save list as pickle-file ...", None, gtk.FILE_CHOOSER_ACTION_SAVE, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("Pickle-file")
        filter.add_pattern("*.pkl")
        dialog.add_filter(filter)
        
        filter = gtk.FileFilter()
        filter.set_name("All files")
        filter.add_pattern("*")
        dialog.add_filter(filter)
        
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            self.tree.foreach(append, "")
            marklist.sort()
            print "Marker:", marklist
            fh = open(dialog.get_filename(),"w")
            pickle.dump(marklist,fh,-1)
            print dialog.get_filename(), 'selected'
            fh.close()
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()
    
    def getMarks(self,rangeToTake=None):
        marklist = {}
        for k in self.etDict.keys():
            for k2 in self.etDict[k].keys():
                #print k, k2
                if not marklist.has_key(k2):
                    marklist[k2] = []
                for xVal in self.etDict[k][k2]:
                    if rangeToTake==None:
                        marklist[k2].append(xVal)
                    else:
                        assert len(rangeToTake)>1, "rangeToTake must be of length 2"
                        if xVal>rangeToTake[0] and xVal<rangeToTake[1]:
                            marklist[k2].append(xVal)
                marklist[k2].sort()
        #print marklist
        return marklist
    
    def add(self, fn):
        iter = self.tree.append(None)
        self.tree.set(iter, 0, os.path.split(fn)[1])  
        self.tree.set(iter, 1, 0)  
        self.tree.set_sort_column_id(0,gtk.SORT_ASCENDING)
        self.etDict[os.path.split(fn)[1]] = EventTable(fn)
        print self.etDict
        for k in self.etDict.keys():
            print "Keys: ", self.etDict[k].keys()
        
             
        