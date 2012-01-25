#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################
# Module-Import #
#################

#eegpy-modules
try:
    import eegpy
    from eegpy.events import EventTable
    from eegpy.misc import FATALERROR
    from eegpy.ui.widgets.windowwidgets import EegpyBaseWin
    from eegpy.ui.icon import image_from_eegpy_stock, eegpy_logo
except ImportError:
    raise FATALERROR('Your installation of EegPy seems to be incomplete.\nMaybe you need to set the PYTHONPATH environment-variable adequatly.')

#from eegpy.filter.filt_misc import filterRecursively

#Third-party
try:
    import numpy
    from scipy.signal import lfilter, butter
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pygtk
    pygtk.require('2.0')
    import gobject
    import gtk
except ImportError:
    raise FATALERROR('GTK cannot be imported.')

#try:
#    from matplotlib.axes import Subplot
#    # uncomment to select /GTK/GTKAgg/GTKCairo
#    from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
#    from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
#    import matplotlib
#    #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, NavigationToolbar
#    from matplotlib.figure import Figure, SubplotParams
#    from matplotlib.axis import Axis
#    import matplotlib.cm
#except ImportError:
#    raise FATALERROR('Error while importing matplotib. Please visit http://matplotlib.sf.net for more information.')

#native python
import sys
import os
import pickle

class EventManager(gtk.Frame):
    _et = None
    _fn = None
    _keylist = None
    
    def __init__(self, label=""):
        gtk.Frame.__init__(self,label)
        self.vbox=gtk.VBox()
        self.tb_box = gtk.HBox()
        self.add(self.vbox)
        self.vbox.pack_start(self.tb_box,expand=False)
        
        self.tb = gtk.Toolbar()
        self.tooltips = gtk.Tooltips()
        self.tb.set_style(gtk.TOOLBAR_ICONS)
        self.add_toolbutton_from_stock(gtk.STOCK_OPEN, 'Load', 'Load an EventTable from a file', 'Private', self.load_et)
        self.add_toolbutton_from_stock(gtk.STOCK_SAVE, 'Save', 'Save the EventTable back to the original file', 'Private', self.save_et, False)
        self.add_toolbutton_from_stock(gtk.STOCK_SAVE_AS, 'Save to', 'Save the EventTable to a file, choose new file', 'Private', self.save_et, True)
        self.tb.insert(gtk.SeparatorToolItem(),-1)
        self.add_toolbutton_eegpy("add_trigger_type", "Add type", "Add a new trigger type", 'Private', self.cb_add_trigger_type, None)
        self.add_toolbutton_eegpy("add_trigger", "Add trigger", "Add a new trigger", 'Private', self.cb_add_trigger, None)
        self.tb_box.pack_start(self.tb,expand=True)
        self.lb_fn = gtk.Label("New EventTable...")
        self.lb_fn.set_max_width_chars(50)
        self.lb_fn.set_justify(gtk.JUSTIFY_RIGHT)
        self.tb_box.pack_end(self.lb_fn, expand=False)
        #HBox für _keylist/triggerlist
        self.pane_kl = gtk.HPaned()
        self.vbox.pack_end(self.pane_kl)
        self.setup_trees()
        self._et = EventTable()
   
    def setup_trees(self):   
        #First: Keys
        self.tvsw_keys = gtk.ScrolledWindow()
        self.tvsw_keys.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.tree_keys = gtk.TreeStore(gobject.TYPE_STRING)
        #self.treeS = gtk.TreeModelSort(self.tree)
        self.tv_keys = gtk.TreeView(self.tree_keys)
        self.tv_keys.get_selection().connect("changed",self.key_selected)
        #self.tv_keys.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        #renderer = gtk.CellRendererText()
        #self.col1 = gtk.TreeViewColumn("File ...", renderer,text=0)
        self.tv_keys.append_column(gtk.TreeViewColumn("Key", gtk.CellRendererText(),text=0))
        #self.tv_keys.show()
        self.tvsw_keys.add(self.tv_keys)
        self.pane_kl.add1(self.tvsw_keys)
        #Second: Triggers
        self.tvsw_tr = gtk.ScrolledWindow()
        self.tvsw_tr.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.tree_tr = gtk.TreeStore(gobject.TYPE_INT)
        #self.treeS = gtk.TreeModelSort(self.tree)
        self.tv_tr = gtk.TreeView(self.tree_tr)
        self.tv_tr.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        #renderer = gtk.CellRendererText()
        #self.col1 = gtk.TreeViewColumn("File ...", renderer,text=0)
        self.tv_tr.append_column(gtk.TreeViewColumn("Timepoint", gtk.CellRendererText(),text=0))
        #self.tv_keys.show()
        #Setting up drag'n'drop
        self.tv_tr.enable_model_drag_source( gtk.gdk.BUTTON1_MASK,
                                                [('INT',0,0)],
                                                gtk.gdk.ACTION_DEFAULT|
                                                gtk.gdk.ACTION_MOVE)
        self.tv_tr.enable_model_drag_dest([('INT',0,0)],
                                             gtk.gdk.ACTION_DEFAULT)

        self.tv_tr.connect("drag_data_get", self.tr_drag_get)
        self.tv_tr.connect("drag_data_received", self.tr_drag_received)
        self.tv_keys.connect("key_press_event", self.cb_key_pressed)
        self.tv_tr.connect("key_press_event", self.cb_key_pressed)

        self.tvsw_tr.add(self.tv_tr)
        self.pane_kl.add2(self.tvsw_tr)
    
    def add_toolbutton_eegpy(self, icon_name, text, tip_text, tip_private, clicked_function, clicked_param1=None):
        iconSize = gtk.ICON_SIZE_SMALL_TOOLBAR
        iconw = eegpy.ui.icon.image_from_eegpy_stock(icon_name)
            
        toolitem = gtk.ToolButton(iconw, text)
        #toolitem = gtk.ToolButton(iconw)
        toolitem.set_icon_widget(iconw)
        toolitem.show_all()
        toolitem.set_tooltip(self.tooltips, tip_text, tip_private)
        toolitem.connect("clicked", clicked_function, clicked_param1)
        #toolitem.connect("scroll_event", clicked_function)
        self.tb.insert(toolitem, -1)
                
    def add_toolbutton_from_stock(self, icon_name, text, tip_text, tip_private, clicked_function, clicked_param1=None):
        iconSize = gtk.ICON_SIZE_SMALL_TOOLBAR
        iconw = gtk.Image()
        iconw.set_from_stock(icon_name, iconSize)
            
        toolitem = gtk.ToolButton(iconw, text)
        #toolitem = gtk.ToolButton(iconw)
        toolitem.set_icon_widget(iconw)
        toolitem.show_all()
        toolitem.set_tooltip(self.tooltips, tip_text, tip_private)
        toolitem.connect("clicked", clicked_function, clicked_param1)
        #toolitem.connect("scroll_event", clicked_function)
        self.tb.insert(toolitem, -1)
    
    def load_et(self,event,data):
        dialog = gtk.FileChooserDialog("Open EventTable from file..", None, gtk.FILE_CHOOSER_ACTION_OPEN, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("eegpy EventTable or similar")
        filter.add_pattern("*.evt")
        filter.add_pattern("*.vmrk")
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
    
    def save_et(self, event, do_save_as = True):
        if do_save_as == False:
            self._et.save(self._fn)
        else:
            dialog = gtk.FileChooserDialog("Save EventTable to file...", None, gtk.FILE_CHOOSER_ACTION_SAVE, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
            dialog.set_default_response(gtk.RESPONSE_OK)
            
            filter = gtk.FileFilter()
            filter.set_name("eegpy EventTable")
            filter.add_pattern("*.evt")
            dialog.add_filter(filter)
            
            filter = gtk.FileFilter()
            filter.set_name("All files")
            filter.add_pattern("*")
            dialog.add_filter(filter)
            
            response = dialog.run()
            if response == gtk.RESPONSE_OK:
                fn = dialog.get_filename()
                print fn, 'selected'
                dialog.destroy()
                self._fn = fn
                #Now save...
                self._et.save(self._fn)
                lbtext = ""
                if len(fn)>40:
                    lbtext = "..."+fn[-37:]
                self.lb_fn.set_text(lbtext) 
                #fh.close()
            else:# response == gtk.RESPONSE_CANCEL:
                dialog.destroy()
                print 'Closed, no files selected'
            pass
            
            

    def set_filename(self,fn):
        print fn, "selected for opening"
        #success = False
        try:
            if not os.path.exists(fn):
                raise ValueError("File doesn't exist")
            self._et = EventTable(fn)
            if len(self._et.keys())==0:
                print self._et.keys()
                raise ValueError("EventTable empty!")
            self._fn = fn
        except ValueError, e:
            print "Error opening EventTable", e
            self._et=None
            self._fn=None
            return False
        lbtext = ""
        if len(fn)>40:
            lbtext = "..."+fn[-37:]
        self.lb_fn.set_text(lbtext)
        self.setup_keylist()
    
    def setup_keylist(self):
        #if self._tv!=None:
        #    try:
        #        self._keylist.hide()
        #        self._keylist.destroy()
        #    except Exception,e:
        #        print "Cannot destroy keylist"
        #TODO: Real functionalityself.tvsw_keys = gtk.ScrolledWindow()
        
        keys = self._et.keys()
        keys.sort()
        self.tree_keys.clear()
        for k in keys:
            iter = self.tree_keys.append(None)
            self.tree_keys.set(iter, 0, k)          
        self.tree_keys.set_sort_column_id(0,gtk.SORT_ASCENDING)
        self.show_all()
    
    def setup_triggerlist(self, key):
        self.tree_tr.clear()
        for tr in self._et[key]:
            #print tr
            iter = self.tree_tr.append(None)
            self.tree_tr.set(iter, 0, int(tr))          
        self.tree_tr.set_sort_column_id(0,gtk.SORT_ASCENDING) 
        
    def key_selected(self,treeselection,*args):
        #print tv, path, col, args, self.tree_keys.get(self.tree_keys.get_iter(path),0)[0]
        self.tv_tr.get_selection().unselect_all()
        #self.tree_tr.clear()
        paths = treeselection.get_selected_rows()[1]
        if len(paths)>0:
            iter = self.tree_keys.get_iter(paths[0])
            key = self.tree_keys.get(iter,0)[0]
            self.setup_triggerlist(key)
        
    
    def cb_add_trigger_type(self,event,data):
        dialog_label = gtk.Dialog("Choose name...", None, gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT, (gtk.STOCK_CANCEL, gtk.RESPONSE_REJECT, gtk.STOCK_OK, gtk.RESPONSE_OK))    
        entry1 = gtk.Entry()
        entry1.set_text("Trigger")
        dialog_label.vbox.pack_start(entry1)
        entry1.show()
        response = dialog_label.run()
        print response
        if response == gtk.RESPONSE_OK:
            trig_name = entry1.get_text()
            print trig_name
        else:
            print "Adding trigger-type aborted by user."
            dialog_label.destroy()
            return False    
        dialog_label.destroy()  
        self.add_trigger_type(trig_name, [])
    
    def cb_add_trigger(self,event,data):
        dialog_label = gtk.Dialog("Add trigger...", None, gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT, (gtk.STOCK_CANCEL, gtk.RESPONSE_REJECT, gtk.STOCK_OK, gtk.RESPONSE_OK))    
        dialog_label.vbox.pack_start(gtk.Label("Timepoint:"))
        sb_time = gtk.SpinButton(gtk.Adjustment(0,0,100000000,1,1000))
        dialog_label.vbox.pack_start(sb_time)
        dialog_label.vbox.show_all()
        response = dialog_label.run()
        print response
        if response == gtk.RESPONSE_OK:
            time = sb_time.get_value()
            print time
        else:
            print "Adding trigger aborted by user."
            dialog_label.destroy()
            return False    
        dialog_label.destroy()  
        self.add_trigger(time)
    
    def add_trigger_type(self,key,ts=[]):
        if not self._et.has_key(key):
            self._et.add_trigger_type(key, ts)
        self.setup_keylist()
        self.tree_tr.clear()
    
    def add_trigger(self,time):
        #find out key
        path = self.tv_keys.get_selection().get_selected_rows()[1][0]
        iter = self.tree_keys.get_iter(path)
        k = self.tree_keys.get(iter,0)[0]
        if self._et.has_key(k):
            self._et.add_trigger(k, time)
        self.setup_triggerlist(k)
    
    def tr_drag_get(self, treeview, context, selection, target_id, etime): 
        pathlist = treeview.get_selection().get_selected_rows()[1]
        model = treeview.get_model()
        iterlist = [model.get_iter(row) for row in pathlist]
        datalist = [model.get(iter,0)[0] for iter in iterlist]
        #print datalist
        selection.set(selection.target,8,pickle.dumps(datalist))
        #print "Drag_get: ", treeview, context, selection, target_id, etime
        
    
    def tr_drag_received(self, treeview, context, x, y, selection, info, etime):
        #print pickle.loads(selection.data)
        datalist = pickle.loads(selection.data)
        self.add_trigger(datalist[0])
        #print "Drag_received:", treeview, context, x, y, selection, info, etime
    
    def cb_key_pressed(self, widget, event, data=None):
        keyname = gtk.gdk.keyval_name(event.keyval)
        #print "Key %s (%d) was pressed in widget %s" % (keyname, event.keyval, str(widget))
        if keyname == "Delete":
            #find out key
            path = self.tv_keys.get_selection().get_selected_rows()[1][0]
            iter = self.tree_keys.get_iter(path)
            k = self.tree_keys.get(iter,0)[0]
            if widget==self.tv_keys:
                self._et.remove(k)
                self.setup_keylist()
                self.tv_keys.get_selection().unselect_all()
                self.tree_tr.clear()
            if widget==self.tv_tr:
                pathlist = self.tv_tr.get_selection().get_selected_rows()[1]
                iterlist = [self.tree_tr.get_iter(row) for row in pathlist]
                datalist = [self.tree_tr.get(iter,0)[0] for iter in iterlist]
                for tr in datalist:
                    self._et.remove(k,tr)
                self.setup_triggerlist(k)
    

class EventTableEditorWin(EegpyBaseWin):
    programName = "eegpy: Frequency-Filtering"
    
    # Konstruktor
    def __init__(self):
        EegpyBaseWin.__init__(self)
        self.inner_pane.set_position(300)
        self.em1 = EventManager("EventTable 1")
        self.em1.tv_tr.get_selection().connect("changed",self.cb_plot_marks)#, "blue")
        self.em2 = EventManager("EventTable 2")
        self.em2.tv_tr.get_selection().connect("changed",self.cb_plot_marks)#, "red")
        self.pane_edit = gtk.HPaned()
        self.upper_hbox.pack_start(self.pane_edit)
        self.pane_edit.add1(self.em1)
        self.pane_edit.pack2(self.em2,False)
        self.pane_edit.set_position(self.get_size()[0]/2)
        #self.setupOptions()
        self.show_all()
        #self.setupGUI()
    
    def setupGUI(self):
        EegpyBaseWin.setupGUI(self)
         
    def cb_plot_marks(self, treeselection, *args):
        #print "Color", color
        self.a.cla()
        pathlist = self.em1.tv_tr.get_selection().get_selected_rows()[1]
        iterlist = [self.em1.tree_tr.get_iter(row) for row in pathlist]
        datalist1 = [self.em1.tree_tr.get(iter,0)[0] for iter in iterlist]
        pathlist = self.em2.tv_tr.get_selection().get_selected_rows()[1]
        iterlist = [self.em2.tree_tr.get_iter(row) for row in pathlist]
        datalist2 = [self.em2.tree_tr.get(iter,0)[0] for iter in iterlist]
        #print datalist1, datalist2
        for i in datalist1:
        #    print i,
            self.a.axvline(i, lw=1, color="blue", ymin=0.5, ymax=1)
        #self.a.plot(datalist1,numpy.zeros(len(datalist1)),"bD")
        #self.a.plot(datalist2,numpy.ones(len(datalist2)),"rD")
        #print ""
        for i in datalist2:
        #    print i,
            self.a.axvline(i, lw=1, color="red", ymin=0, ymax=0.5)
        #print ""
#        if len(datalist1) == 1:
#            self.a.set_xlim(datalist1[0]-1000,datalist1[0]+1000)
#        elif len(datalist2)==1:
#            self.a.set_xlim(datalist2[0]-1000,datalist2[0]+1000)
#        else:
#            self.a.autoscale_view()
#        elif:
#            xlim0 = max(min(datalist1),min(datalist2))-500
#            xlim1 = min(max(datalist1),max(datalist2))+500
#            if xlim1<xlim0:
#                xlim0 = min(min(datalist1),min(datalist2))-500
#                xlim1 = max(max(datalist1),max(datalist2))+500
#            self.a.set_xlim(xlim0,xlim1)
        #self.a.set_xlim(numpy.array(datalist1+datalist2).min()-1000,numpy.array(datalist1+datalist2).max()+1000)
        self.a.set_ylim(0,1)
        self.a.set_yticks([])
        self.canvas.draw()
        

        

def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    etew = EventTableEditorWin()
    main()