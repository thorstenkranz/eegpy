# -*- coding: utf-8 -*-

from time import sleep
import os.path
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, HGroup, VGroup, \
        HSplit, Handler, ButtonEditor, ListEditor, SetEditor, TabularEditor
from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, ToolBar, Action
from enthought.pyface.api import ImageResource, FileDialog, YES, CANCEL
from enthought.traits.ui.tabular_adapter import TabularAdapter
from matplotlib.figure import Figure
from scipy import *
import wx
import eegpy
from message_box import Message
from etmanager import EventTableManager

import eegpy
from eegpy import EventTable
from mpl_figure_editor import MPLFigureEditor
from message_box import Message

class Marker(HasTraits):
    t = Int
    name = Str

class MarkerlistAdapter ( TabularAdapter ):

    columns = [ ( 'Time',    't' ), 
                ( 'Name',     'name' )]
                
    font                      = 'Courier 10'
    #age_alignment             = Constant( 'right' )
    #MarriedPerson_age_image   = Property
    #MarriedPerson_bg_color    = Color( 0xE0E0FF )
    #MarriedPerson_spouse_text = Property
    #Person_spouse_text        = Constant( '' )
    
    #def _get_MarriedPerson_age_image ( self ):
    #    if self.item.age < 18:
    #        return 'red_flag'
    #    return None
    #    
    #def _get_MarriedPerson_spouse_text ( self ):
    #    return self.item.partner.name

#tabular_editor = TabularEditor(
#    adapter    = MarkerlistAdapter(),
#    operations = [ 'delete' ],
#    multi_select = True,
#    activated = "goto_marker",
#    dclicked = "goto_marker",
#    editable=False
    #images     = [ ImageResource( 'red_flag', search_path = search_path ) ]
#)

class MarkerTab(HasTraits):
    """ Object used to display the results.
    """
    
    record_mark = Bool(False)
    name_new = Str("Mark")
    markers = List( Marker ) 
    update_marks = Int(0)
    active_row = Int(0)
    
    remove_f2 = Button("Remove")
    load_evt = Action(name = "Load",
                        action = "_load",
                        toolip = "Load markers from EventTable",
                        image = ImageResource("images/load_32.png")
                        )
    save_evt = Action(name = "EventTable",
                        action = "_save_to",
                        toolip = "Save markers as EventTable",
                        image = ImageResource("images/save_32.png")
                        )
    save_ascii = Action(name = "ASCII",
                        action = "_save_to_ascii",
                        toolip = "Save markers as ASCII-file",
                        image = ImageResource("images/save_32.png")
                        )
    toolbar = ToolBar(load_evt, save_evt, save_ascii)
    
    def __init__(self):
        HasTraits.__init__(self)
        #self.markers.append(Marker(t=10,name="Test10"))
        #self.markers.append(Marker(t=100,name="Test100"))
    
    #x = Float(50, label="X", desc="X position of the center")
    #y = Float(50, label="Y", desc="Y position of the center")
    
    def cmp_markers(self,x,y):
        return x.t-y.t
    
    def append(self,t,name=None,notify=True,do_sort=True):
        if name==None:
            name=str(self.name_new)
        self.markers.append(Marker(t=t,name=name))
        if do_sort:
            self.markers.sort(cmp=self.cmp_markers)
        if notify:
            self.update_marks+=1
        #self.goto_marker = self.markers[-1]
        
    def get_marks(self,start,stop):
        rv = []
        for m in self.markers:
            if m.t>=start and m.t<=stop:
                rv.append(m)
        return rv
    
    def _load(self):
        print "Load from EventTable"
        extension = "evt"
        wildcard = "EventTable|*.evt|VA MarkerFile|*.vmrk"
        fileDialog = FileDialog(action='open', title='Load EventTable',
                                     wildcard=wildcard)
        fileDialog.open()
        if fileDialog.path == '' or fileDialog.return_code == CANCEL:
            return False
        else:
            print "Opening", fileDialog.path
            #et = eegpy.load(str(fileDialog.path))
            et = eegpy.EventTable(str(fileDialog.path))
            for k in et.keys():
                for t in et[k]:
                    #print k,t
                    self.append(t,k,False,False)
            self.markers.sort(cmp=self.cmp_markers)
            
            
    def _save_to(self):
        print "Save to EventTable"
        extension = "evt"
        wildcard = "*.evt"
        fileDialog = FileDialog(action='save as', title='Save As',
                                     wildcard=wildcard)
        fileDialog.open()
        if fileDialog.path == '' or fileDialog.return_code == CANCEL:
            return False
        else:
            extLen = len(extension)
            if extLen and fileDialog.path[-extLen-1:] != '.' + extension:
                fileDialog.path += '.' + extension               
        #print "fc.fn:", fileDialog.path
        #TODO: Check if file exists and join EventTables
        et = EventTable()
        for m in self.markers:
            try:
                et.add_trigger(m.name,m.t)
            except ValueError, e:
                et.add_trigger_type(m.name, [m.t])
        try:
            et.save(fileDialog.path) 
        except Exception, e:
            mb = Message(message="Error while writing EventTable to file.\n%s"%str(e))
        return True
    
    def _save_to_ascii(self):
        print "Save to ASCII"
        pass
    
    @on_trait_change("active_row")
    def go_to(self):
        print "Going to mark", self.markers[self.active_row].name, "at", self.markers[self.active_row].t
    
        
    traits_view = View( Group(
                   VGroup(
                    Item("record_mark",show_label=True,label = "Record marks:"),
                    Item("name_new",show_label=True,label = "Name for new marks:"),
                   ),
                   Item( 'markers', id = 'table', 
                         editor = TabularEditor(
                                        adapter    = MarkerlistAdapter(),
                                        operations = [ 'delete' ],
                                        multi_select = True,
                                        activated = "goto_marker",
                                        activated_row = "active_row",
                                        dclicked = "goto_marker",
                                        editable=False
                                        #images     = [ ImageResource( 'red_flag', search_path = search_path ) ]
                                  ) ), 
                       show_labels = False
                   ),
                 toolbar=toolbar,
               )
    
        

    

        
        
            
        
        

