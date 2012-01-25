# -*- coding: utf-8 -*-

from time import sleep
import os.path
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, HGroup, VGroup, \
        HSplit, Handler, ButtonEditor, ListEditor, SetEditor, TabularEditor
from enthought.traits.ui.menu import NoButtons, OKButton, ToolBar, Action
from enthought.traits.ui.tabular_adapter import TabularAdapter
from mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from scipy import *
import wx
import eegpy
from message_box import Message

class EventTableListEntry(HasTraits):
    fn = Str
    short_fn = Str
    marker_count = Int
    et = Instance( eegpy.EventTable )
    
    def _fn_changed(self):
        print "EventTableListEntry._fn_changed"
        print self.fn
        self.short_fn = os.path.split(self.fn)[1]
        self.et = eegpy.EventTable(str(self.fn))
        self.marker_count = len(self.et.get_all_events())
    
class EventTableAdapter ( TabularAdapter ):

    columns = [ ( 'Filename',    'short_fn' ), 
                ( '# of events',     'marker_count' )]
                
    font = 'Courier 10'

tabular_editor = TabularEditor(
    adapter    = EventTableAdapter(),
    operations = [ ],#'delete' ],
    multi_select = True,
    #images     = [ ImageResource( 'red_flag', search_path = search_path ) ]
)

class EventTableManager(HasTraits):
    """ Manage the EventTables used for the display. 
    Also controls where new marks are saved to.
    """
    evt_filenames = List(EventTableListEntry)
    evts = List(Instance(eegpy.EventTable))

    view = View( VGroup(
                    Group(
                    Item( 'evt_filenames', editor = tabular_editor,  
                       show_label=False,
                    ),
                    label="EventTables",
                 ),
                 ),
               )
    
    def append(self,et_fn):
        for etle in self.evt_filenames:
            if etle.fn == et_fn:
                return False
        self.evt_filenames.append(EventTableListEntry(fn=et_fn))
        #self.evt_filenames.sort(cmp=lambda x,y: cmp(x.short_fn,y.short_fn))
        self.evts.append(eegpy.EventTable(str(et_fn)))
        
    def get_marks(self,start,stop):
        rv = []
        for et in self.evts:
            for k in et.keys():
                for t in et[k]:
                    if t>=start and t<=stop:
                        rv.append((k,t))
        return rv
    
  

        
        
            
        
        

