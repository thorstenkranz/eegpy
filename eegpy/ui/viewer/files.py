# -*- coding: utf-8 -*-

from time import sleep
import os.path
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, HGroup, \
        HSplit, Handler, ButtonEditor, ListEditor, SetEditor, ColorEditor
from enthought.traits.ui.menu import NoButtons, OKButton, ToolBar, Action
from enthought.traits.ui.ui_traits import Image
from enthought.pyface.api import ImageResource
from mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from scipy import *
import wx
import eegpy
from message_box import Message
from etmanager import EventTableManager


class FileSelectionHandler(Handler):    
    def object_fn_changed(self,info):
        print "fn has changed,", info.object.fn
        #print self.controler
        if os.path.exists(info.object.fn):
            if os.path.splitext(info.object.fn)[1] in [".f32",".dat"]:
                if len(info.object.eeg_filename1)==0:
                    info.object.eeg_filename1 = os.path.split(info.object.fn)[1]
                    info.object.long_eeg_filename1 = info.object.fn
                elif len(info.object.eeg_filename2)==0:
                    info.object.eeg_filename2 = os.path.split(info.object.fn)[1]
                    info.object.long_eeg_filename2 = info.object.fn
                else:
                    dialog = Message(message="Maximum number of eeg-files reached.")
                    dialog.edit_traits()
                    return False
                info.object.fn = os.path.split(os.path.abspath(info.object.fn))[0]
            elif os.path.splitext(info.object.fn)[1] in [".evt",".vmrk"]:
                #info.object.evt_filenames.append(os.path.split(info.object.fn)[1])
                #info.object.long_evt_filenames.append(info.object.fn)
                info.object.etm.append(info.object.fn)
            elif os.path.isdir(info.object.fn):
                pass
            else:
                dialog = Message(message="Unknown file extension!")
                dialog.edit_traits()
                return False
            #info.object.controler.update()
        else:
            dialog = Message(message="File not found!")
            print "fn:", info.object.fn
            dialog.edit_traits()
            return False
        
    def object_remove_f1_changed(self,info):
        print "FileSelection._remove_f1_fired"
        info.object.eeg_filename1 = ""
        info.object.long_eeg_filename1 = ""
        
    def object_remove_f2_changed(self,info):
        print "FileSelection._remove_f2_fired"
        info.object.eeg_filename2 = ""
        info.object.long_eeg_filename2 = ""


class FileSelection(HasTraits):
    """ Object used to display the results.
    """
    
    fn = File()
    eeg_filename1 = Str
    long_eeg_filename1 = Str
    remove_f1 = Button("Remove")
    eeg_filename2 = Str
    long_eeg_filename2 = Str
    remove_f2 = Button("Remove")
    
    etm = Instance(EventTableManager,())
    
    def _fn_default(self):
        tmp_fn="/media/Extern/public"
        if not os.path.exists(tmp_fn):
            tmp_fn="."
        fn = File(tmp_fn, label="File", desc="Select filename") 
        return fn
        

    traits_view = View( Item('fn', 
                      style='custom',
                      show_label=False,
                 ),
                 Group(
                    HGroup(
                        Item('eeg_filename1', style="readonly", springy=True, label="Main file:"),
                        Item('remove_f1',style="custom",image=ImageResource("images/up.png"),width=-100,resizable=False,show_label=False),
                    ),
                    HGroup(
                        Item('eeg_filename2', style="readonly", springy=True, label="Compare to:"),
                        Item('remove_f2',width=-100,resizable=False,show_label=False),
                    ),
                    label="EEG-files",
                 ),
                 Item('etm',style="custom",show_label=False
                 ),
                 handler=FileSelectionHandler(),
               )
    

        
        
            
        
        

