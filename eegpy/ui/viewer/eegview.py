# -*- coding: utf-8 -*-

import os
import sys
from threading import Thread
from time import sleep
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, HGroup, \
        HSplit, Handler, ButtonEditor
from enthought.traits.ui.menu import NoButtons, ToolBar, Action
from enthought.traits.ui.key_bindings import KeyBinding, KeyBindings
from mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from scipy import *
import wx

import eegpy
from control import ControlPanel
#from markertab import Marker

#! User interface objects
#!------------------------
#!
#! These objects store information for the program to interact with the
#! user via traitsUI.



#! Threads and flow control
#!--------------------------
#!
#! There are three threads in this application:
#! 
#!  * The GUI event loop, the only thread running at the start of the program.
#!
#!  * The acquisition thread, started through the GUI. This thread is an
#!    infinite loop that waits for the camera to be triggered, retrieves the 
#!    images, displays them, and spawns the processing thread for each image
#!    recieved.
#!
#!  * The processing thread, started by the acquisition thread. This thread
#!    is responsible for the numerical intensive work of the application.
#!    it processes the data and displays the results. It dies when it is done.
#!    One processing thread runs per shot acquired on the camera, but to avoid
#!    accumulation of threads in the case that the processing takes longer than
#!    the time lapse between two images, the acquisition thread checks that the
#!    processing thread is done before spawning a new one.
#! 

def process(image, results_obj):
    """ Function called to do the processing """
    X, Y = indices(image.shape)
    x = sum(X*image)/sum(image)
    y = sum(Y*image)/sum(image)
    width = sqrt(abs(sum(((X-x)**2+(Y-y)**2)*image)/sum(image))) 
    results_obj.x = x
    results_obj.y = y
    results_obj.width = width

#class AcquisitionThread(Thread):
#    """ Acquisition loop. This is the worker thread that retrieves images 
#        from the camera, displays them, and spawns the processing job.
#    """
#    wants_abort = False
#
#    def process(self, image):
#        """ Spawns the processing job.
#        """
#        try:
#            if self.processing_job.isAlive():
#                self.display("Processing to slow")
#                return
#        except AttributeError:
#            pass
#        self.processing_job = Thread(target=process, args=(image,
#                                                            self.results))
#        self.processing_job.start()
#
#    def run(self):
#        """ Runs the acquisition loop.
#        """
#        self.display('Camera started')
#        n_img = 0
#        while not self.wants_abort:
#            n_img += 1
#            img =self.acquire(self.experiment)
#            self.display('%d image captured' % n_img)
#            self.image_show(img)
#            self.process(img)
#            sleep(1)
#        self.display('Camera stopped')

#! The GUI elements
#!------------------
#!

key_bindings = KeyBindings(
    KeyBinding( binding1    = 'Ctrl-left',
                description = 'Move left',
                method_name = 'move_left' ),
    KeyBinding( binding1    = 'Ctrl-right',
                description = 'Move right',
                method_name = 'move_right' ),
    KeyBinding( binding1    = 'Ctrl-1',
                description = 'Mark as sleep-stage 1',
                method_name = 'mark_ss1' ),
    KeyBinding( binding1    = 'Ctrl-2',
                description = 'Mark as sleep-stage 2',
                method_name = 'mark_ss2' ),
    KeyBinding( binding1    = 'Ctrl-3',
                description = 'Mark as sleep-stage 3',
                method_name = 'mark_ss3' ),
    KeyBinding( binding1    = 'Ctrl-4',
                description = 'Mark as sleep-stage 4',
                method_name = 'mark_ss4' ),
    KeyBinding( binding1    = 'Ctrl-r',
                description = 'Mark as REM-sleep',
                method_name = 'mark_REM' ),
    KeyBinding( binding1    = 'Ctrl-a',
                description = 'Mark as awake',
                method_name = 'mark_awake' ),
    KeyBinding( binding1    = 'Ctrl-m',
                description = 'Mark as sleep-stage m',
                method_name = 'mark_movementtime' ),
)


class MainWindowHandler(Handler):
    
    window = Any
    
    def __init__(self,window=None):
        Handler.__init__(self)
        self.window = window
    
    def close(self, info, is_OK):
        #if ( info.object.panel.acquisition_thread 
        #                and info.object.panel.acquisition_thread.isAlive() ):
        #    info.object.panel.acquisition_thread.wants_abort = True
        #    while info.object.panel.acquisition_thread.isAlive():
        #        sleep(0.1)
        #    wx.Yield()
        return True
    
    def move_left ( self, info=None ):
        self.window.panel._move_left_fired()
        
    def move_right ( self, info=None ):
        self.window.panel._move_right_fired()
    
    def mark_ss1 ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="Stage 1")
    
    def mark_ss2 ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="Stage 2")
    
    def mark_ss3 ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="Stage 3")
    
    def mark_ss4 ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="Stage 4")
    
    def mark_REM ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="REM")
    
    def mark_awake ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="awake")
    
    def mark_movementtime ( self, info=None ):
        self.window.panel.marker.append(self.window.panel.inputs.start,name="MT")


class MainWindow(HasTraits):
    """ The main window, here go the instructions to create and destroy
        the application.
    """
    figure = Instance(Figure)
    panel = Instance(ControlPanel)
    mpl_editor = MPLFigureEditor()
    
    initial_fn = Str("/media/Extern/public")
    
#    def __init__(self,fn=None):
#        HasTraits.__init__(self)
#        if fn!=None:
#            self.initial_fn = fn
        
    
    def _figure_default(self):
        figure = Figure()
        return figure

    def _panel_default(self):
        panel = ControlPanel(figure=self.figure, initial_fn=self.initial_fn)
        #print "connect"
        #panel.figure.canvas.mpl_connect('button_press_event', panel._click_cb)
        #print "connected"
        panel._do_plot_fired()
        return panel
    
    view = View(HSplit(  Item('figure',  editor=mpl_editor,
                                                        dock='vertical'),
                        Item('panel', style="custom", width=100),
                    show_labels=False, 
                    ),
                resizable=True, 
                height=1.0, width=1.0,
                #handler=my_handler,
                buttons=NoButtons,
                
                key_bindings=key_bindings,
                title="eegpy viewer 1.0a",
                )
    
                

if __name__ == '__main__':
    try:
        import psyco
        psyco.full()
        print "Speeding things up with psyco"
    except ImportError:
        pass
    if len(sys.argv)>1:
        #print sys.argv[1], "wurde als Argument übergeben."
        if os.path.exists(sys.argv[1]):
            win = MainWindow(initial_fn = sys.argv[1])
        else:
            win = MainWindow()
    else:
        win = MainWindow()
    hand = MainWindowHandler(win)
    win.configure_traits(handler=hand)
    
