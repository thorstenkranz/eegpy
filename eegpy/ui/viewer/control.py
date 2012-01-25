# -*- coding: utf-8 -*-
#! The imports
#!-------------
#!
#! The MPLFigureEditor is imported from last example.

from threading import Thread
from time import sleep
import os.path
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, HGroup, VGroup, \
        HSplit, Handler, ButtonEditor, ListEditor, CheckListEditor, spring
from enthought.traits.ui.key_bindings import KeyBinding, KeyBindings
from enthought.traits.ui.menu import NoButtons, ToolBar, Action
from enthought.pyface.api import ImageResource
from mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
from scipy import *
from scipy.signal import hilbert, detrend
import numpy as N
import wx
import eegpy
from files import FileSelection
from markertab import MarkerTab
from formula_eval import FormulaTab
from message_box import Message

from eegpy.tools import schlaf_ged
from eegpy.filter.freqfilt import filtfilt_low,filtfilt_high,filtfilt_band,filtfilt_bandstop

#color_seq = ["k","r","b","c","m","y"]

input_presets = {
    "Oddball_eeg":[None,1100,5,None],
    "20 s":[None,1000,20,None],
    "10 s":[None,1000,10,None],
    "1 s":[None,1000,1,None],
}

class Inputs(HasTraits):
    """ Object that contains the parameters that control the experiment, 
        modified by the user.
    """
    start = Int(0, label="Start", desc="Startpoint from where to plot")
    length = Int(1100, label="N/o samples.", desc="Number of datapoints to get")
    step = Int(5, label="Step", desc="Step size between the samples")
    offset = Int(50, label="Scaling", desc="Vertical distance between channels")
    preset = Enum("Oddball_eeg","20 s","10 s","1 s")
    timescale = Enum("samples","seconds","minutes","hours")
    
    linewidth = Range(0.1,2.0,0.5)
    color1 = RGBColor((0,0,0))
    color2 = RGBColor((0,0,1))
    xtickswidth = Int(1000, label="X-Ticks", desc="Horizontal distance between ticks on x-axis")

    
    show_trial_length = Bool(False, label="Show")
    trial_0 = Int(-500, label="t0", desc="Beginning of trial")
    trial_1 = Int(1000, label="t0", desc="End of trial")
    trial_alpha = Enum(0.3,0.1,0.2,0.3,0.4,0.5,0.8,1)
    #color_trial = RGBColor((0,0,1))
    
    traits_view = View(VGroup(
                       VGroup(
                            Item('preset'),
                            Item('start'),
                            Item('length'),
                            Item('step'),
                            Item('offset'),
                            Item('timescale'),
                       ),
                       VGroup(
                            Item('linewidth'),
                            Item('color1', style='custom', springy=True,resizable=False,show_label=True,label="Color 1"),
                            Item('color2', style='custom', springy=True,resizable=False,show_label=True,label="Color 2"),
                            Item('xtickswidth',label="Grid-spacing"),
                            VGroup(
                                 Item('show_trial_length'),
                                Item('trial_0'),
                                Item('trial_1'),
                                Item('trial_alpha'),
                                #Item('color_trial', style='custom', springy=True,resizable=False,show_label=True,label="Color for Trial"),
                                label="Plot trial-length around markers",
                            ),
                             label="Plot-parameters",
                        ),
                 ),
            ) 
    
    def _preset_changed(self):
        slso =  input_presets[str(self.preset)]
        if slso[0] != None:
            self.start = slso[0]
        if slso[1] != None:
            self.length = slso[1]
        if slso[2] != None:
            self.step = slso[2]
        if slso[3] != None:
            self.offset = slso[3]



        
       
class ControlTab(HasTraits):
    """ This object is the core of the traitsUI interface. Its view is
        the right panel of the application, and it hosts the method for
        interaction between the objects and the GUI.
    """
    move_step = Enum(1,0.2,0.5,1,2,5,10)
    scale_step = Enum(10,1,5,10,20,50,100)
    channels = List()
    ch_editor = CheckListEditor(values=[],cols=2)
    num_channels = Int(0)
    patient_name = Str()
    select_all = Button("all")
    select_none = Button("none")
    select_good = Button("good (db)")
    mlcount = Int(0)
    mrcount = Int(0)
    mucount = Int(0)
    mdcount = Int(0)
    #_sls = [0,1100,5] #StartLengthStride
    #acquisition_thread = Instance(AcquisitionThread)
    
    move_left = Action(name = "Left",
                        action = "_trigger_move_left",
                        toolip = "Move to the left",
                        image = ImageResource("images/left_32.png")
                        )
    move_right = Action(name = "Right",
                        action = "_trigger_move_right",
                        toolip = "Move to the right",
                        image = ImageResource("images/right_32.png")
                        )
    scale_up = Action(name = "Scale up",
                        action = "_trigger_scale_up",
                        toolip = "Move to the right",
                        image = ImageResource("images/up_32.png")
                        )
    scale_down = Action(name = "Scale down",
                        action = "_trigger_scale_down",
                        toolip = "Move to the right",
                        image = ImageResource("images/down_32.png")
                        )
                        
    toolbar = ToolBar(move_left, move_right,scale_up,scale_down)
    
    
    traits_view = View(VGroup(
                  #Item('do_plot', show_label=False ),
                    VGroup(
                      Item('move_step',label="Move by:"),
                      Item('scale_step',label="Scale by"),
                    ),
                    VGroup(
                      Group(
                          Item('channels',
                               show_label=False,
                               style='custom',
                               editor=ch_editor,
                               #springy=True,
                               #height=-500,
                          ),
                          scrollable=True,
                          springy=True,
                      ),
                      HGroup(
                        spring,
                        Item('select_all',show_label=False),
                        Item('select_none',show_label=False),
                        Item('select_good',show_label=False),
                        spring,
                      ),
                      label = "Choose channels",
                      #height=300,
                    ),
                  #Item('results_string',show_label=False, 
                  #                      springy=True, style='custom' ),
                  #springy=True,
                  #scrollable=True,
                  #height=-500,
                  dock='tab',
                  springy=True
                  #scrollable=True
                  ),
               toolbar=toolbar,
               scrollable=True
               #height=-500,
               )
               #key_bindings=key_bindings,
               #handler=ControlHandler,
               
    def _trigger_move_left(self):
        #print "Trigger Left"
        self.mlcount+=1
    def _trigger_move_right(self):
        self.mrcount+=1
        #print "Trigger Right"
    def _trigger_scale_up(self):
        self.mucount+=1
    def _trigger_scale_down(self):
        self.mdcount+=1

    def _select_none_fired(self):
        self.channels = []
 
    def _select_all_fired(self):
        self.channels = range(self.num_channels)

    def _select_good_fired(self):
        try:
            gc = schlaf_ged.get_good_channels(self.patient_name)
            print "Got", gc, "from database"
            self.channels = gc
        except ValueError, ve:
            #Zeige Fehlermeldung
            print ve
            print self.patient_name
            try:
                short_patient_name = self.patient_name.split("_")[0]
                gc = schlaf_ged.get_good_channels(short_patient_name)
                print "Got", gc, "from database"
                self.channels = gc
            except ValueError,ve:
                print ve
                print short_patient_name
                message = Message(message="Cannot find good channels for subject %s" % self.patient_name)
                message.edit_traits()
        except Exception, e:
            #Zeige Fehlermeldung
            print e
            pass
        
class ControlPanel(HasTraits):
    """ This object is the core of the traitsUI interface. Its view is
        the right panel of the application, and it hosts the method for
        interaction between the objects and the GUI.
    """
    #fn = ""
    reader1 = None#eegpy.F32(fn,"r")
    reader2 = None
    _ts_factor = 1.
    controltab = Instance(ControlTab, ())
    inputs = Instance(Inputs, ())
    #camera = Instance(Camera, ())
    figure = Instance(Figure)
    file_selection = Instance(FileSelection, ())
    marker = Instance(MarkerTab,())
    formulastab = Instance(FormulaTab,())
    #start_stop_acquisition = Button("Start/Stop acquisition")
    do_plot = Button("Plot")
    #move_left = Button("<")
    #move_right = Button(">")
    scale_up = Button("up")
    scale_down = Button("down")
    results_string =  String()
    ovData = None
    subplots = None
    
    #Visibility of tabs
    show_markerstab = Bool(True)
    show_formulastab = Bool(True)
    show_filestab = Bool(True)
    show_settingstab = Bool(True) 
    
    initial_fn = Str("")
    
    lines = {}
    mpl_connected = Bool(False)
    
    update_count = 0 #Every time lines are updated instead of plotted, increased by one. 
    #Triggers clean redraw if > 10
    #_sls = [0,1100,5] #StartLengthStride
    #acquisition_thread = Instance(AcquisitionThread)
    
    #For formula evaluation
    chs_for_formulas = None
    formula_labels = None
    formulas = []
    data_for_formulas = None
    
    view = View(Group(
                    Item('controltab', style='custom', show_label=False,label="Control",resizable=True,springy=True),
                    #),
                    Group(
                      Group(
                        Item('inputs', style='custom', show_label=False),
                        label="Input",),
                    label='Settings', dock="tab",visible_when="show_settingstab"),
                    Group(
                     Item('file_selection', style='custom', show_label=False,),
                     dock="tab", visible_when='show_filestab',label="Files"
                    ),
                    Group(
                    Item('marker', style='custom', show_label=False,),
                    dock="tab", visible_when='show_markerstab',label="Markers"
                    ),
                    Group(
                    Item('formulastab', style='custom', show_label=False,),
                    dock="tab", visible_when='show_formulastab',label="Formulas"
                    ),
                    Group(
                      Group(
                        Item('show_settingstab', show_label=True),
                        Item('show_filestab', show_label=True),
                        Item('show_markerstab', show_label=True),
                        Item('show_formulastab', show_label=True),
                        label="Which tabs should be shown?",),
                    label='View', dock="tab"), 
                   layout='tabbed'),
                )
                   #key_bindings=key_bindings,
                   #handler=ControlHandler,
    
    #@on_trait_change("inputs.[start,length,step,offset]")
    #def do_test(self):
    #    print "Test2"
        
    def _file_selection_default(self):
        file_selection = FileSelection(fn=self.initial_fn)
        #file_selection.on_trait_change(self.update)
        return file_selection
    
    @on_trait_change("file_selection.[long_eeg_filename1,long_eeg_filename2]")    
    def update(self):
        #Create readers that do not exist
        #reader_exists = N.zeros((len(self.file_selection.long_eeg_filenames)),N.bool)
        print "update", self.file_selection.long_eeg_filename1
        if len(self.file_selection.long_eeg_filename1) > 0:
            if self.reader1 == None or (self.reader1 != None and self.reader1.fn != self.file_selection.long_eeg_filename1):
                self.reader1 = eegpy.F32(self.file_selection.long_eeg_filename1,"r")
                #self.channels = []
                #for ch in self.reader1.channel_names:
                self.controltab.ch_editor.values = [(i,cn) for i,cn in enumerate(self.reader1.channel_names)]
                self.controltab.channels = range(self.reader1.num_channels)
                self.controltab.num_channels = self.reader1.num_channels
                self.controltab.patient_name = self.file_selection.eeg_filename1[:-4]
                self.inputs.start=0
                self.formulas=[]
        else:
            self.reader1=None
            self.subplots = ""
            self.lines = {}
            self.figure.clf()
            wx.CallAfter(self.figure.canvas.draw)
        if len(self.file_selection.long_eeg_filename2) > 0:
            if self.reader2 == None or (self.reader1 != None and self.reader2.fn != self.file_selection.long_eeg_filename2):
                self.reader2 = eegpy.F32(self.file_selection.long_eeg_filename2,"r")
            #Check if dimesions are equal
            if self.reader2 != None:
                if not self.reader2.shape==self.reader1.shape:
                    self.reader2 = None
                    self.file_selection.long_eeg_filename2 = ""
                    self.file_selection.eeg_filename2 = ""
                    dialog = Message(message="The selected file doesn't have the same dimensions")
                    dialog.edit_traits()
        else:
            self.reader2=None
        
        if not self.mpl_connected:
            self.figure.canvas.mpl_connect('button_press_event', self.canvas_cb)
            self.mpl_connected = True
        self._do_plot_fired()
    
    @on_trait_change("inputs.[color1,color2,linewidth]")
    def force_redraw(self):
        self.update_count=10e50
        self._do_plot_fired()
    
    @on_trait_change("[inputs.[start,length,step,offset,show_trial_length,trial_0,trial_1,trial_alpha],marker.update_marks]")
    def _do_plot_fired(self):
        """ Plot with the given values
        """
        #print "_do_plot_fired"
        if self.reader1 == None:
            return False
        
        #print "_do_plot_fired2"
        
        #Control settings in self.inputs
        start,length,step = self.inputs.start, self.inputs.length, self.inputs.step
        if start<0:
            start=0
        if start+length>self.reader1.num_datapoints:
            start = self.reader1.num_datapoints-length
        if length*step>self.reader1.num_datapoints:
            step = self.reader1.num_datapoints/length
        if start+step*length>self.reader1.numDatapoints:
            start=self.reader1.numDatapoints-step*length
        self.inputs.start, self.inputs.length, self.inputs.step = start, length, step
        if not self.inputs.offset>0:
            self.inputs.offset=self.controltab.scale_step
        #Dann DAten holen
        #tStart = time.time()
        #color_seq = self.inputs.color_seq
        #try:
        #    print self.subplots
        #except Exception,e:
        #    self._setup_subplots(self.figure)
        #if self.subplots == None:
        self._setup_subplots(self.figure)
        #self.figure.axes[0].clear()
        self.figure.axes[0].cla()
        if self.reader1 != None:
            #print "Plot for reader1"
            self._get_data(self.reader1)
            #Dann Daten plotten
            self._plot_data("main",self.inputs.color1) 
            #print "Color 1:", self.file_selection.color1
        if self.reader2 != None:
            self._get_data(self.reader2)
            #Dann Daten plotten
            self._plot_data("compare",self.inputs.color2) 
        self.update_count+=1
        self.plot_markers()
        self.decorate_plot()
    
    @on_trait_change("controltab.mlcount")    
    def _move_left_fired(self):
        """Move view to the left
        """
        self.inputs.start-=int(self.inputs.length*self.inputs.step*self.controltab.move_step)
        #self._do_plot_fired()
    
    @on_trait_change("controltab.mrcount")    
    def _move_right_fired(self):
        """Move view to the left
        """
        self.inputs.start+=int(self.inputs.length*self.inputs.step*self.controltab.move_step)
        #self._do_plot_fired()
        
    @on_trait_change("marker.active_row")
    def _goto_mark(self):
        """React to changes of activated row in the markers tab
        """
        self.inputs.start = self.marker.markers[self.marker.active_row].t #- 50
    
    @on_trait_change("controltab.mucount")    
    def _scale_up_fired(self):
        """Increase the scale
        """
        self.inputs.offset-=int(self.controltab.scale_step)
        self._do_plot_fired()
    
    @on_trait_change("controltab.mdcount")    
    def _scale_down_fired(self):
        """Increase the scale
        """
        self.inputs.offset+=int(self.controltab.scale_step)
        self._do_plot_fired()
        
    @on_trait_change("inputs.timescale")    
    def _set_ts_factor(self):
        try:
            Fs= self.reader1.Fs
            if str(self.inputs.timescale)=="samples":
                self._ts_factor=1.
            elif str(self.inputs.timescale)=="seconds":
                self._ts_factor=1./Fs
            elif str(self.inputs.timescale)=="minutes":
                self._ts_factor=1./60/Fs
            elif str(self.inputs.timescale)=="hours":
                self._ts_factor=1./3600/Fs
        except Exception, e:
            print type(e), e
            pass

    def _get_data(self, reader):
        #Zeit stoppen
        #tStart = time.time()
        #Test des lw
        #tStart = time.time()
        self.chList = self.controltab.channels#range(self.reader1.num_channels)
        self.chList.sort()
        if self.chList == []:
            return False
        start=self.inputs.start
        length=self.inputs.length
        step=self.inputs.step
        #if self.scale_strideMulti.get()>0:
        #    stride = stride*self.scale_strideMulti.get()
        #print "Start/Laenge/Schritt? Dauer: ", time.time()-tStart, "s"
        #tStart=time.time()    
        self.data = reader.getData(start,length,step,self.chList)

        #Data for formulas
        #print self.formulas
        #if len(self.formulas>0):
        try:
            self.data_for_formulas = reader.getData(start,length,step,self.chs_for_formulas)
        except:
            self.data_for_formulas = None

        #self.get_times()
        self.ts = [(start+n*step) for n in range(length)]
        pass
        
    def _plot_data(self,id,color):
        if self.chList == []:
            print "Keine Kanaele"
            return False
        
        offset = self.inputs.offset
        
        #self.figure.axes[0].clear()
        self.figure.axes[0].grid(True,which="minor",color='#AAAAAA', linestyle='--', linewidth=0.5)
        
        #Daten umschreiben mit Offset
        #print "Starte umschreiben"
        plotdata = zeros((self.data.shape[0],self.data.shape[1]+len(self.formulas)),"d")
        i=0
        #tStart = time.time()
        for d in range(self.data.shape[1]):
            #plotdata.append([(v-i*self.offset) for v in d])
            plotdata[:,d] = self.data[:,d]-i*offset
            i+=1
        for d in range(len(self.formulas)):
            try:
                formula_res_ar = eval(self.formulas[d])
                #print formula_res_ar, formula_res_ar.shape
                plotdata[:,d+self.data.shape[1]] = formula_res_ar-i*offset
            except Exception, e:
                self.formulastab.formulas[d].error = str(e)
                plotdata[:,d+self.data.shape[1]] = -i*offset
            i+=1
            
        
        
        #try: 
        #    for i in range(plotdata.shape[1]):
                #print self.lines[i],
                #print self.lines[i].get_xdata().mean(),
                #print self.lines[i].get_ydata().mean(),
        #        self.lines[id][i].set_xdata(self.ts)
        #        self.lines[id][i].set_ydata(plotdata[:,i])
        #        self.lines[id][i]._invalid = True
                #print self.lines[i].get_xdata().mean(),
                #print self.lines[i].get_ydata().mean()
            #print "Update", i
        #except Exception,e:
            #print "plot", i, e   
        self.lines[id] = self.figure.axes[0].plot(self.ts,plotdata,color=color,linewidth=self.inputs.linewidth,antialiased=True)
        major_xticks = N.arange(self.ts[0],self.ts[-1],((self.ts[-1]-self.ts[0])/5/100)*100)
        major_xticklabels = ["+%i"%tx for tx in major_xticks-self.ts[0]]
        major_xticklabels[0] = "%.3f %s"%(self.ts[0]*self._ts_factor,str(self.inputs.timescale))
        self.figure.axes[0].set_xticks(major_xticks,minor=False)
        self.figure.axes[0].set_xticklabels(major_xticklabels,minor=False)
        self.figure.axes[0].set_xticks(N.arange(self.ts[0],self.ts[-1],self.inputs.xtickswidth),minor=True)
            
    def plot_markers(self):
        """Plots all Markers in the viewed part of the eeg"""
        markersColors = [ '#9E7448' , '#658049' , '#C03DC6' , '#945563' , '#194EB3' , '#718EB6' , '#B03C42' , '#C38591' , '#A2357B' , '#7569A2' , '#7E8DA2' , '#282760' , '#92C08A' , '#372F88' , '#7E3789' , '#768966' , '#A97C21' , '#ACB3C2' , '#1B7266' , '#40569E' , '#BCA4C4' , '#A8996A' , '#584539' , '#9B3F4B' , '#595F58' , '#53798E' , '#6C7183' , '#72AD71' , '#1EA041' , '#7DABBA' , '#24371B' , '#6F7122' , '#9D20B5' , '#593851' , '#9C1A2D' , '#AC57C5' , '#378B48' , '#4223B2' , '#75B32F' , '#B4664E' , '#2E6C85' , '#903EA9' , '#3ABD19' , '#9F6CC5' , '#C01619' , '#6FAD8D' , '#40C753' , '#7F6199' , '#23853D' , '#258816' , '#34518F' , '#4B168D' , '#8CBB95' , '#1B8036' , '#76959F' , '#A6759E' , '#A4C059' , '#498283' , '#A26F5C' , '#50BB35' , '#228796' , '#905A45' , '#AB3694' , '#1B8EBF' , '#92B96F' , '#4D6E4E' , '#701431' , '#5BC547' , '#37495C' , '#486E1B' , '#86449B' , '#B4B68E' , '#958115' , '#20415F' , '#214D4F' , '#2AA42C' , '#5E52B9' , '#178E39' , '#A5166D' , '#A49F79' , '#49C06A' , '#6D881B' , '#701A67' , '#9090AC' , '#377832' , '#51C53B' , '#B53936' , '#334C1D' , '#A147C5' , '#C2C628' , '#363265' , '#633AC5' , '#75AB9C' , '#333079' , '#4B2945' , '#773A57' , '#51B3B2' , '#84A598' , '#274B15' , '#BD4268']
        xrange = self.figure.axes[0].get_xlim()
        #print xrange
        
        markers = self.marker.get_marks(xrange[0],xrange[1])#self.inputs.start,self.inputs.start+self.inputs.length*self.inputs.step)
        ylim_ax = self.figure.axes[0].get_ylim()
        y_diff = ylim_ax[1]-ylim_ax[0]
        try:
            y_pos=0 # ylim_ax[0]+0.8*y_diff
            for m in markers:
                #TODO: change when other timeunits are supported
                self.figure.axes[0].axvline(m.t, lw=2, alpha=0.5, color=markersColors[hash(m.name)%len(markersColors)])
                #text_y = self.figure.axes[0].transData.inverted().transform((self.figure.axes[0].transAxes.transform([0,textSwitch])))[1]
                #print "text_y", text_y, self.figure.axes[0].get_ylim()
                self.figure.axes[0].text(m.t,y_pos,m.name,fontsize=10,bbox=dict(facecolor="white",edgecolor="red",alpha=0.5))
                if self.inputs.show_trial_length:
                    self.figure.axes[0].axvspan(m.t+self.inputs.trial_0,m.t+self.inputs.trial_1,alpha=self.inputs.trial_alpha,color=markersColors[hash(m.name)%len(markersColors)])
                #textSwitch-=0.025
                #if textSwitch<0:
                #    textSwitch=0.8
                y_pos-=0.02*y_diff
                if y_pos<(ylim_ax[0]+0.1*y_diff):
                    y_pos=0 #ylim_ax[0]+0.8*y_diff
        except Exception, e:
            #print "Keine Markers geladen..."
            print e
            pass
         
        evt_markers = self.file_selection.etm.get_marks(xrange[0],xrange[1])#self.inputs.start,self.inputs.start+self.inputs.length*self.inputs.step)
        try:
            for m in evt_markers:
                
                #TODO: change when other timeunits are supported
                self.figure.axes[0].axvline(m[1], lw=2, alpha=0.5, color=markersColors[hash(m[0])%len(markersColors)])
                #text_y = self.figure.axes[0].transData.inverted().transform((self.figure.axes[0].transAxes.transform([0,textSwitch])))[1]
                #print "text_y", text_y
                self.figure.axes[0].text(m[1],y_pos,m[0],fontsize=10,bbox=dict(facecolor="white",edgecolor="blue",alpha=0.5))
                if self.inputs.show_trial_length:
                    self.figure.axes[0].axvspan(m[1]+self.inputs.trial_0,m[1]+self.inputs.trial_1,alpha=self.inputs.trial_alpha,color=markersColors[hash(m[0])%len(markersColors)])
                #textSwitch-=0.025
                #if textSwitch<0:
                #    textSwitch=0.8
                y_pos-=0.02*y_diff
                if y_pos<(ylim_ax[0]+0.1*y_diff):
                    y_pos=0 #ylim_ax[0]+0.8*y_diff
        except Exception, e:
            #print "Keine Markers geladen..."
            print e
            pass
        
        #if self.trigManager != None:
            #print "rangeToTake:", [float(x)*self.tsFactor for x in xrange]
        #    trigs = self.trigManager.getMarks(rangeToTake=[float(x)*self.tsFactor for x in xrange])
        
        #try: 
        #    textSwitch=0
        #    for i in trigs.keys():
        #        for j in trigs[i]:
                    #print (j-self.timesOffset)/self.tsFactor
                    #if (j-self.timesOffset)/self.tsFactor>xrange[0] and (j-self.timesOffset)/self.tsFactor<xrange[1]:
                        #print (j-self.timesOffset)/self.tsFactor
        #            try:
        #                colorIdx=int(i[-3:])%97
        #            except Exception,e:
        #                colorIdx=0
        #            self.a.axvline((j-self.timesOffset)/self.tsFactor, lw=2, alpha=0.5, color=markersColors[colorIdx])
        #            self.a.text((j-self.timesOffset)/self.tsFactor,0+textSwitch*self.offset/2,i,fontsize=6)
        #            textSwitch=textSwitch-1
        #except Exception, e:
            #print "Keine Markers geladen..."
        #    pass
        
               
        
    def decorate_plot(self):
        #ytics verändern
        offset = self.inputs.offset
        if  offset > 0:
            ytics = [-(y*offset) for y in range(len(self.chList)+len(self.formulas))]
            yticsNames = [self.reader1.channel_names[i] for i in self.chList] + [f.name for f in self.formulastab.formulas]
            self.figure.axes[0].set_yticks(ytics)
            self.figure.axes[0].set_yticklabels(yticsNames)
        else:
            self.figure.axes[0].set_yticks([])
        
        
        #xlim/ylim anpassen 
        minx = min(self.ts)
        maxx = max(self.ts)
        #print minx, maxx
        self.figure.axes[0].set_xlim((minx, maxx))
        if offset > 0: 
            self.figure.axes[0].set_ylim(0-(len(self.chList)+len(self.formulas)-1)*offset-offset,offset)
        else:
            self.figure.axes[0].set_ylim(-offset,offset)
        
        self._plot_overview_data()
        #sys.stdout.flush()
        wx.CallAfter(self.figure.canvas.draw)
        #print "Fertig show"  
    
    def _plot_overview_data(self):
        #self.canvas.hide()
        if self.ovData == None:
            tmpData = self.reader1.getOverviewData(1000,range(0,len(self.reader1.channel_names),2))
            self.ovData = tmpData.mean(axis=1) 
            self.ovTs = range(0,self.reader1.numDatapoints,self.reader1.numDatapoints/tmpData.shape[0])[:1000]
        
        self.figure.axes[1].clear()
        self.figure.axes[1].plot(self.ovTs, self.ovData)
        self.figure.axes[1].set_xticks([])
        self.figure.axes[1].set_yticks([])
        self.figure.axes[1].autoscale_view(tight=True)
        #set xmin and xmax for rectangle
        xmin = self.inputs.start
        xmax = self.inputs.start+self.inputs.step*self.inputs.length
        #plot rectangle
        self.figure.axes[1].axvspan(xmin,xmax,fc="g",alpha=0.4)
    
    def _setup_subplots(self,figure):
        #self.figure.clear()
#        if False:#self.showAnalysis:
#            self.a = self.f.add_subplot(212)
#            self.a2 = self.f.add_subplot(414)
#            self.a3 = self.f.add_subplot(211)
#            self.subplAxes = self.f.get_axes()
#            self.subplAxes[0].set_position([0.06,0.05,0.91,0.50])
#            self.subplAxes[2].set_position([0.06,0.55,0.91,0.35])
#            self.subplAxes[1].set_position([0.06,0.96,0.91,0.02])
#            self.subplAxes[1].set_xticks([])
#            self.subplAxes[1].set_yticks([])
#        else:
        if not self.subplots == "plot":
            self.subplots = "plot"
            figure.clf()
            figure.add_axes([0.09, 0.05, 0.9, 0.9])
            figure.add_axes([0.09,0.96,0.9,0.02])
            figure.axes[1].set_xticks([])
            figure.axes[1].set_yticks([])
            
    def canvas_cb(self,event=None):
        #print "event:", event
        if event.inaxes == self.figure.axes[1]:
            #print "In Übersicht"            
            self.inputs.start = int(event.xdata)
            self._do_plot_fired()
        elif event.inaxes == self.figure.axes[0]:
            #print "In Übersicht"            
            if self.marker.record_mark:
                self.marker.append(int(event.xdata))
            self._do_plot_fired()
        
    @on_trait_change("[formulastab.formulas,formulastab.update_formulas]")
    def _create_formula_chs(self):
        """ Create pseudo-channels to plot the results of the formulas
        """
        if self.reader1 == None:
            return False
        print "formulas changed"
        self.chs_for_formulas, cnames_for_formulas = self.formulastab.get_all_channels(self.reader1.channel_names)
        self.formula_labels = [f.name for f in self.formulastab.formulas]
        self.formulas = []
        for formula in self.formulastab.formulas:
            formula = formula.formula
            for icn, cname in enumerate(cnames_for_formulas):
                #print icn,cname,self.chs_for_formulas,cnames_for_formulas
                formula = formula.replace(cname,"self.data_for_formulas[:,%i]"%icn)
            self.formulas.append(formula)
        self._do_plot_fired()


                
    @on_trait_change("[show_markerstab,show_filestab]")
    def redraw_window(self):
        #TODO: implement redraw of window
        print "redraw should be called"
        
        
