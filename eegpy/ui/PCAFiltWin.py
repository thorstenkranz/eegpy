#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################
# Module-Import #
#################

#eegpy-modules
try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.filter import pcafilt, icafilt
    #from eegpy.filter.pcafilt import #
    from eegpy.ui.widgets.iowidgets import VarChooser, ComponentChooser
    from eegpy.ui.widgets.plotwidgets import MultichannelSubplot
    from eegpy.ui.widgets.windowwidgets import EegpyWin
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

#native python
import sys

class PCAFiltWin(EegpyWin):
    programName = "eegpy: PCA-Filtering"
    x=None
    x2 = None #For PCA
    y=None
    compCh = None
    
    # Konstruktor
    def __init__(self, x_in=123456, x_out=123456):
        gobject.threads_init()
        EegpyWin.__init__(self, x_in, x_out)
        self.setupOptions()
        self.show_all()
    
    def delete_event(self, widget, event, data=None):
        if self.compCh != None:
            self.compCh.destroy()
        EegpyWin.delete_event(self, widget, event, data)
        return False
 
    def cb_Buttons(self, widget):
        if widget == self.varChIn:
            #try:
                self.log("Working on array with shape"+str(self.varChIn.get_var().shape))
                self.showFilter()
            #except AttributeError:
            #    self.log("Given Variable is not an array!!!")
        elif widget == self.varChOut:
            try:
                if self.y!=None:
                    self.log("Saving results...")
                    self.varChOut.set_var(self.y)
                    self.compCh.hide()
                    self.cb_quit(widget)
            except AttributeError:
                self.log("Given Variable is not an array!!!")
        pass    

    def cb_Options(self, widget):
        if widget == self.bOpt["showCA"]:
            self.x = self.varChIn.get_var()
            if self.fTypeW["PCA"].get_active():
                self.x2 = pcafilt.unmix(self.x)
            elif self.fTypeW["ICA"].get_active():
                self.x2 = icafilt.unmix(self.x)
            if self.compCh==None:
                self.compCh = ComponentChooser(self.x2)
            else:
                self.compCh.show_all()
        
    def setupOptions(self):
        self.bTable = gtk.Table(2, 2, True)
        #SpinBoxes for Frequencies and order
        self.bOpt = {}
        self.bOpt["showCA"] = gtk.Button("Analyze")
        self.bTable.attach(self.bOpt["showCA"],0,2,0,1)
        self.bOpt["showCA"].connect("clicked",self.cb_Options)
     
        #Radio-Buttons
        self.fTypeW = {}
        self.fTypeW["PCA"] = gtk.RadioButton(None,"PCA")
        self.bTable.attach(self.fTypeW["PCA"],0,1,1,2)
        self.fTypeW["ICA"] = gtk.RadioButton(self.fTypeW["PCA"],"ICA")
        self.bTable.attach(self.fTypeW["ICA"],1,2,1,2)
        
        self.upper_hbox.pack_start(self.bTable, expand=False, padding=1)
        
        #self.bTable2 = gtk.Table(1, 1, True)        
        #self.btAnal = gtk.Button("Show/Hide Analysis")
        #self.btAnal.connect("clicked", self.cb_showAnalysis)
        #self.bTable2.attach(self.btAnal,0,1,0,1)
        #self.hbox.pack_start(self.bTable2, False, padding=50)
    
    def setupSubplots(self):
        self.f.clear()
        self.a = self.f.add_subplot(MultichannelSubplot(self.f,1,1,1))  
    
    def showFilter(self):
        if self.compCh != None:
            if self.fTypeW["PCA"].get_active():
                self.y = pcafilt.mix(self.x2,self.compCh.choose_mask)
            elif self.fTypeW["ICA"].get_active():
                self.y = icafilt.mix(self.x2,self.compCh.choose_mask)
            self.a.cla()
            self.a.plotChannels(self.x,"k")
            self.a.plotChannels(self.y,"r")
            #for i in range(17):
            #    self.a.plot(self.x[:,i])
            self.canvas.draw()
        

def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    Wert = 9
    ar = numpy.zeros((3,5,7),"d")
    try:
        import eegpy.formats.edf
        r = eegpy.formats.edf.EdfReader("/home/thorsten/Documents/Downloads/Newtest17-256.bdf")
        ar2 = r.getData(0,2000)
        ar3 = r.getData(0,1000)
    except Exception:
        pass
    pw = PCAFiltWin(ar2)
    main()