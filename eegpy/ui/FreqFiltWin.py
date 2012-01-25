#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################
# Module-Import #
#################

#eegpy-modules
try:
    import eegpy
    from eegpy.misc import FATALERROR
    from eegpy.filter.freqfilt import filtfilt, butter
    from eegpy.ui.widgets.iowidgets import VarChooser
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

class FreqFiltWin(EegpyWin):
    programName = "eegpy: Frequency-Filtering"
    
    # Konstruktor
    def __init__(self, x_in=123456, x_out=123456):
        EegpyWin.__init__(self, x_in, x_out)
        self.setupOptions()
        self.show_all()
        #self.setupGUI()
         
    def cb_Buttons(self, widget):
        if widget == self.varChIn:
            try:
                self.log("Working on array with shape"+str(self.varChIn.get_var().shape))
                self.showFilter()
            except AttributeError:
                self.log("Given Variable is not an array!!!")
        elif widget == self.varChOut:
            try:
                if self.y!=None:
                    self.log("Saving results...")
                    self.varChOut.set_var(self.y)
                    self.cb_quit(widget)
            except AttributeError:
                self.log("Given Variable is not an array!!!")
        pass
    
    def cb_fOpt(self, widget):
        #print "cb_fType called"
        myDict= {}
        k=""
        for i in self.fTypeW.keys():
            if self.fTypeW[i] == widget:
                myDict=self.fTypeW
                k=i
        for i in ["order","low","high","Fs"]:
            if self.fOpt[i] == widget:
                myDict=self.fOpt
                k=i
        if myDict == self.fTypeW:
            #Do things like deactivating SpinButtons
            pass
        elif myDict == self.fOpt:
            #Change Spinboxes
            if k == "Fs":
                self.fOpt["high"].set_adjustment(gtk.Adjustment(self.fOpt["high"].get_value(),0,self.fOpt[k].get_value()/2,1,10))
                #print k
        if k == "":
            print "None of the RadioButtons was chosen."
            return False
        
    def setupOptions(self):
        self.bTable = gtk.Table(2, 4, True)
        #SpinBoxes for Frequencies and order
        self.fOpt = {}
        self.fOpt["order"] = gtk.SpinButton(gtk.Adjustment(3,1,6,1,1))
        self.bTable.attach(self.fOpt["order"],0,1,0,1)
        self.fOpt["low"] = gtk.SpinButton(gtk.Adjustment(0.0,0.0,5000.0,1,10))
        self.bTable.attach(self.fOpt["low"],1,2,0,1)
        self.fOpt["high"] = gtk.SpinButton(gtk.Adjustment(500.0,1.0,5000.0,1,10))
        self.bTable.attach(self.fOpt["high"],2,3,0,1)
        self.fOpt["Fs"] = gtk.SpinButton(gtk.Adjustment(1000.0,1.0,10000.0,1,100))
        self.bTable.attach(self.fOpt["Fs"],3,4,0,1)
        for k in self.fOpt.keys():
            self.fOpt[k].set_numeric(True)
            self.fOpt[k].connect("value_changed",self.cb_fOpt)
        #Radio-Buttons
        self.fTypeW = {}
        self.fTypeW["low"] = gtk.RadioButton(None,"Lp")
        self.bTable.attach(self.fTypeW["low"],0,1,1,2)
        self.fTypeW["high"] = gtk.RadioButton(self.fTypeW["low"],"Hp")
        self.bTable.attach(self.fTypeW["high"],1,2,1,2)
        self.fTypeW["bs"] = gtk.RadioButton(self.fTypeW["low"],"Bs")
        self.bTable.attach(self.fTypeW["bs"],2,3,1,2)
        self.fTypeW["bp"] = gtk.RadioButton(self.fTypeW["low"],"Bp")
        self.bTable.attach(self.fTypeW["bp"],3,4,1,2)
        for k in self.fTypeW.keys():
            self.fTypeW[k].connect("clicked", self.cb_fOpt)
        
        self.upper_hbox.pack_start(self.bTable, expand=False, padding=1)
        
        #self.bTable2 = gtk.Table(1, 1, True)        
        #self.btAnal = gtk.Button("Show/Hide Analysis")
        #self.btAnal.connect("clicked", self.cb_showAnalysis)
        #self.bTable2.attach(self.btAnal,0,1,0,1)
        #self.hbox.pack_start(self.bTable2, False, padding=50)
        
    def setupSubplots(self):
        self.f.clear()
        
        self.a = self.f.add_subplot(211)
        self.a2 = self.f.add_subplot(212)      
    
    def showFilter(self):
        k="low"
        for i in self.fTypeW.keys():
            if self.fTypeW[i].get_active():
                k=i
        self.x = self.varChIn.get_var()
        Wn = []
        if k in ["low","bs","bp"]:
            Wn.append(float(self.fOpt["low"].get_value())/self.fOpt["Fs"].get_value()*2)
        if k in ["high","bs","bp"]:
            Wn.append(float(self.fOpt["high"].get_value())/self.fOpt["Fs"].get_value()*2)
        #print "butter", int(self.fOpt["order"].get_value()), Wn, k
        b,a = butter(int(self.fOpt["order"].get_value()),Wn,k)
        #print "filtfilt"
        self.y = filtfilt(b,a,self.x)
        #print "plot"
        self.canvas.hide()
        self.a.cla()
        self.a.plot(self.x[:,0])
        self.a.plot(self.y[:,0])
        self.a2.cla()
        self.a2.psd(self.x[:,0],Fs=self.fOpt["Fs"].get_value())
        self.a2.psd(self.y[:,0],Fs=self.fOpt["Fs"].get_value())
        self.a2.autoscale_view(tight=True)
        self.canvas.show()
        

def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    ffw = FreqFiltWin()
    Wert = 9
    ar = numpy.zeros((3,5,7),"d")
    try:
        import eegpy.formats.edf
        r = eegpy.formats.edf.EdfReader("/home/thorsten/Documents/Downloads/Newtest17-256.bdf")
        ar2 = r.getData(0,2000)
        ar3 = r.getData(0,1000)
    except Exception:
        pass
    main()