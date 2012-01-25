# -*- coding: utf-8 -*-
import pygtk
pygtk.require('2.0')
import gobject
import gtk
import pickle
from Eeg.Formats import f32
from Eeg.Formats import bvFormat
import eegpy
from eegpy.events import EventTable
from eegpy.filter.freqfilt import filtfilt
#import pylab
import os.path
import sys, time
import numpy as n
from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.axes import Subplot
from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
from matplotlib.figure import Figure, SubplotParams
from matplotlib.axis import Axis
from scipy.interpolate import interp1d
from scipy.signal import butter


class Marker:
    """For marking"""
    def __init__(self, plot=None):
        self.plot = plot
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_default_size(200,300)
        self.window.set_title("Making marks")

        # Set a handler for delete_event that immediately
        # exits GTK.
        self.window.connect("delete_event", self.delete_event)
        self.window.connect("destroy", self.delete_event)
        
        self.mainBox = gtk.VBox()
        self.window.add(self.mainBox)
        
        self.tvScrolledWin = gtk.ScrolledWindow()
        self.tvScrolledWin.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
        self.tree = gtk.TreeStore(gobject.TYPE_INT)
        #self.treeS = gtk.TreeModelSort(self.tree)
        self.treeV = gtk.TreeView(self.tree)
        self.treeV.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        renderer = gtk.CellRendererText()
        self.col1 = gtk.TreeViewColumn("Marker at ...", renderer,text=0)
        self.treeV.append_column(self.col1)
        self.treeV.show()
        self.tvScrolledWin.add(self.treeV)
        self.tvScrolledWin.show_all()
        #self.hbox.pack_start(self.tvScrolledWin)
        self.mainBox.pack_start(self.tvScrolledWin)
        
        self.cbx_pauseMarking = gtk.CheckButton("pause marking")
        self.mainBox.pack_start(self.cbx_pauseMarking,expand=False)
        
        self.btRemove = gtk.Button("Remove marked")
        #self.btRemove.set_size_request(-1, 20)
        self.btRemove.connect("clicked", self.cb_remove)
        self.mainBox.pack_start(self.btRemove,expand=False)
        self.btSave = gtk.Button("Save to EventTable")
        #self.btRemove.set_size_request(-1, 20)
        self.btSave.connect("clicked", self.cb_save)
        self.mainBox.pack_start(self.btSave,expand=False)
        
        self.window.show_all()
        
        
    
    def delete_event(self, widget, event, data=None):
        self.window.hide()
        return True
    
    def cb_remove(self, widget):
        def remove(model, path, iter):
            model.remove(iter)
            return False     # keep the foreach going
        
        pathlist = self.treeV.get_selection().get_selected_rows()[1]
        iterlist = [self.tree.get_iter(row) for row in pathlist]
        for row in iterlist:
            self.tree.remove(row)
        self.plot.plot_data()
            #print row
    
    def cb_save(self, widget):
        marklist = []
        def append(model, path, iter, user_data):
            marklist.append(self.tree.get(iter,0)[0])
        dialog_label = gtk.Dialog("My dialog", None, gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT, (gtk.STOCK_CANCEL, gtk.RESPONSE_REJECT, gtk.STOCK_OK, gtk.RESPONSE_OK))    
        entry1 = gtk.Entry()
        entry1.set_text("Marks")
        dialog_label.vbox.pack_start(entry1)
        entry1.show()
        response = dialog_label.run()
        print response
        if response == gtk.RESPONSE_OK:
            trig_name = entry1.get_text()
            print trig_name
        else:
            print "Saving aborted by user."
            dialog_label.destroy()
            return False    
        dialog_label.destroy()
        
        dialog = gtk.FileChooserDialog("Save list as pickle-file ...", None, gtk.FILE_CHOOSER_ACTION_SAVE, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("EventTables")
        filter.add_pattern("*.evt")
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
            et = EventTable({trig_name:marklist})
            et.save(dialog.get_filename())
            #fh = open(dialog.get_filename(),"w")
            #pickle.dump(marklist,fh,-1)
            print dialog.get_filename(), 'selected'
            #fh.close()
        elif response == gtk.RESPONSE_CANCEL:
            print 'Closed, no files selected'
        dialog.destroy()
    
    def getMarks(self,rangeToTake=None):
        marklist = []
        def append(model, path, iter, user_data):
            xVal = self.tree.get(iter,0)[0]
            if rangeToTake==None:
                marklist.append(xVal)
            else:
                assert len(rangeToTake)>1, "rangeToTake must be of length 2"
                if xVal>rangeToTake[0] and xVal<rangeToTake[1]:
                    marklist.append(xVal)
        self.tree.foreach(append, "")
        marklist.sort()
        return marklist
    
    def add(self, tp):
        if not self.cbx_pauseMarking.get_active():
            iter = self.tree.append(None)
            self.tree.set(iter, 0,tp)  
            self.tree.set_sort_column_id(0,gtk.SORT_ASCENDING) 
             

class MarkerWithAverage(Marker):
    """Erweiterung der Marker-Klasse. Ermöglicht das erstellen von Mittelungen um die Markierungen herum."""
    f = None
    a = None
    toolbar = None
    canvas = None
    data = None
    combo = None
    pbar = None
    
    def __init__(self, plot=None):
        Marker.__init__(self,plot)
        self.window.resize(200,500)
        #SPinbutton für Breite des Fensters
        self.mainBox.pack_start(gtk.VSeparator(),expand=False)
        self.hboxWidth=gtk.HBox()
        self.mainBox.pack_start(self.hboxWidth,expand=False)
        self.hboxWidth.pack_start(gtk.Label("Breite"))
        self.sbWidth = gtk.SpinButton(gtk.Adjustment(600,0,10000,10,100))
        self.hboxWidth.pack_end(self.sbWidth,expand=False)
        self.btCalc = gtk.Button("Calculate Avg")
        self.btCalc.connect("clicked", self.cb_calcAvg)
        self.mainBox.pack_start(self.btCalc, expand=False)
        
        #Setup Bandpass-Filter
        sr = self.plot.reader.samplRate
        [b,a]=butter(3,[3./(sr/2),13./(sr/2)], btype="band")
        self.bp_filt = lambda x: filtfilt(b,a,x)
        
        self.window.show_all()
    
    def cb_calcAvg(self, widget):
        marks = self.getMarks()
        assert len(marks) > 0, "No marks set yet!"
        width = self.sbWidth.get_value()
        self.data = n.zeros((width,self.plot.reader.numChannels),"d")
        
        #print "Try to aquire lock..."
        self.plot.f32Lock.acquire()
        for mark in marks:
            #print mark
            self.data += self.plot.reader.getData(int(mark-width/2),int(width))
        self.plot.f32Lock.release()
        self.data/=len(marks)
        self.data = self.bp_filt(self.data)
        #print self.data
        self.plotAvg()

    def cb_combo(self, widget):
        cnum = self.combo.get_active()
        self.plotAvg(cnum)
        
    def cb_findMatches(self, widget):
        start=int(self.plot.sbStartpoint.get_value())
        length=int(self.plot.sbNDataP.get_value())
        stride=int(self.plot.sbStride.get_value())
        data2 = self.plot.reader.getData(start,length*stride,1)
        corr = self.correlate2d(self.data, data2)
        self.canvas.hide()
        self.a.clear()
        #xs = n.arange(-int(self.sbWidth.get_value()/2), -int(self.sbWidth.get_value()/2)+corr.shape[0])
        self.a.plot(corr)
        #self.a.set_xticks([])
        self.a.set_yticks([])
        for lab in self.a.get_xticklabels():
            lab.set_fontsize(6)
        self.canvas.show()
        self.findArtifacts()

                    
    def cb_subtract(self, widget):
        assert self.data != None, "Data must not be None!"    
        start=int(self.plot.sbStartpoint.get_value())
        stride=int(self.plot.sbStride.get_value()) #Stride ermitteln...
        
        stData = self.data.copy()
        #for i in range(stData.shape[1]):
        #    stData[:,i] *= n.hanning(stData.shape[0])
        #print stData
        chList = self.plot.chList
        arts = [x/stride for x in self.findArtifacts()]
        #arts = self.findArtifacts()
        print "Arts",arts
        for i in arts:
            #startSt = (i - stData.shape[0]/2)/stride
            #endSt = (i + stData.shape[0]/2)/stride
            #print "start/endSt:", startSt, endSt
            #startCorrect = 0
            #endCorrect = 0
            #if startSt<0:
            #    startCorrect = -startSt
            #    startSt=0
            #if endSt>self.plot.data.shape[0]:
            #    endCorrect = endSt-self.plot.data.shape[0]
            #    endSt = self.plot.data.shape[0]
            #print "start/endSt:", startSt, endSt
            for j in range(self.plot.data.shape[1]):
                #xs = n.arange((i - stData.shape[0]/2),(i + stData.shape[0]/2))
                #artInterp = interp1d(xs,stData[:,j])
                #print n.arange((i - stData.shape[0]/2),(i + stData.shape[0]/2))
                startSt = (i - stData.shape[0]/2)/stride
                endSt = (i + stData.shape[0]/2)/stride
                #print "start/endSt:", startSt, endSt
                if startSt<0:
                    startSt=0
                if endSt>self.plot.data.shape[0]:
                    endSt = self.plot.data.shape[0]
                #print "start/endSt:", startSt, endSt
                for k in range(startSt,endSt):
                    #print xs[0], xs[-1], k
                    #print self.plot.data[k,j],artInterp(k), self.plot.data[k,j] - artInterp(k)
                    self.plot.data[k,j] -= stData[k-startSt,self.plot.chList[j]]
                #print "---", self.plot.data[startSt:endSt,j].shape, stData[startCorrect:-endCorrect:stride,chList[j]].shape
                #self.plot.data[startSt:endSt,j] -= stData[startCorrect*stride:-endCorrect*stride:stride,chList[j]]
        self.plot.plotData()
     
    def cb_findAll(self, widget):
        """Searches for all BCG-artifacts in file and adds them to the list"""
        #create ProgressBar
        if self.pbar==None:
            self.pbar = gtk.ProgressBar()
            self.pbar.set_text("Processing...")
            self.mainBox.pack_start(self.pbar,expand=False)
        self.pbar.show()
        
        positions = []
        #loop through file
        self.plot.f32Lock.acquire()
        #Tree leeren
        self.tree.clear()
        #x = self.plot.reader[:].copy() #Uses F32-class!!!
        
        ##################
        if self.cbx_useLogic.get_active():
            delWidth=500
        else:
            delWidth=300
        corr = n.array([],"d")
        bsl_ch = n.zeros((self.plot.reader.numChannels),n.bool)
        for ch in self.plot.chList:
            bsl_ch[ch] = True
        startpoints = range(0,self.plot.reader.numDatapoints,10000)
        for i,sp in enumerate(startpoints):
            data2 = self.plot.reader.getData(sp,10000)
            corr = n.r_[corr,self.correlate2d(self.data[:,bsl_ch], data2[:,bsl_ch])]
            self.pbar.set_fraction(float(i)/len(startpoints)*0.7) #At the end, 70% completed
            while gtk.events_pending():
                gtk.main_iteration_do(False)
        print "corr", corr.shape
        try: #Finde Artefakte
            corrC = corr#.copy()
            #self.a.clear()
            #self.a.plot(corrC[::10])
            alocs = []
            allMax = 0
            i=0
            #First round
            while True:#for i in range():
                aloc = corrC.argmax()
                if i == 0:
                    allMax = corrC[aloc]
                #print "corrC[aloc]/allMax", corrC[aloc]/allMax
                if corrC[aloc]/allMax>0.3:
                    alocs.append(aloc)
                else:
                    break
                #Bereich um lokales Maximum auslöschen
                width = delWidth/2
                startDel = aloc-width
                endDel = aloc+width
                if startDel<0:
                    startDel=0
                if endDel>corrC.shape[0]:
                    endDel=corrC.shape[0]
                corrC[startDel:endDel] = n.zeros((endDel-startDel),"d")
                i+=1
            alocs.sort()
            #Second round: correction
            dPos_dT = n.diff(n.array(alocs))
            print "dPos shape 1", dPos_dT.shape
            dPos_dT = dPos_dT[dPos_dT<1500]
            print "dPos shape 2", dPos_dT.shape
            #dPos_dT = dPos_dT[dPos_dT>500]
            #print "dPos shape 3", dPos_dT.shape
            pulse = int(dPos_dT.mean())
            print "Pulse:", pulse
            #print minx, maxx
            #print "Positionen der Artefakte:", alocs
            print len(alocs), alocs
            positions_filt = []
            for i,pos in enumerate(alocs):
                self.pbar.set_fraction(0.7+0.3*(float(i)/len(alocs)))
                while gtk.events_pending():
                    gtk.main_iteration_do(False)
                if i==0:
                    positions_filt.append(pos)
                elif pos-positions_filt[-1]>500:
                    positions_filt.append(pos)
                    try:
                        while (positions_filt[-1]-positions_filt[-2])>1.6*pulse:
                            #Setze
                            step=1
                            aloc = positions_filt[-2]+pulse
                            step=2
                            aloc = corrC[aloc-100:aloc+100].argmax()+aloc-100
                            step=3
                            positions_filt.insert(-1,aloc)
                    except IndexError,e:
                        print "IndexError:", e, "Step", step
        except Exception, e:
            print "Fehler beim finden der Artefakte:", e
        positions_filt.sort()
        
        for i,pos in enumerate(positions_filt):
            self.add(pos)
        #remove ProgressBar
        self.pbar.destroy()
        self.pbar = None   
        print "positions_filt: ", len(positions_filt), positions_filt  
        self.plot.f32Lock.release()    
    
    def cb_saveRemoved(self, widget):
        def filterToFile(fn, winsize=5500):
            #f32Lock holen
            self.plot.f32Lock.acquire()
            #create ProgressBar
            if self.pbar==None:
                self.pbar = gtk.ProgressBar()
                self.pbar.set_text("Processing...")
                self.mainBox.pack_start(self.pbar,expand=False)
            self.pbar.show()
            #EEG-writer
            out = f32.F32WriterAdvanced(fn,self.plot.reader.channel_names)
            stData = self.data.copy()
            startpoints = range(0,self.plot.reader.numDatapoints,winsize)
            for i,sp in enumerate(startpoints):
                length = winsize
                if sp+length>self.plot.reader.numDatapoints:
                    length=self.plot.reader.numDatapoints-sp-1
                tmpdata = self.plot.reader.getData(sp,length)
                arts = self.getMarks(rangeToTake=[sp,sp+length])
                for a in [x-sp for x in arts]:
                    for ch in range(tmpdata.shape[1]):
                        startSt = (a - stData.shape[0]/2)
                        endSt = (a + stData.shape[0]/2)
                        #print "start/endSt:", startSt, endSt
                        if startSt<0:
                            startSt=0
                        if endSt>tmpdata.shape[0]:
                            endSt = tmpdata.shape[0]
                        #print "start/endSt:", startSt, endSt
                        for k in range(startSt,endSt):
                            #print xs[0], xs[-1], k
                            #print self.plot.data[k,j],artInterp(k), self.plot.data[k,j] - artInterp(k)
                            #print "k,ch,startSt",k,ch, startSt
                            tmpdata[k,ch] -= stData[k-startSt,ch]
                out.appendData(tmpdata)
                self.pbar.set_fraction(float(i)/len(startpoints))
                self.pbar.set_text("Processed %i of %i"%(i,len(startpoints)))
                while gtk.events_pending():
                    gtk.main_iteration_do(False)
            out.writeNumDatapoints()
            del out
            #f32Lock freigeben
            self.plot.f32Lock.release()
            #remove ProgressBar
            self.pbar.destroy()
            self.pbar = None 
                        
                    
        dialog = gtk.FileChooserDialog("Save filtered EEG-File...", None, gtk.FILE_CHOOSER_ACTION_SAVE, (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        dialog.set_default_response(gtk.RESPONSE_OK)
        
        filter = gtk.FileFilter()
        filter.set_name("F32-EEG-File")
        filter.add_pattern("*.f32")
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
            filterToFile(fn)
            #fh.close()
        else:# response == gtk.RESPONSE_CANCEL:
            dialog.destroy()
            print 'Closed, no files selected'
        pass
            
        
    def plotAvg(self,cnum=10):
        if self.a == None:
            self.addCanvas()
        self.canvas.hide()
        self.a.clear()
        self.a.plot(self.data[:,cnum])
        self.a.set_xticks([])
        #self.a.set_yticks([])
        #self.a.set_ylim(self.data.min(),self.data.max())
        #print self.a.get_xticklabels()[0]
        #
        #print [t.get_text() for t in self.a.get_yticklabels()]
        for tick in self.a.yaxis.get_major_ticks():
            tick.set_pad(0)
        for lab in self.a.get_yticklabels():
            lab.set_fontsize(6)
             
        #self.a.set_yticklabels([t.get_text() for t in self.a.get_yticklabels()], fontsize=8) 
        self.canvas.show()
        #In entry1: schreibe min bzw. max / std()
        msd = abs(self.data[:,cnum]).max()/self.data[:,cnum].std()
        self.entry1.set_text("%.2f"%msd)
        
    def addCanvas(self):
        #Canvas
        self.f = Figure(figsize=(5,4), dpi=100, subplotpars=SubplotParams(left=0.06, top=0.95, right=0.97, bottom=0.1,hspace=0))
        self.a = self.f.add_subplot(111)
        #self.a.yaxis.set_pad(0)
        #p.rc("ytick.major", pad=0 )
        self.canvas = FigureCanvas(self.f)
        #self.canvas.show()
        self.mainBox.pack_start(self.canvas)
        #self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        #self.toolbar = NavigationToolbar( self.canvas, self.window )
        self.combo = gtk.combo_box_new_text()
        #self.combo.set_value_in_list(True,False)
        #self.combo.set_popdown_strings(self.plot.reader.channel_names)
        for i in self.plot.reader.channel_names:
            self.combo.append_text(i)
        self.combo.set_active(0)
        self.combo.connect("changed", self.cb_combo)
        self.mainBox.pack_start(self.combo,expand=False)
        self.entry1 = gtk.Entry(100)
        self.entry1.set_editable(False)
        self.entry1.set_text("---")
        self.mainBox.pack_start(self.entry1,expand=False)
        self.btCorr = gtk.Button("Correlate")
        self.btCorr.connect("clicked",self.cb_findMatches)
        self.mainBox.pack_start(self.btCorr,expand=False)
        self.btSubtr = gtk.Button("Subtract artifacts")
        self.btSubtr.connect("clicked",self.cb_subtract)
        self.mainBox.pack_start(self.btSubtr,expand=False)
        self.btFindAll = gtk.Button("Find all matches")
        self.btFindAll.connect("clicked",self.cb_findAll)
        self.mainBox.pack_start(self.btFindAll,expand=False)
        self.cbx_useLogic = gtk.CheckButton("use additional logic")
        self.cbx_useLogic.set_active(True)
        self.mainBox.pack_start(self.cbx_useLogic,expand=False)
        self.btSaveRemoved = gtk.Button("Save file without BCG")
        self.btSaveRemoved.connect("clicked",self.cb_saveRemoved)
        self.mainBox.pack_start(self.btSaveRemoved,expand=False)
        self.mainBox.show_all()
    
    def correlate2d(self, x,y):
        assert len(x.shape)==2, "Array x must be 2d"
        assert len(y.shape)==2, "Array y must be 2d"
        #print x.shape, y.shape
        assert x.shape[1] == y.shape[1], "Number of channels must be equal in both arrays."
        
        #Bandpass on array y (x is self.data, already filtered)
        y = self.bp_filt(y)
        #weights = n.array([0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        #weights = n.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
        #weights = n.zeros((x.shape[1]),"d")
        #for i in self.plot.chList:
        #    weights[i] = 1
        #print weights
        #weights = weights / (weights.mean()*weights.shape[0])
        corr = n.correlate(x[:,0],y[:,0],"same")
        for i in range(1,x.shape[1]):
            corr += n.correlate(x[:,i],y[:,i],"same")
        corr/=x.shape[1]
        #corr /= x.shape[1]
        #print corr.shape
        return corr
    
        #for
    
    
#        def sqrError(x,y):
#            err = 0
#            for i in range(x.shape[0]):
#                err+=(x[i]-y[i])**2
#            return err
#                
#        corr = x[:,0].copy()
#        for i in range(x.shape[0]):
#            corr[i] = sqrError(x[i,:],y[100,:])
#        return corr

    def findArtifactsIn(self,data2,delWidth=500):
        corr = self.correlate2d(self.data, data2)
        try: #Finde Artefakte
            corrC = corr#.copy()
            #self.a.clear()
            #self.a.plot(corrC[::10])
            alocs = []
            allMax = 0
            i=0
            while True:#for i in range():
                aloc = corrC.argmax()
                if i == 0:
                    allMax = corrC[aloc]
                #print "corrC[aloc]/allMax", corrC[aloc]/allMax
                if corrC[aloc]/allMax>0.2:
                    alocs.append(aloc)
                else:
                    break
                #Bereich um lokales Maximum auslöschen
                width = delWidth/2
                startDel = aloc-width
                endDel = aloc+width
                if startDel<0:
                    startDel=0
                if endDel>corrC.shape[0]:
                    endDel=corrC.shape[0]
                corrC[startDel:endDel] = n.zeros((endDel-startDel),"d")
                i+=1
            #print minx, maxx
            #print "Positionen der Artefakte:", alocs
            return alocs
        except Exception, e:
            print "Fehler beim finden der Artefakte:", e
            
    def findArtifacts(self):
        start=int(self.plot.sbStartpoint.get_value())
        length=int(self.plot.sbNDataP.get_value())
        stride=int(self.plot.sbStride.get_value())
        data2 = self.plot.reader.getData(start,length*stride,1)
        
        self.plot.canvas.hide()
        alocs = self.findArtifactsIn(data2)
        minx = min(self.plot.ts)
        maxx = max(self.plot.ts)
        for aloc in alocs:
            self.plot.a.axvline((start+aloc)/self.plot.tsFactor, lw=2, alpha=0.5, color="#FF0000")
        self.plot.subplAxes[0].set_xlim((minx, maxx))
        self.plot.canvas.show()
        return alocs
        

        
        