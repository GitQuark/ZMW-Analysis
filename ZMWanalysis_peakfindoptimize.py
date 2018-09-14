# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 11:55:59 2014

@author: Robert Henley
"""
import sys
import os
import numpy as np
from zmwanalysiswidget import *
from filterdialog import *
from threshdialogui import *
import pyqtgraph as pg
import h5py
from scipy.ndimage import filters
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import imreg_dft as ird
from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from peakdetect import peakdet
import peakFindingZMW as pfz
from scipy.signal import savgol_filter



class GUIForm(QtGui.QMainWindow):
    
    def __init__(self, master=None):
        QtGui.QMainWindow.__init__(self,master)
        self.ui = Ui_ZmwAnalysisWidget()
        self.ui.setupUi(self)
        self.filterdialog = filterpopup()
        self.analysisdialog = analysispopup()
        
#        self.plotlist = [[],[],[],[],[]]
        self.seqplotlist = [[]] 
        self.firingplotlist = []


        self.ui.actionLoad.triggered.connect(self.getfile)
        self.ui.actionZ_Project.triggered.connect(self.zproject)
        self.ui.actionView_Stack.triggered.connect(self.viewstack)
        self.ui.actionClear_Settings.triggered.connect(self.clearsettings)
        self.ui.actionFilter.triggered.connect(self.showfilterdialog)
        self.ui.actionCrop.triggered.connect(self.crop)
        self.ui.actionBase_Call.triggered.connect(self.showanalysisdialog)
#        self.ui.actionCheck_Controls.triggered.connect(self.check_controls)
        self.filterdialog.acceptsignal.connect(self.filterstack)
        self.analysisdialog.acceptsignal.connect(self.analyze)

        
        self.p1 = self.ui.imagePlot
        self.p1.roi.scaleSnap = self.p1.roi.translateSnap = True
        self.p1.roi.setSize([6,17])
        colors = [
            (0, 0, 0),
            (45, 5, 61),
            (84, 42, 55),
            (150, 87, 60),
            (208, 171, 141),
            (255, 255, 255)
        ]
        cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=colors)
        self.p1.setColorMap(cmap)
        self.direc = []
        self.roip1 = self.roip2 = []
        self.xfilt = self.yfilt = self.zfilt = self.filtertype = []
        self.p1.getRoiPlot().enableAutoRange(False)



    def getfile(self):

        if self.direc==[]:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',os.getcwd(),("*.h5"))[0]
            self.direc=os.path.dirname(str(self.datafilename))
            print(self.datafilename)
        else:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self,'Open file',self.direc,("*.h5"))[0]
            self.direc=os.path.dirname(str(self.datafilename))
            print(self.datafilename)


        self.Load()

    def Load(self):

        if self.seqplotlist != [[]]  :
            for x in self.seqplotlist:
                self.p1.getRoiPlot().removeItem(x)
                
        if self.firingplotlist != []  :
            for x in self.firingplotlist:
                self.p1.getRoiPlot().removeItem(x)

        self.seqplotlist = [[]]
        self.firingplotlist = []
       
        self.stack = h5py.File(self.datafilename)
        self.stack = np.array(self.stack["images"]).astype(float)
        self.stack = self.backgroundsubtract(self.stack)
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
        self.originalstack = self.stack
        if self.roip1 == []:
            pass
        else:
            self.stack = self.stack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]           
            
#        self.stack -= float(self.stack.mean())
#        self.stack /= float(self.stack.std())
#        baseline = np.median(self.stack[self.stack < np.median(self.stack)])
#        self.background = self.stack[np.max(self.stack,(1,2)) < baseline]
#        self.stack = self.subtract_background(self.stack)
#        self.stack -= np.mean(self.stack)
#        self.stack = self.stack/self.stack.std()
#        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        self.p1.roi.setSize([6,17])
        
#        if self.xfilt != []:
#            self.filter_stack(self.xfilt,self.yfilt,self.zfilt,self.filtertype)
#            self.analyze(self.blueThresh,self.redThresh,
#                               self.bluePeakThresh,self.redPeakThresh)

    def backgroundsubtract(self,stack):
        avg = self.stack.mean((1,2))
        self.background = stack[avg < avg.mean()-.5*avg.std()]
        self.background = np.median(self.background,0)
        stack -= self.background
        
        return stack

    def zproject(self):
        self.p1.setImage(np.transpose(256*(1+self.zpro),(1,0)))

    def viewstack(self):
        self.p1.setImage(np.transpose((1+self.stack), (0,2,1)))
        self.p1.play(rate=60)
        
    def showfilterdialog(self):
        self.filterdialog.show()
        
    def filterstack(self,xfilt,yfilt,zfilt,filtertype):
        if self.roip1 == []:
            pass
        else:
            self.stack = self.originalstack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]  
        self.xfilt = xfilt
        self.yfilt = yfilt
        self.zfilt = zfilt        
        self.filtertype = filtertype
        if filtertype == 0:
            self.stack = filters.median_filter(
                self.stack, size = (xfilt,yfilt,zfilt))
        if filtertype ==1:
            self.stack = filters.gaussian_filter(
                self.stack, sigma = (xfilt,yfilt,zfilt))
        self.zpro = np.std(self.stack,0)
        self.viewstack()
        
    def clearsettings(self):
        self.roip1 = []
        self.xfilt = []

    def closedialog(self):
        self.filterdialog.destroy()
        
    def crop(self):
        self.roip1,self.roip2 =  np.array(self.p1.roi.pos()).astype(int), np.array(self.p1.roi.pos()).astype(int)+ np.array(self.p1.roi.size()).astype(int)
        self.stack = self.stack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]
        self.p1.roi.setPos(0,0)
        self.p1.roi.setSize((int(self.p1.roi.size()[0]),int(self.p1.roi.size()[1])))
#        self.stack = self.subtract_background(self.stack)
#        self.stack -= np.mean(self.stack)
#        self.stack = self.stack/self.stack.std()
#        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        self.zpro = np.std(self.stack,0)
        self.viewstack()
        
    def showanalysisdialog(self):
        self.analysisdialog.show()
        
    def definecontrols(self):
        dntpnames = ["dCTP","dATP","dGTP","dTTP"]
        dntps = [[],[],[],[]]
        zpro = [[],[],[],[]]
        
        dntpdirec = r'C:/Users/Meni/Robert/'
        fn = 'dntps.h5'
        hfile = h5py.File(dntpdirec + fn)
        for i,x in enumerate(dntpnames):
            fn = 'dntps.h5'
            zpro[i] = np.array(hfile[x]).astype(float)
            zpro[i] -=zpro[i].mean()
            zpro[i] /=zpro[i].std()
        
        composite = zpro[0] + zpro[1] + zpro[2] + zpro[3]
        ly, lx = composite.shape
        
        p0 = [composite.mean(1).max(), (ly/2)-4, (ly/2)+4, 1.]
        coeff1, var_matrix1 = curve_fit(dubgauss,np.linspace(0,ly-1, ly), composite.mean(1), p0=p0)
        p0 = [self.zpro.mean(1).max(), (ly/2)-4, (ly/2)+4, 1.]
        coeff2, var_matrix2 = curve_fit(dubgauss,np.linspace(0,ly-1, ly), self.zpro.mean(1), p0=p0)
        shifty = np.mean((coeff2[1], coeff2[2]))-np.mean((coeff1[1], coeff1[2]))
        
        p0 = [composite.mean(1).max(), lx/2, 1.]
        coeff1, var_matrix1 = curve_fit(gauss,np.linspace(0,lx-1, lx), composite.mean(0), p0=p0)
        p0 = [self.zpro.mean(1).max(), lx/2, 1.]
        coeff2, var_matrix2 = curve_fit(gauss,np.linspace(0,lx-1, lx), self.zpro.mean(0), p0=p0)
        shiftx = coeff2[1]-coeff1[1]
        
        self.czpro = [[],[],[],[]]
        for i,x in enumerate(dntps):
            self.czpro[i] = zpro[i]
            self.czpro[i] = ird.transform_img(self.czpro[i], tvec=[shifty,shiftx])
        
        
    def analyze(self,blueThresh,redThresh,bluesubThresh, redsubThresh):
        self.definecontrols()
        
        seqdf = pd.DataFrame({'base':[],'times':[]})
        if self.seqplotlist != [[]]  :
            for x in self.seqplotlist:
                self.p1.getRoiPlot().removeItem(x)
            self.seqplotlist = []
                
        if self.firingplotlist != []  :
            for x in self.firingplotlist:
                self.p1.getRoiPlot().removeItem(x)
            self.firingplotlist = []
        time1, red = self.p1.roiCurve1.getData()
        time2, blue = self.p1.roiCurve2.getData()
#        blue = savgol_filter(blue,5, 3)
        redstart, redend, redlimit = pfz.findpeaks(red,redThresh)
        bluestart, blueend, bluelimit = pfz.findpeaks(blue,blueThresh)  
        redstart,redend,bluestart,blueend = pfz.cleanpeaks(red,blue,redstart,
                                redend,bluestart,blueend,redlimit,bluelimit)
        
        df = pfz.findSubPeaks(red,blue,redstart,redend,bluestart,blueend, 
                              redsubThresh,bluesubThresh)

        for i,x in enumerate(df.ident):
            if x == 0:
                color = 'r'
                data = red
            else:
                color = 'b'
                data = blue
            ep = int(df.etimes[i])
            sp = int(df.stimes[i])
#            curve = pg.PlotDataItem(x = np.linspace(sp,ep, ep-sp+1), y = data[sp:ep+1],pen = color)
#            self.firingplotlist.append(curve)
#            self.p1.getRoiPlot().addItem(curve)
            try:
                if pd.isnull(df.peaks[i]):
                    substack = np.median(self.stack[sp:ep],0)
                    if x == 0:
                        cscore = np.dot(substack.ravel(),self.czpro[0].ravel())
                        ascore = np.dot(substack.ravel(),self.czpro[1].ravel())
                        if cscore > ascore:
                            call = 'C'
                            curve = pg.PlotDataItem(x = np.linspace(sp,ep, ep-sp+1), y = data[sp:ep+1],pen = 'r',autoDownsample = True)

                        else:
                            call = 'A'
                            curve = pg.PlotDataItem(x = np.linspace(sp,ep, ep-sp+1), y = data[sp:ep+1],pen = (255,165,0),autoDownsample = True)                                                                                            
                    if x == 1:
                        gscore = np.dot(substack.ravel(),self.czpro[2].ravel())
                        tscore = np.dot(substack.ravel(),self.czpro[3].ravel())
                        if gscore > tscore:
                            call = 'G'
                            curve = pg.PlotDataItem(x = np.linspace(sp,ep, ep-sp+1), y = data[sp:ep+1],pen = 'g',autoDownsample = True)                                                                                            
                        else:
                            call = 'T'
                            curve = pg.PlotDataItem(x = np.linspace(sp,ep, ep-sp+1), y = data[sp:ep+1],pen = 'b',autoDownsample = True)                                                                                            
                    seqdf = seqdf.append(pd.DataFrame({'base':[call],'times':[sp]}),
                             ignore_index = True)
                    self.firingplotlist.append(curve)
#                    self.p1.getRoiPlot().addItem(curve)
                else:
#                    point = pg.PlotDataItem(df.mins[i], pen = None, symbol = 'o', symbolBrush = color)
#                    self.p1.getRoiPlot().addItem(point)
#                    self.firingplotlist.append(point) 
                    for j in range(len(df.mins[i])+1):
                            if j == 0:
                                ssp = sp
                                sep = int(df.mins[i][0][0])
                            elif j == len(df.mins[i]-1):
                                ssp = int(df.mins[i][-1][0])
                                sep = ep
                            else:
                                ssp = int(df.mins[i][j-1][0])
                                sep = int(df.mins[i][j][0])
                            substack = np.median(self.stack[ssp:sep+1],0)
                            if x == 0:
                                cscore = np.dot(substack.ravel(),self.czpro[0].ravel())
                                ascore = np.dot(substack.ravel(),self.czpro[1].ravel())
                                if cscore > ascore:
                                    call = 'C'
                                    curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'r',autoDownsample = True)                                    
                                else:
                                    call = 'A'
                                    curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = (255,165,0),autoDownsample = True)                                                                                                    
                            if x == 1:
                                gscore = np.dot(substack.ravel(),self.czpro[2].ravel())
                                tscore = np.dot(substack.ravel(),self.czpro[3].ravel())
                                if gscore > tscore:
                                    call = 'G'
                                    curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'g',autoDownsample = True)                                                                                            
                                else:
                                    call = 'T'
                                    curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'b',autoDownsample = True)                                                                                            
                                    
                            seqdf = seqdf.append(pd.DataFrame({'base':[call],'times':[ssp]}),
                                     ignore_index = True)
                            self.firingplotlist.append(curve)
#                            self.p1.getRoiPlot().addItem(curve)
            except ValueError:
#                point = pg.PlotDataItem(df.mins[i], pen = None, symbol = 'o', symbolBrush = color)
#                self.p1.getRoiPlot().addItem(point)
#                self.firingplotlist.append(point) 
                for j in range(len(df.mins[i])+1):
                        if j == 0:
                            ssp = sp
                            sep = int(df.mins[i][0][0])
                        elif j == len(df.mins[i]-1):
                            ssp = int(df.mins[i][-1][0])
                            sep = ep
                        else:
                            ssp = int(df.mins[i][j-1][0])
                            sep = int(df.mins[i][j][0])
                        substack = np.median(self.stack[ssp:sep+1],0)
                        if x == 0:
                            cscore = np.dot(substack.ravel(),self.czpro[0].ravel())
                            ascore = np.dot(substack.ravel(),self.czpro[1].ravel())
                            if cscore > ascore:
                                call = 'C'
                                curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'r',autoDownsample = True)                                
                            else:
                                call = 'A'
                                curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = (255,165,0),autoDownsample = True)                                                                
                        if x == 1:
                            gscore = np.dot(substack.ravel(),self.czpro[2].ravel())
                            tscore = np.dot(substack.ravel(),self.czpro[3].ravel())
                            if gscore > tscore:
                                call = 'G'
                                curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'g',autoDownsample = True)                                                                                            
                            else:
                                call = 'T'
                                curve = pg.PlotDataItem(x = np.linspace(ssp,sep, sep-ssp+1), y = data[ssp:sep+1],pen = 'b',autoDownsample = True)                                                                                            
                        seqdf = seqdf.append(pd.DataFrame({'base':[call],'times':[ssp]}),
                                 ignore_index = True) 
                        self.firingplotlist.append(curve)
#                        self.p1.getRoiPlot().addItem(curve)                        
                        
        predictedseq = seqdf.base.str.cat()
        print(predictedseq)

        fn = self.datafilename[:-3] + '_seq.fasta'
        predictedseq = Seq.Seq(predictedseq, generic_dna)
        predictedseq = SeqRecord(predictedseq, id = os.path.split(fn)[1],
                                 description = str((blueThresh,redThresh,bluesubThresh, redsubThresh)))
        SeqIO.write(predictedseq, fn, "fasta")
        print(seqdf)
        
        for i in self.firingplotlist:
            self.p1.getRoiPlot().addItem(i)    
        
        
def dubgauss(x, *p):
    A, mu1, mu2, sigma = p
    return A*np.exp(-(x-mu1)**2/(2.*sigma**2))+A*np.exp(-(x-mu2)**2/(2.*sigma**2))
    
def gauss(x, *p):
    A, mu1, sigma = p
    return A*np.exp(-(x-mu1)**2/(2.*sigma**2))
                
                
class filterpopup(QtGui.QWidget):
    acceptsignal = QtCore.pyqtSignal(int,int,int,int)
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        
        self.ui.cancelButton.clicked.connect(self.close)
        self.ui.okButton.clicked.connect(self.accept)

        
    def accept(self):
        xfilt = int(self.ui.xRangeEntry.text())
        yfilt = int(self.ui.yRangeEntry.text())
        zfilt = int(self.ui.zRangeEntry.text())
        filtertype = self.ui.filterBox.currentIndex()
        self.hide()
        self.acceptsignal.emit(zfilt,xfilt,yfilt,filtertype)
        
        
    def close(self):
        self.hide()

class analysispopup(QtGui.QWidget):
    acceptsignal = QtCore.pyqtSignal(float,float,float,float)
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.ui = analysisUi_Dialog()
        self.ui.setupUi(self)
        
        self.ui.cancelButton.clicked.connect(self.close)
        self.ui.okButton.clicked.connect(self.accept)

        
    def accept(self):
        blueThresh = float(self.ui.blueThresh.text())
        redThresh = float(self.ui.redThresh.text())
        bluePeakThresh = float(self.ui.bluePeakThresh.text())
        redPeakThresh = float(self.ui.redPeakThresh.text())
        
        self.hide()
        self.acceptsignal.emit(blueThresh,redThresh,
                               bluePeakThresh,redPeakThresh)
        
        
    def close(self):
        self.hide()
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())