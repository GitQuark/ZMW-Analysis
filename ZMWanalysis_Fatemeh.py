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
import peakFindingZMW


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
        self.ui.actionCheck_Controls.triggered.connect(self.checkcontrols)
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
                self.p1.getRoiPlot().autoRange()
                
        if self.firingplotlist != []  :
            for x in self.firingplotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.p1.getRoiPlot().autoRange()           

        self.seqplotlist = [[]]
        self.firingplotlist = []
       
        self.stack = h5py.File(self.datafilename)
        self.stack = np.array(self.stack["images"]).astype(float)
#        self.originalstack = self.stack
        if self.roip1 == []:
            self.originalstack = self.stack
        else:
            self.stack = self.stack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]           
            
#        self.stack -= float(self.stack.mean())
#        self.stack /= float(self.stack.std())
#        baseline = np.median(self.stack[self.stack < np.median(self.stack)])
#        self.background = self.stack[np.max(self.stack,(1,2)) < baseline]
        self.stack = self.backgroundsubtract(self.stack)
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
#        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        
        if self.xfilt != []:
            self.filterstack(self.xfilt,self.yfilt,self.zfilt,self.filtertype)
            self.analyze(self.blueThresh,self.redThresh,
                               self.bluePeakThresh,self.redPeakThresh)

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
        self.stack = self.originalstack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]
        self.p1.roi.setPos(0,0)
        self.p1.roi.setSize((int(self.p1.roi.size()[0]),int(self.p1.roi.size()[1])))
        self.stack = self.backgroundsubtract(self.stack)
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
#        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        self.zpro = np.std(self.stack,0)
        self.viewstack()
        
    def showanalysisdialog(self):
        self.analysisdialog.show()
        
        
    def analyze(self, blueThresh,redThresh,bluePeakThresh,redPeakThresh):
        self.blueThresh = blueThresh
        self.redThresh = redThresh
        self.bluePeakThresh = bluePeakThresh
        self.redPeakThresh = redPeakThresh
        if self.seqplotlist != [[]]  :
            for x in self.seqplotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.p1.getRoiPlot().autoRange()
                
        if self.firingplotlist != []  :
            for x in self.firingplotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.p1.getRoiPlot().autoRange()           

        self.seqplotlist = [[]]
        self.firingplotlist = []
        
        dntpnames = ["dCTP","dATP","dGTP","dTTP"]
        cdf = pd.DataFrame({'a':[],'t':[],'g':[],'c':[]})
        dntps = [[],[],[],[]]
        zpro = [[],[],[],[]]
        n = len(dntps)
        
        dntpdirec = r'C:/Users/Meni/Desktop/Fatemeh_dNTPs'
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
    
        seqdf = self.peakdetection(blueThresh,redThresh,bluePeakThresh,redPeakThresh)
        predictedseq = seqdf.base.str.cat()

        fn = self.datafilename[:-3] + '_seq.fasta'
        predictedseq = Seq.Seq(predictedseq, generic_dna)
        predictedseq = SeqRecord(predictedseq, id = os.path.split(fn)[1])
        SeqIO.write(predictedseq, fn, "fasta")


    def checkcontrols(self):
        dntpnames = ["dCTP","dATP","dGTP","dTTP"]
        dntps = [[],[],[],[]]
        zpro = [[],[],[],[]]
        czpro = [[],[],[],[]]
        n = len(dntps)
        
        for i,x in enumerate(dntpnames):
            fn = [filename for filename in os.listdir(self.direc) if filename.startswith(x)
                    and filename.endswith('.h5')]
            hfile = h5py.File(self.direc + os.sep + fn[-1])
#            print self.direc + os.sep+ fn[-1]
            dntps[i] = np.array(hfile["images"]).astype(float)
            if self.roip1 != []:
                dntps[i] = dntps[i][:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]
#            dntps[i] -= self.background

#            if self.filtertype != []:
    
            dntps[i] = dntps[i][-1500:]
#            dntps[i] = (dntps[i]/dntps[i].std((1,2))[0])
#            dntps[i] -= dntps[i].mean()
#            dntps[i] += self.background.mean()
            zpro[i] = np.median(dntps[i],0)
#            zpro[i] = np.ma.masked_where(zpro[i] < 0, zpro[i])
 
#        zpro[0] -= np.mean(zpro[0])
#        zpro[1] -= np.mean(zpro[1])
#        zpro[2] -= np.mean(zpro[2])
#        zpro[3] -= np.mean(zpro[3])

        
        fig = plt.figure()
        plt.subplot(1,3,1)   
        cmap = plt.cm.get_cmap("Reds")
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        im1 = plt.imshow(zpro[0], cmap = my_cmap)
        
        cmap = plt.cm.get_cmap("Greens")
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        im2 = plt.imshow(zpro[1], cmap = my_cmap)
        
        cmap = plt.cm.get_cmap("Blues")
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        plt.imshow(zpro[2], cmap = my_cmap)
        
        cmap = plt.cm.get_cmap("PuRd")
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        im4 = plt.imshow(zpro[3], cmap = my_cmap)
        fig.show()
        
        plt.subplot(1,3,2)   
        plt.imshow(self.zpro)
        
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
        
        for i,x in enumerate(dntps):
            czpro[i] = np.median(dntps[i],0)
            czpro[i] = ird.transform_img(czpro[i], tvec=[shifty,shiftx])
            czpro[i] -= self.background
            czpro[i] = czpro[i]/czpro[i].std()
            czpro[i] -= czpro[i].mean()
        
        ccomposite = czpro[0] + czpro[1] + czpro[2] + czpro[3]
        plt.subplot(1,3,3)   
        plt.imshow(ccomposite)
        
    def peakdetection(self,blueThresh,redThresh,bluePeakThresh,redPeakThresh):
        markers = ['C','A','G','T']
        colors = ['r','b']
        score = ascore = gscore = tscore = spseries = pd.Series()
        df = pd.DataFrame({"ident":[],'stimes':[],'etimes':[],
                           'peaks':[],'mins':[]})
        seqdf = pd.DataFrame({'base':[],'times':[]})
        time1, intensity1 = self.p1.roiCurve1.getData()
        time2, intensity2 = self.p1.roiCurve2.getData()
        intensities = [intensity1,intensity2]
        self.firingplotlist = []
        peakseries = pd.Series()
        minseries = pd.Series()
        for ind,intensity in enumerate(intensities):
            if ind == 0:
                thresh = redThresh
                noise = redPeakThresh                
            if ind ==1:
                thresh = blueThresh
                noise = bluePeakThresh
            
            firing = np.where(((intensity > thresh) & (intensity > intensities[ind-1])))[0]
            startandend = np.diff(firing)
            startpoints = np.insert(startandend, 0, 2)
            endpoints = np.insert(startandend, -1, 2)
            startpoints = np.where(startpoints>1)[0]
            endpoints = np.where(endpoints>1)[0]
            startpoints = firing[startpoints]
            endpoints = firing[endpoints]
            df = df.append(pd.DataFrame({"ident":[ind]*len(startpoints),
            'stimes':startpoints,'etimes':endpoints}),ignore_index=True)
    
            for i,x in enumerate(startpoints):
                sp = startpoints[i] 
                ep = endpoints[i]+1
                curve = pg.PlotDataItem(x = np.linspace(sp, ep,ep-sp),
                                   y = intensity[sp:ep], 
                                   pen = pg.mkPen(colors[ind], width = 2))
                self.firingplotlist.append(curve)
                
                self.p1.getRoiPlot().addItem(curve)
                try:
                    peaks, mins = peakdet(v = intensity[sp:ep],delta = noise, x = np.arange(sp,ep))
                    if len(peaks) == 0 or len(mins) == 0:
                        peaks = np.NAN
                        substack = np.mean(self.stack[sp:ep],0)
                        if ind == 0:
                            cscore = np.dot(substack.ravel(),self.czpro[0].ravel())
                            ascore = np.dot(substack.ravel(),self.czpro[1].ravel())
                            if cscore > ascore:
                                call = 'C'
                            else:
                                call = 'A'
                        if ind == 1:
                            gscore = np.dot(substack.ravel(),self.czpro[2].ravel())
                            tscore = np.dot(substack.ravel(),self.czpro[3].ravel())
                            if gscore > tscore:
                                call = 'G'
                            else:
                                call = 'T'
                        seqdf = seqdf.append(pd.DataFrame({'base':[call],'times':[sp]}),
                                 ignore_index = True)                         
                    else:
#                        point = pg.PlotDataItem(peaks, pen = None, symbol = 'o', symbolBrush = 'g')
#                        self.p1.getRoiPlot().addItem(point)
#                        self.firingplotlist.append(point)
#                        point = pg.PlotDataItem(mins, pen = None, symbol = 'o', symbolBrush = 'r')
#                        self.p1.getRoiPlot().addItem(point)
#                        self.firingplotlist.append(point)
                        for i,x in enumerate(peaks):
                            if i == 0:
                                ssp = sp
                                sep = int(mins[i][0])
                            elif i == len(peaks)-1:
                                ssp = int(mins[i-1][0])
                                sep = ep
                            else:
                                ssp = int(mins[i-1][0])
                                sep = int(mins[i][0])
                            substack = np.mean(self.stack[ssp:sep+1],0)
                            if ind == 0:
                                cscore = np.dot(substack.ravel(),self.czpro[0].ravel())
                                ascore = np.dot(substack.ravel(),self.czpro[1].ravel())
                                if cscore > ascore:
                                    call = 'C'
                                else:
                                    call = 'A'
                            if ind == 1:
                                gscore = np.dot(substack.ravel(),self.czpro[2].ravel())
                                tscore = np.dot(substack.ravel(),self.czpro[3].ravel())
                                if gscore > tscore:
                                    call = 'G'
                                else:
                                    call = 'T'
                            seqdf = seqdf.append(pd.DataFrame({'base':[call],'times':[ssp]}),
                                     ignore_index = True)                        
                    peakseries = peakseries.append(pd.Series([peaks]))
                    minseries = minseries.append(pd.Series([mins]))
                except:
                    ValueError

        seqdf = seqdf.sort(['times','base'])
        for i,x in enumerate(seqdf.index):
            base = seqdf.base[i]
            if base == 'C':
                color = 'r'
                intensity = intensities[0][int(seqdf.times[i])]
            if base == 'A':
                color = 'y'
                intensity = intensities[0][int(seqdf.times[i])]
            if base == 'G':
                color = 'g'
                intensity = intensities[1][int(seqdf.times[i])]
            if base == 'T':
                color = 'b'
                intensity = intensities[1][int(seqdf.times[i])]
            text = pg.TextItem(base,color = color)
            seqplot = self.p1.getRoiPlot().addItem(text)
            if i == 0:
                self.seqplotlist = [text]
            else:
                self.seqplotlist.append(text)
            text.setPos(seqdf.times[i],intensity)
        
        print(seqdf.base.str.cat())
        return seqdf
                    
        
        
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