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
import pyqtgraph as pg
from operator import itemgetter
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


class GUIForm(QtGui.QMainWindow):
    
    def __init__(self, master=None):
        QtGui.QMainWindow.__init__(self,master)
        self.ui = Ui_ZmwAnalysisWidget()
        self.ui.setupUi(self)
        self.filterdialog = filterpopup()
        
#        self.plotlist = [[],[],[],[],[]]
        self.seqplotlist = [[]] 


        self.ui.actionLoad.triggered.connect(self.getfile)
        self.ui.actionZ_Project.triggered.connect(self.zproject)
        self.ui.actionView_Stack.triggered.connect(self.viewstack)
        self.ui.actionFilter.triggered.connect(self.showfilterdialog)
        self.ui.actionCrop.triggered.connect(self.crop)
        self.ui.actionBase_Call.triggered.connect(self.analyze)
        self.ui.actionCheck_Controls.triggered.connect(self.checkcontrols)
        self.filterdialog.acceptsignal.connect(self.filterstack)

        
        self.p1 = self.ui.imagePlot
        self.p1.roi.scaleSnap = self.p1.roi.translateSnap = True
        self.p1.roi.setSize([4,13])
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

#        if self.plotlist != [[],[],[],[],[]]  :
#            for x in self.plotlist:
#                self.p1.getRoiPlot().removeItem(x)
#                self.plotlist = [[],[],[],[],[]]
#                self.p1.getRoiPlot().autoRange()

        if self.seqplotlist != [[]]  :
            for x in self.seqplotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.p1.getRoiPlot().autoRange()
                
#        self.plotlist = [[],[],[],[]]
        self.seqplotlist = [[]]
#        self.roip1 = self.roip2 = []
#        self.xfilt = self.yfilt =  self.zfilt = self.filtertype = []
       
        self.stack = h5py.File(self.datafilename)
        self.stack = np.array(self.stack["images"]).astype(float)
        if self.roip1 == []:
            self.originalstack = self.stack
        else:
            self.stack = self.stack[:,self.roip1[1]:self.roip2[1]+1,
                                    self.roip1[0]:self.roip2[0]+1]           
            
#        self.stack -= float(self.stack.mean())
#        self.stack /= float(self.stack.std())
        avg = self.stack.mean((1,2))
        self.background = self.stack[avg < avg.mean()-.5*avg.std()]
        self.background = np.median(self.background,0)
        self.stack -= self.background
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        
        if self.xfilt != []:
            self.filterstack(self.xfilt,self.yfilt,self.zfilt,self.filtertype)
            self.analyze()

    def zproject(self):
        self.p1.setImage(np.transpose(256*(1+self.zpro),(1,0)))
#        item = pg.ImageItem(np.transpose(256*(1+self.zpro),(1,0)))
#        item.setLookupTable(self.table)
#        self.p1.addItem(item)

    def viewstack(self):
#        self.p1.clear()
        self.p1.setImage(np.transpose(40*(2+self.stack), (0,2,1)))
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

    def closedialog(self):
        self.filterdialog.destroy()
        
    def crop(self):
        self.roip1,self.roip2 =  np.array(self.p1.roi.pos()).astype(int), np.array(self.p1.roi.pos()).astype(int)+ np.array(self.p1.roi.size()).astype(int)
        self.stack = self.originalstack[:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]
        self.p1.roi.setPos(0,0)
        self.p1.roi.setSize((int(self.p1.roi.size()[0]),int(self.p1.roi.size()[1])))
        avg = self.stack.mean((1,2))
        self.background = self.stack[avg < avg.mean()-.5*avg.std()]
        self.background = np.median(self.background,0)
        self.stack -= self.background
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
#        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        self.zpro = np.std(self.stack,0)
        self.viewstack()
        
    def analyze(self):
#        if self.plotlist != [[],[],[],[],[]]  :
#            for x in self.plotlist:
#                self.p1.getRoiPlot().removeItem(x)
#                self.plotlist = [[],[],[],[],[]] 

        if self.seqplotlist != [[]]  :
            for x in self.seqplotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.p1.getRoiPlot().autoRange()
                self.seqplotlist = [[]] 

        dntpnames = ["dCTP","dATP","dGTP","dTTP"]
        dntps = [[],[],[],[]]
        zpro = [[],[],[],[]]
        n = len(dntps)
        
        lz,lx,ly =  self.stack.shape
        seq = self.stack.reshape(lz,lx*ly)

        for i,x in enumerate(dntpnames):
            fn = [filename for filename in os.listdir(self.direc) if filename.startswith(x)
                    and filename.endswith('.h5')]
            hfile = h5py.File(self.direc + os.sep + fn[-1])
            dntps[i] = np.array(hfile["images"]).astype(float)
            if self.roip1 != []:
                dntps[i] = dntps[i][:,self.roip1[1]:self.roip2[1],
                                    self.roip1[0]:self.roip2[0]]
                                    
                                    
#            if self.filtertype != []:
    
            dntps[i] = dntps[i][-1500:]
            zpro[i] = np.median(dntps[i],0)
 
        
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
            self.czpro[i] = np.median(dntps[i],0)
            self.czpro[i] = ird.transform_img(self.czpro[i], tvec=[shifty,shiftx])
#            czpro[i] -= self.background
            self.czpro[i] = self.czpro[i]/self.czpro[i].std()
            self.czpro[i] -= self.czpro[i].mean()
    
        predictedseq = self.peakdetection()
        predictedseq = predictedseq.str.cat()
#        cscore = ascore = gscore = tscore = np.array(())
#        scores = [cscore,ascore,gscore,tscore]
#        total =[]
#
#        for i,x in enumerate(scores):
#            x = np.dot(seq,czpro[i].ravel())
#            x = (x - x.mean())/x.std()
#            if i==0:
#                total = x**2
#            else:       
#                total = total+x**2
#            
#            color = pg.intColor(i, hues = 4)
#            scores[i] = x
#            
#            self.plotlist[i] = self.p1.getRoiPlot().plot(x, pen = color)
#        self.p1.getRoiPlot().setRange(yRange=[-4,10])
#            
#        scores = pd.DataFrame(np.transpose(scores))
##        tot = scores[0]+scores[1] + scores[2] + scores[3]
##        tot = self.stack.mean((1,2))
#
#        predictedseq = np.zeros(len(scores[0]))
#        predictedseq[np.array(scores.idxmax(axis=1)==0)] = 1
#        predictedseq[np.array(scores.idxmax(axis=1)==1)] = 2
#        predictedseq[np.array(scores.idxmax(axis=1)==2)] = 3
#        predictedseq[np.array(scores.idxmax(axis=1)==3)] = 4
##        
##        baseline = tot[tot < tot.mean()]
##        baseline = np.median(baseline)
##        predictedseq[np.array(tot < baseline)] = 0         
##        predictedseq[scores[0]<1 and scores[1]<1 and scores[2]<1 and scores[3]<1] = 0
#
#        scores = scores.T
#        predictedseq[np.array((scores<1).all())] = 0
#        
#        for i,x in enumerate(predictedseq):
#            if i != len(predictedseq)-1 and i !=0:
#                if x != 0 and (x != predictedseq[i-1] or x != predictedseq[i+1]):
#                    if x == 0:
#                        predictedseq[i] = 0
#                    elif predictedseq[i-1] ==predictedseq[i+1]:
#                        predictedseq[i] = predictedseq[i-1]
#                     
#        self.plotlist[4] = self.p1.getRoiPlot().plot(predictedseq, pen = pg.mkPen('w', width = 3))
#        
#        deletearray = np.zeros(len(predictedseq))
#
#        for i,x in enumerate(predictedseq):
#            if i != len(predictedseq)-1:
#                if x != 0 and (x != predictedseq[i-1] or x != predictedseq[i+1]):
#                    deletearray[i] = 1
#        
#        predictedseq = predictedseq[deletearray!=1]
#        deletearray = np.zeros(len(predictedseq))
#        
#        for i,x in enumerate(predictedseq):
#            if i != len(predictedseq)-1:
#                if predictedseq[i] == predictedseq[i+1]:
#                    deletearray[i] = 1
#        
#        predictedseq = predictedseq[deletearray!=1]
#        predictedseq = predictedseq[predictedseq!=0]
#        
#        predictedseq = predictedseq.astype('str')
#        
#        for i,x in enumerate(predictedseq):
#            if x == '1.0':
#                predictedseq[i] = 'C'
#            if x == '2.0':
#                predictedseq[i] = 'A'
#            if x == '3.0':
#                predictedseq[i] = 'G'
#            if x == '4.0':
#                predictedseq[i] = 'T'
#                
#        predictedseq = "".join(predictedseq)
#        print(predictedseq)
#                
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
        
    def peakdetection(self):
        cscore = ascore = gscore = tscore = spseries = pd.Series()
        markers = ['C','A','G','T']
        time, intensity = self.p1.roiCurve1.getData()
        noise = np.std(intensity[intensity < np.mean(intensity)])
        thresh = np.median(intensity) + 2*noise
        
        firing = np.where(intensity > thresh)[0]
        startandend = np.diff(firing)
        startpoints = np.insert(startandend, 0, 2)
        endpoints = np.insert(startandend, -1, 2)
        startpoints = np.where(startpoints>1)[0]
        endpoints = np.where(endpoints>1)[0]
        startpoints = firing[startpoints]
        endpoints = firing[endpoints]
        
        peakseries = pd.Series()
        minseries = pd.Series()
        for i,x in enumerate(startpoints):
            sp = startpoints[i] 
            ep = endpoints[i]+1
#            self.p1.getRoiPlot().plot(x = np.arange(sp, ep), y = intensity[sp:ep], pen = 'b')
            try:
                peaks, mins = peakdet(v = intensity[sp:ep],delta = noise, x = np.arange(sp,ep))
                if len(peaks) == 0:
                    peaks = np.NAN
                if len(mins) == 0:
                    peaks = np.NAN
#                else:
#                    self.p1.getRoiPlot().plot(peaks, pen = None, symbol = 'o', symbolBrush = 'g')
#                    self.p1.getRoiPlot().plot(mins, pen = None, symbol = 'o', symbolBrush = 'r')
                peakseries = peakseries.append(pd.Series([peaks]))
                minseries = minseries.append(pd.Series([peaks]))
        #        p1.plot(df[i], pen = None, symbol = 'o', symbolBrush = 'b')
            except:
                ValueError
                
            if np.size(peaks) < 4:
                substack = np.mean(self.stack[sp:ep],0)
                cscore = cscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[0].ravel())),ignore_index=True) 
                ascore = ascore.append(pd.Series(np.dot(substack.ravel(),self.czpro[1].ravel())),ignore_index=True)
                gscore = gscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[2].ravel())),ignore_index=True)
                tscore = tscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[3].ravel())),ignore_index=True)
                spseries = spseries.append(pd.Series(sp),ignore_index = True)
            else:
                for j,y in enumerate(peaks):
                    if j == 0:
                        ssp = sp
                        sep = int(mins[j][0])
                    elif j == len(peaks)-1:
                        ssp = int(mins[j-1][0])
                        sep = ep
                    else:
                        ssp = int(mins[j-1][0])
                        sep = int(mins[j][0])
                    substack = np.mean(self.stack[ssp:sep+1],0)
                    cscore = cscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[0].ravel())),ignore_index=True) 
                    ascore = ascore.append(pd.Series(np.dot(substack.ravel(),self.czpro[1].ravel())),ignore_index=True)
                    gscore = gscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[2].ravel())),ignore_index=True)
                    tscore = tscore.append(pd.Series(np.dot(substack.ravel(),self.czpro[3].ravel())),ignore_index=True)
                    spseries = spseries.append(pd.Series(ssp),ignore_index = True)
                
        scoredf = pd.DataFrame({'C':cscore,'A':ascore,'G':gscore,'T':tscore,'sp':spseries})
        sequence = scoredf[['A','C','G','T']].idxmax(axis=1)
        print(sequence.str.cat())
        for i,x in enumerate(sequence):
            if x == 'C':
                color = 'r'
            if x == 'A':
                color = 'y'
            if x == 'G':
                color = 'g'
            if x == 'T':
                color = 'b'
            text = pg.TextItem(x,color = color)
            seqplot = self.p1.getRoiPlot().addItem(text)
            if i == 0:
                self.seqplotlist = [text]
            else:
                self.seqplotlist.append(text)
            text.setPos(spseries[i],intensity[spseries[i]])
        
        return sequence
                    
        
        
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
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())