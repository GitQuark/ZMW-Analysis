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
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image, ImageSequence


class GUIForm(QtGui.QMainWindow):
    
    def __init__(self, master=None):
        QtGui.QMainWindow.__init__(self,master)
        self.ui = Ui_ZmwAnalysisWidget()
        self.ui.setupUi(self)
        self.filterdialog = filterpopup()
        
        self.plotlist = [[],[],[],[],[]]


        self.ui.actionLoad.triggered.connect(self.getfile)
        self.ui.actionZ_Project.triggered.connect(self.zproject)
        self.ui.actionView_Stack.triggered.connect(self.viewstack)
        self.ui.actionFilter.triggered.connect(self.showfilterdialog)
        self.ui.actionCrop.triggered.connect(self.crop)
        self.ui.actionBase_Call.triggered.connect(self.analyze)
        self.ui.actionCheck_Controls.triggered.connect(self.checkcontrols)
        self.filterdialog.acceptsignal.connect(self.filterstack)

        
        self.p1 = self.ui.imagePlot
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


    def getfile(self):

        if self.direc==[]:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',os.getcwd(),("*.tif"))[0]
            self.direc=os.path.dirname(str(self.datafilename))
            print(self.datafilename)
        else:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self,'Open file',self.direc,("*.tif"))[0]
            self.direc=os.path.dirname(str(self.datafilename))


        self.Load()

    def Load(self):

        if self.plotlist != [[],[],[],[],[]]  :
            for x in self.plotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.plotlist = [[],[],[],[],[]]
                self.p1.getRoiPlot().autoRange()
                
        self.plotlist = [[],[],[],[]]
        self.roip1 = self.roip2 = []
        self.xfilt = self.yfilt =  self.zfilt = self.filtertype = []

        if str(os.path.splitext(self.datafilename)[1])== u'.tif':
            im=Image.open(self.datafilename)
            self.imagesize=np.array(im.size)
            self.stack = []
            for frame in ImageSequence.Iterator(im):
                self.stack.append(frame)
            self.stack = [np.array(frame.getdata()) for frame in ImageSequence.Iterator(im)]
            self.stack =np.reshape(self.stack,[len(self.stack),self.imagesize[1],self.imagesize[0]])
#            self.stack = float(self.stack)
        
        if str(os.path.splitext(self.datafilename)[1])== u'.h5':
            self.stack = h5py.File(self.datafilename)
            self.stack = np.array(self.stack["images"]).astype(float)
            
        self.originalstack = self.stack
            
#        self.stack -= float(self.stack.mean())
#        self.stack /= float(self.stack.std())
        avg = self.stack.mean((1,2))
        self.background = self.stack[avg < avg.mean()-.5*avg.std()]
        self.background = np.median(self.background,0)
#        self.stack -= self.background
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()

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
        self.viewstack()

    def closedialog(self):
        self.filterdialog.destroy()
        
    def crop(self):
        self.roip1,self.roip2 =  self.p1.roi.pos(), self.p1.roi.pos() + self.p1.roi.size()
        self.stack = self.stack[:,self.roip1[1]:self.roip2[1]+1,
                                    self.roip1[0]:self.roip2[0]+1]
        self.p1.roi.setPos(0,0)
        self.p1.roi.setSize((int(self.p1.roi.size()[0]),int(self.p1.roi.size()[1])))
        self.viewstack()
        
    def analyze(self):
        if self.plotlist != [[],[],[],[],[]]  :
            for x in self.plotlist:
                self.p1.getRoiPlot().removeItem(x)
                self.plotlist = [[],[],[],[],[]] 

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
#            print self.direc + os.sep+ fn[-1]
            dntps[i] = np.array(hfile["images"]).astype(float)
            if self.roip1 != []:
                dntps[i] = dntps[i][:,self.roip1[1]:self.roip2[1]+1,
                                    self.roip1[0]:self.roip2[0]+1]

#            if self.filtertype != []:
    
            dntps[i] = dntps[i][100:900]
            dntps[i] -= np.mean(dntps[i])
            dntps[i] = (dntps[i]/dntps[i].std())
            zpro[i] = np.median(dntps[i],0)
 
        zpro[0] -= np.mean(zpro[0])
        zpro[1] -= np.mean(zpro[1])
        zpro[1] -= np.mean(zpro[2])
        zpro[2] -= np.mean(zpro[3])
    
        cscore = ascore = gscore = tscore = np.array(())
        scores = [cscore,ascore,gscore,tscore]
        total =[]

        for i,x in enumerate(scores):
            x = np.dot(seq,zpro[i].ravel())
            x = (x - x.mean())/x.std()
            if i==0:
                total = x**2
            else:       
                total = total+x**2
            
            color = pg.intColor(i, hues = 4)
            scores[i] = x
            
            self.plotlist[i] = self.p1.getRoiPlot().plot(x-np.median(x), pen = color)
        self.p1.getRoiPlot().setRange(yRange=[-4,10])
            
        scores = pd.DataFrame(np.transpose(scores))
        tot = abs(scores[0]+scores[1]-scores[2]-scores[3])

        predictedseq = np.zeros(len(tot))
        predictedseq[np.array(scores.idxmax(axis=1)==0)] = 1
        predictedseq[np.array(scores.idxmax(axis=1)==1)] = 2
        predictedseq[np.array(scores.idxmax(axis=1)==2)] = 3
        predictedseq[np.array(scores.idxmax(axis=1)==3)] = 4
        
        predictedseq[np.array(tot < np.median(predictedseq) + np.std(predictedseq))] = 0
        
        self.plotlist[4] = self.p1.getRoiPlot().plot(predictedseq, pen = pg.mkPen('w', width = 3))

    def checkcontrols(self):
        dntpnames = ["dCTP","dATP","dGTP","dTTP"]
        dntps = [[],[],[],[]]
        zpro = [[],[],[],[]]
        n = len(dntps)
        
        for i,x in enumerate(dntpnames):
            fn = [filename for filename in os.listdir(self.direc) if filename.startswith(x)
                    and filename.endswith('.h5')]
            hfile = h5py.File(self.direc + os.sep + fn[-1])
#            print self.direc + os.sep+ fn[-1]
            dntps[i] = np.array(hfile["images"]).astype(float)
            if self.roip1 != []:
                dntps[i] = dntps[i][:,self.roip1[1]:self.roip2[1]+1,
                                    self.roip1[0]:self.roip2[0]+1]

#            if self.filtertype != []:
    
            dntps[i] = dntps[i][100:900]
            dntps[i] -= np.mean(dntps[i])
            dntps[i] = (dntps[i]/dntps[i].std())
            zpro[i] = np.median(dntps[i],0)
            zpro[i] = np.ma.masked_where(zpro[i] < 0, zpro[i])
 
#        zpro[0] -= np.mean(zpro[0])
#        zpro[1] -= np.mean(zpro[1])
#        zpro[2] -= np.mean(zpro[2])
#        zpro[3] -= np.mean(zpro[3])

        
        fig = plt.figure()
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
        im3 = plt.imshow(zpro[2], cmap = my_cmap)
        
        cmap = plt.cm.get_cmap("PuRd")
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)
        im4 = plt.imshow(zpro[3], cmap = my_cmap)
        fig.show()
            
            
                
                
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