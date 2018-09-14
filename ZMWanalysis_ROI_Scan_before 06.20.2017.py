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
import h5py
from scipy.ndimage import filters
import pandas as pd
import peakFindingZMW as pfz



class GUIForm(QtWidgets.QMainWindow):
    fileSelectSignal = QtCore.pyqtSignal(str)     

    def __init__(self, master=None):

        QtWidgets.QMainWindow.__init__(self,master)
        self.ui = Ui_ZmwAnalysisWidget()
        self.ui.setupUi(self)
        self.filterdialog = filterpopup()
        self.dataProcessThread = dataProcesser()

        self.ui.actionLoad.triggered.connect(self.getfile)
        self.ui.actionZ_Project.triggered.connect(self.zproject)
        self.ui.actionView_Stack.triggered.connect(self.viewstack)
        self.ui.actionFilter.triggered.connect(self.showfilterdialog)
        self.ui.actionCrop.triggered.connect(self.crop)
        self.ui.actionBase_Call.triggered.connect(self.analyze)

        self.filterdialog.acceptsignal.connect(self.filterstack)
        self.dataProcessThread.dataLoadSignal.connect(self.Load)
        self.fileSelectSignal.connect(self.dataProcessThread.loadFile)

        self.p1 = self.ui.imagePlot
#        self.p1.roi.scaleSnap = self.p1.roi.translateSnap = True
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

    def getfile(self):

        if self.direc == []:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',os.getcwd(),("*.h5"))[0]
            self.direc = os.path.dirname(str(self.datafilename))
            print(self.datafilename)
        else:
            self.datafilename = QtWidgets.QFileDialog.getOpenFileName(self,'Open file',self.direc,("*.h5"))[0]
            self.direc=os.path.dirname(str(self.datafilename))
            print(self.datafilename)
            
        QtGui.QApplication.processEvents()
        self.fileSelectSignal.emit(self.datafilename)

    def Load(self,stack):

        self.stack = stack

        avg = self.stack.mean((1,2))
        self.background = self.stack[avg < avg.mean()-.5*avg.std()]
        self.background = np.median(self.background,0)
        self.stack -= self.background
        self.stack -= np.mean(self.stack)
        self.stack = self.stack/self.stack.std()
        self.stack += 100
        self.zpro = np.std(self.stack,0)
        self.zproject()
        

    def zproject(self):
        self.p1.setImage(np.transpose(256*(1+self.zpro),(1,0)))

    def viewstack(self):
        self.p1.setImage(np.transpose(40*(2+self.stack), (0,2,1)))
        self.p1.play(rate=60)
        
    def showfilterdialog(self):
        self.filterdialog.show()
        
    def filterstack(self,xfilt,yfilt,zfilt,filtertype):
        self.p1.play(0)
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
        self.zpro = np.std(self.stack,0)
        self.zproject()
        self.zpro = np.std(self.stack,0)
        self.viewstack()
        
    def analyze(self):
        if self.p1.image.ndim == 2:
            self.viewstack()
        self.p1.play(0)
        df = pd.DataFrame()
        peakCount = []
        
        roiSPosx, roiSPosy = self.p1.roi.pos()
        rows = 22
        columns = 8
        with pg.ProgressDialog("Scanning...", 0, rows*columns) as dlg:        
            for i in range(columns):
                for j in range(rows):
                    time, data = self.p1.roiCurve1.getData()
                    df[str((i,j))] = pd.Series(data)                
                    roiPosx, roiPosy = self.p1.roi.pos()
                    roiPosx += 7.5
                    self.p1.roi.setPos([roiPosx,roiPosy])
                    dlg += 1
                roiPosx = roiSPosx
                roiPosy += 20
                self.p1.roi.setPos([roiPosx,roiPosy])
    
            self.win = pg.GraphicsWindow(title="pyqtgraph example: Linked Views")
            self.win.resize(800,800)
            self.win.scene().sigMouseClicked.connect(self.plotzoom)
            
            plotlist = [[]]
            count = 0

        with pg.ProgressDialog("Analyzing...", 0, 20*8) as dlg:               
            for i,x in enumerate(df.keys()):
                if i%10 == 0:
                    self.win.nextRow()
                plotTemp = self.win.addPlot()
                plotTemp.plot(df[x],pen = pg.intColor(i, 100))
                plotTemp.hideAxis('bottom')
                plotTemp.hideAxis('left')
                plotlist.append(plotTemp)
                
                try:
                    datastart, dataend, datalimit = pfz.findpeaks(df[x],12)
                    count += len(datastart)
                    peakCount.append(len(datastart))
                    point = pg.PlotDataItem(x = datastart, y = df[x][datastart], 
                                            pen = None, symbol = 'o', 
                                            symbolBrush = pg.intColor(i, 100),
                                            symbolSize=6)
                    plotTemp.addItem(point)
                except IndexError:
                    pass
                except ValueError:
                    pass
                dlg += 1                    
                
            
        self.p1.roi.setPos([roiSPosx,roiSPosy])
        print('Events per file per ZMW = ' + str(count/(len(df.keys()))))
        filebase = str(os.path.splitext(self.datafilename)[0])
        np.savetxt(filebase + 'peakCount.txt', peakCount, delimiter = ',')
    
    def plotzoom(self, event):
        pass
        items = self.win.scene().itemsNearEvent(event)
        p1 = pg.plot()
        for i in items:
            print(str(type(i)))
            if str(type(i)) == "<class 'pyqtgraph.graphicsItems.PlotCurveItem.PlotCurveItem'>" or str(type(i)) == "<class 'pyqtgraph.graphicsItems.ScatterPlotItem.ScatterPlotItem'>":
                zoomitem = i
                p1.addItem(zoomitem)
                
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

class dataProcesser(QtCore.QThread):

    dataLoadSignal = QtCore.pyqtSignal(object) 

    def __init__(self):
        QtCore.QThread.__init__(self, parent=None)        

    def loadFile(self,fn):
        stack = h5py.File(fn)
        stack = np.array(stack["images"]).astype(float)
        self.dataLoadSignal.emit(stack)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())
