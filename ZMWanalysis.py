# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 11:55:59 2014

@author: Robert Henley
"""
import sys
import os
import numpy as np
from Utility import *
from zmwanalysiswidget import *
from filterdialog import *
from threshdialogui import *
import pyqtgraph as pg
import h5py
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


class GUIForm(QtWidgets.QMainWindow):

    def __init__(self, master=None):
        QtWidgets.QMainWindow.__init__(self, master)
        self.ui = Ui_ZmwAnalysisWidget()
        self.ui.setupUi(self)
        self.filter_dialog = FilterPopup()
        self.analysis_dialog = AnalysisPopup()

        #        self.plotlist = [[],[],[],[],[]]
        self.seqplotlist = [[]]
        self.firingplotlist = []

        self.ui.actionLoad.triggered.connect(self.getfile)
        self.ui.actionZ_Project.triggered.connect(self.z_project)
        self.ui.actionView_Stack.triggered.connect(self.view_stack)
        self.ui.actionClear_Settings.triggered.connect(self.clear_settings)
        self.ui.actionFilter.triggered.connect(self.show_filter_dialog)
        self.ui.actionCrop.triggered.connect(self.crop)
        self.ui.actionBase_Call.triggered.connect(self.show_analysis_dialog)
        self.ui.actionCheck_Controls.triggered.connect(self.check_controls)
        self.filter_dialog.acceptsignal.connect(self.filter_stack)
        self.analysis_dialog.acceptsignal.connect(self.analyze)

        self.image_plot = self.ui.imagePlot
        self.image_plot.roi.scaleSnap = self.image_plot.roi.translateSnap = True
        self.image_plot.roi.setSize([6, 17])
        colors = [
            (0, 0, 0),
            (45, 5, 61),
            (84, 42, 55),
            (150, 87, 60),
            (208, 171, 141),
            (255, 255, 255)
        ]
        cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=colors)
        self.image_plot.setColorMap(cmap)
        self.direc = []
        self.roi_image_plot = self.roip2 = []
        self.x_filter = self.yfilt = self.zfilt = self.filtertype = []
        self.xyz_filter = [[], [], []]
        self.blueThresh = None
        self.redThresh = None
        self.bluePeakThresh = None
        self.redPeakThresh = None

        self.data_file_name = None
        self.background = None
        self.stack = None
        self.zpro = None

        self.p1.roi.scaleSnap = self.p1.roi.translateSnap = True

    def getfile(self):

        if not self.direc:
            self.direc = os.getcwd()

        self.data_file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.direc, "*.h5")[0]
        self.direc = os.path.dirname(str(self.data_file_name))
        print(self.data_file_name)

        self.load(self.data_file_name, self.image_plot)

    def load(self, data_file_name, image_plot):

        if self.seqplotlist != [[]]:
            for x in self.seqplotlist:
                image_plot.getRoiPlot().removeItem(x)
                image_plot.getRoiPlot().autoRange()

        if self.firingplotlist:
            for x in self.firingplotlist:
                image_plot.getRoiPlot().removeItem(x)
                image_plot.getRoiPlot().autoRange()

        self.seqplotlist = [[]]
        self.firingplotlist = []

        stack = h5py.File(data_file_name)
        if 'images' in list(stack.keys()):
            stack = np.array(stack["images"]).astype(float)
        else:
            print('Entry "images" not in loaded h5 file')
            return

        # self.original_stack = self.stack
        if not self.roi_image_plot:
            self.original_stack = stack
        else:
            # First dimension is fraME NUMBER
            stack = stack[:, self.roi_image_plot[1]:self.roip2[1], self.roi_image_plot[0]:self.roip2[0]]

        # self.stack -= float(self.stack.mean())
        # self.stack /= float(self.stack.std())
        # baseline = np.median(self.stack[self.stack < np.median(self.stack)])
        # self.background = self.stack[np.max(self.stack,(1,2)) < baseline]
        stack = subtract_background(stack, self.background)
        stack -= np.mean(stack)
        stack = stack / stack.std()
        # self.stack += 100
        self.zpro = np.std(stack, 0)
        self.stack = stack
        z_project(self.image_plot, self.zpro)

        if self.x_filter:
            self.stack, self.zpro = filter_stack(self.stack, self.xyz_filter, self.filtertype, self.image_plot)
            self.analyze(self.blueThresh, self.redThresh,
                         self.bluePeakThresh, self.redPeakThresh)

        return stack

    def show_filter_dialog(self):
        self.filter_dialog.show()

    def view_stack(self):
        view_stack(self.stack, self.image_plot)

    def z_project(self):
        z_project(self.image_plot, self.zpro)

    def filter_stack(self):
        filter_stack(self.stack, self.xyz_filter, self.filtertype, self.image_plot)

    def clear_settings(self):
        self.roi_image_plot = []
        self.xyz_filter = [[], [], []]

    def close_dialog(self):
        self.filter_dialog.destroy()

    def crop(self):
        upper_left = self.image_plot.roi.pos().astype(int)
        roi_size = np.array(self.image_plot.roi.size()).astype(int)
        lower_right = upper_left + roi_size
        x0, y0 = upper_left
        x1, y1 = lower_right
        stack = self.original_stack[:, y0:y1, x0:x1]
        self.image_plot.roi.setPos(0, 0)
        self.image_plot.roi.setSize((int(roi_size[0]), int(roi_size[1])))
        stack = subtract_background(stack, self.background)
        stack -= np.mean(stack)
        stack = stack / stack.std()
        #        self.stack += 100
        self.zpro = np.std(stack, 0)
        z_project(self.image_plot, self.zpro)
        self.zpro = np.std(stack, 0)
        view_stack(stack, self.image_plot)
        self.stack = stack

    def show_analysis_dialog(self):
        self.analysis_dialog.show()

    def analyze(self, blue_thresh, red_thresh, blue_peak_thresh, red_peak_thresh):
        self.blueThresh = blue_thresh
        self.redThresh = red_thresh
        self.bluePeakThresh = blue_peak_thresh
        self.redPeakThresh = red_peak_thresh
        if self.seqplotlist != [[]]:
            for name in self.seqplotlist:
                self.image_plot.getRoiPlot().removeItem(name)
                self.image_plot.getRoiPlot().autoRange()

        if self.firingplotlist:
            for name in self.firingplotlist:
                self.image_plot.getRoiPlot().removeItem(name)
                self.image_plot.getRoiPlot().autoRange()

        self.seqplotlist = [[]]
        self.firingplotlist = []

        dntp_names = ["dCTP", "dATP", "dGTP", "dTTP"]
        cdf = pd.DataFrame({'a': [], 't': [], 'g': [], 'c': []})
        dntps = [[], [], [], []]
        zpro = [[], [], [], []]

        # dntpdirec = r'D:/Vivek/ZMW Data/02.27.2017 - 20kbp SMRT - looks good'
        # fn = 'dntps.h5'
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd(), '*.h5')
        if file_name == ('', ''):  # Interaction canceled
            return
        file = h5py.File(file_name)
        for idx, name in enumerate(dntp_names):
            zpro[idx] = np.array(file[name]).astype(float)
            zpro[idx] -= zpro[idx].mean()
            zpro[idx] /= zpro[idx].std()

        composite = zpro[0] + zpro[1] + zpro[2] + zpro[3]
        if len(composite.shape) == 2:
            ly, lx = composite.shape
        else:  # Should only be 3
            lz, ly, lx = composite.shape

        shift_x, shift_y = calculate_xy_shift(lx, ly, composite, zpro)
        self.zpro = zpro
        self.czpro = [[], [], [], []]
        for idx, name in enumerate(dntps):
            self.czpro[idx] = zpro[idx]
            self.czpro[idx] = ird.transform_img(self.czpro[idx], tvec=[shift_y, shift_x])

        seqdf = self.peak_detection(blue_thresh, red_thresh, blue_peak_thresh, red_peak_thresh)
        predictedseq = seqdf.base.str.cat()

        fn = self.data_file_name[:-3] + '_seq.fasta'
        predictedseq = Seq.Seq(predictedseq, generic_dna)
        predictedseq = SeqRecord(predictedseq, id=os.path.split(fn)[1])
        SeqIO.write(predictedseq, fn, "fasta")

    def check_controls(self):
        dntp_names = ["dCTP", "dATP", "dGTP", "dTTP"]
        dntps = [[], [], [], []]
        zpro = [[], [], [], []]
        czpro = [[], [], [], []]

        for idx, name in enumerate(dntp_names):
            fn = [filename for filename in os.listdir(self.direc) if filename.startswith(name)
                  and filename.endswith('.h5')]
            hfile = h5py.File(self.direc + os.sep + fn[-1])
            #            print self.direc + os.sep+ fn[-1]
            dntps[idx] = np.array(hfile["images"]).astype(float)
            if self.roi_image_plot:
                dntps[idx] = dntps[idx][:, self.roi_image_plot[1]:self.roip2[1],
                             self.roi_image_plot[0]:self.roip2[0]]
            #            dntps[i] -= self.background

            #            if self.filtertype != []:

            dntps[idx] = dntps[idx][-1500:]
            #            dntps[i] = (dntps[i]/dntps[i].std((1,2))[0])
            #            dntps[i] -= dntps[i].mean()
            #            dntps[i] += self.background.mean()
            zpro[idx] = np.median(dntps[idx], 0)
        #            zpro[i] = np.ma.masked_where(zpro[i] < 0, zpro[i])

        #        zpro[0] -= np.mean(zpro[0])
        #        zpro[1] -= np.mean(zpro[1])
        #        zpro[2] -= np.mean(zpro[2])
        #        zpro[3] -= np.mean(zpro[3])

        fig = plt.figure()
        plt.subplot(1, 3, 1)

        def create_color_plot(data, color):
            color_map = plt.cm.get_cmap(color)
            my_color_map = color_map(np.arange(color_map.N))
            my_color_map[:, -1] = np.linspace(0, 1, color_map.N)
            my_color_map = ListedColormap(my_color_map)
            image = plt.imshow(data, cmap=my_color_map)
            return image

        im1 = create_color_plot(zpro[0], 'Reds')
        im2 = create_color_plot(zpro[1], 'Greens')
        im3 = create_color_plot(zpro[2], 'Blues')
        im4 = create_color_plot(zpro[3], 'PuRd')

        fig.show()

        plt.subplot(1, 3, 2)
        plt.imshow(self.zpro)

        composite = zpro[0] + zpro[1] + zpro[2] + zpro[3]
        ly, lx = composite.shape

        shiftx, shifty = calculate_xy_shift(lx, ly, composite, zpro)

        for idx, name in enumerate(dntps):
            czpro[idx] = np.median(dntps[idx], 0)
            czpro[idx] = ird.transform_img(czpro[idx], tvec=[shifty, shiftx])
            czpro[idx] -= self.background
            czpro[idx] = czpro[idx] / czpro[idx].std()
            czpro[idx] -= czpro[idx].mean()

        composite = czpro[0] + czpro[1] + czpro[2] + czpro[3]
        plt.subplot(1, 3, 3)
        plt.imshow(composite)

    def peak_detection(self, blueThresh, redThresh, bluePeakThresh, redPeakThresh):
        markers = ['C', 'A', 'G', 'T']
        colors = ['r', 'b']
        score = ascore = gscore = tscore = spseries = pd.Series()
        df = pd.DataFrame({
            "ident": [],
            'stimes': [],
            'etimes': [],
            'peaks': [],
            'mins': []
        })
        seqdf = pd.DataFrame({'base': [], 'times': []})
        time1, intensity1 = self.image_plot.roiCurve1.getData()
        time2, intensity2 = self.image_plot.roiCurve2.getData()
        intensities = [intensity1, intensity2]
        self.firingplotlist = []
        peak_series = pd.Series()
        min_series = pd.Series()
        for idx, intensity in enumerate(intensities):
            if idx == 0:
                thresh = redThresh
                noise = redPeakThresh
            elif idx == 1:
                thresh = blueThresh
                noise = bluePeakThresh
            else:
                raise IndexError("Intensity index out of range")

            firing = np.where(((intensity > thresh) & (intensity > intensities[idx - 1])))[0]
            start_and_end = np.diff(firing)
            # Why is two added to the end of each?
            start_points = np.insert(start_and_end, 0, 2)
            end_points = np.insert(start_and_end, -1, 2)

            start_points = np.where(start_points > 1)[0]
            end_points = np.where(end_points > 1)[0]

            start_points = firing[start_points]
            end_points = firing[end_points]
            df = df.append(pd.DataFrame({
                "ident": [idx] * len(start_points),
                'stimes': start_points,
                'etimes': end_points}), ignore_index=True)

            for idx, x in enumerate(start_points):
                sp = start_points[idx]
                ep = end_points[idx] + 1
                curve = pg.PlotDataItem(x=np.linspace(sp, ep, ep - sp),
                                        y=intensity[sp:ep],
                                        pen=pg.mkPen(colors[idx], width=2))
                self.firingplotlist.append(curve)

                self.image_plot.getRoiPlot().addItem(curve)

                try:
                    peaks, mins = peakdet(v=intensity[sp:ep], delta=noise, x=np.arange(sp, ep))
                    if len(peaks) == 0 or len(mins) == 0:
                        peaks = np.NAN
                        substack = np.mean(self.stack[sp:ep], 0)
                        call = get_call(substack, self.czpro, idx)
                        seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [sp]}), ignore_index=True)
                    else:
                        # point = pg.PlotDataItem(peaks, pen = None, symbol = 'o', symbolBrush = 'g')
                        # self.p1.getRoiPlot().addItem(point)
                        # self.firingplotlist.append(point)
                        # point = pg.PlotDataItem(mins, pen = None, symbol = 'o', symbolBrush = 'r')
                        # self.p1.getRoiPlot().addItem(point)
                        # self.firingplotlist.append(point)
                        for idx, x in enumerate(peaks):
                            if idx == 0:
                                ssp = sp
                                sep = int(mins[idx][0])
                            elif idx == len(peaks) - 1:
                                ssp = int(mins[idx - 1][0])
                                sep = ep
                            else:
                                ssp = int(mins[idx - 1][0])
                                sep = int(mins[idx][0])
                            substack = np.mean(self.stack[ssp:sep + 1], 0)
                            call = get_call(substack, self.czpro, idx)
                            seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [ssp]}), ignore_index=True)
                    peak_series = peak_series.append(pd.Series([peaks]))
                    min_series = min_series.append(pd.Series([mins]))
                except Exception as e:
                    raise ValueError

        seqdf = seqdf.sort(['times', 'base'])
        base_colors = {
            'C': 'r',
            'A': 'y',
            'G': 'g',
            'T': 'b'
        }
        for idx, x in enumerate(seqdf.index):
            base = seqdf.base[idx]
            color = base_colors.get(base)
            if base in ['C', 'A']:
                intensity = intensities[0][int(seqdf.times[idx])]
            else:
                intensity = intensities[1][int(seqdf.times[idx])]
            text = pg.TextItem(base, color=color)
            seqplot = self.image_plot.getRoiPlot().addItem(text)
            if idx == 0:
                self.seqplotlist = [text]
            else:
                self.seqplotlist.append(text)
            text.setPos(seqdf.times[idx], intensity)

        print(seqdf.base.str.cat())
        return seqdf


class FilterPopup(QtWidgets.QWidget):
    acceptsignal = QtCore.pyqtSignal(int, int, int, int)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
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
        self.acceptsignal.emit(zfilt, xfilt, yfilt, filtertype)

    def close(self):
        self.hide()


class AnalysisPopup(QtWidgets.QWidget):
    acceptsignal = QtCore.pyqtSignal(float, float, float, float)

    def __init__(self):
        QtWidgets.QWidget.__init__(self)
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
        self.acceptsignal.emit(blueThresh, redThresh,
                               bluePeakThresh, redPeakThresh)

    def close(self):
        self.hide()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = GUIForm()
    myapp.show()
    sys.exit(app.exec_())
