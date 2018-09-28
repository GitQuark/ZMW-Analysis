import os
import h5py
import numpy as np
import pandas as pd
import pyqtgraph as pg
import imreg_dft as ird

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from scipy.ndimage import filters
from scipy.optimize import curve_fit
from peakdetect import peakdet

from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
# Notes:
# dntp - deoxyneucleotide triphosphate, ATGC only
# dATP, dCTP, etc - Video from each CCD, one for each DNA base
# +---+
# | T | Green Illuminated
# | G | D dimmer than T
# | A | Red illuminated
# | C | C dimmer than A
# +---+
# dNTP_#### is image capture from one of four CCD
# dntps  - plural of dntp; all of the bases in file
# dntpss - standard deviation over all time - 2D file; no time
# A to C difference 2-3px; same for T and G
# ChildProcessError 40 microns/side  with ZMW array number depending on fabrication
#
# Z project -
# Base call - Detect the bases


def subtract_background(stack):
    avg = stack.mean((1, 2))
    background = stack[avg < avg.mean() - .5 * avg.std()]
    background = np.median(background, 0)
    stack -= background
    return stack, background


def process_stack(stack):
    stack, background = subtract_background(stack)
    stack -= np.mean(stack)
    stack = stack / stack.std()
    zpro = np.std(stack, 0)
    return stack, background, zpro


def crop(original_stack, roi):
    upper_left = roi.pos().astype(int)
    roi_size = np.array(roi.size()).astype(int)
    lower_right = upper_left + roi_size
    x0, y0 = upper_left
    x1, y1 = lower_right
    stack = original_stack[:, y0:y1, x0:x1]
    roi.setPos(0, 0)
    roi.setSize((int(roi_size[0]), int(roi_size[1])))
    stack, _, zpro = process_stack(stack)
    return stack, zpro


def reset_roi_plot(roi_plot, seq_plot_list, firing_plot_list):
    if seq_plot_list != [[]]:
        for x in seq_plot_list:
            roi_plot.removeItem(x)
    if firing_plot_list:
        for x in firing_plot_list:
            roi_plot.removeItem(x)
    roi_plot.autoRange()


def check_controls(directory, roi_image_plot, roip2, background):
    dntp_names = ["dCTP", "dATP", "dGTP", "dTTP"]
    dntps = [[], [], [], []]
    zpro = [[], [], [], []]
    czpro = [[], [], [], []]

    for idx, name in enumerate(dntp_names):
        file_name = [filename for filename in os.listdir(directory) if filename.startswith(name)
                     and filename.endswith('.h5')]
        file_location = os.path.join(directory, file_name[-1])
        hfile = h5py.File(file_location)
        #            print self.direc + os.sep+ fn[-1]
        dntps[idx] = np.array(hfile["images"]).astype(float)
        if roi_image_plot:
            dntps[idx] = dntps[idx][:, roi_image_plot[1]:roip2[1], roi_image_plot[0]:roip2[0]]
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

    fig = plt.figure(zpro)
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
    plt.imshow(zpro)

    composite = zpro[0] + zpro[1] + zpro[2] + zpro[3]
    ly, lx = composite.shape

    shiftx, shifty = calculate_xy_shift(lx, ly, composite, zpro)

    for idx, name in enumerate(dntps):
        czpro[idx] = np.median(dntps[idx], 0)
        czpro[idx] = ird.transform_img(czpro[idx], tvec=[shifty, shiftx])
        czpro[idx] -= background
        czpro[idx] = czpro[idx] / czpro[idx].std()
        czpro[idx] -= czpro[idx].mean()

    composite = czpro[0] + czpro[1] + czpro[2] + czpro[3]
    plt.subplot(1, 3, 3)
    plt.imshow(composite)


def peak_detection(blueThresh, redThresh, bluePeakThresh, redPeakThresh, image_plot, stack, czpro):
    markers = ['C', 'A', 'G', 'T']
    colors = ['r', 'b']
    score = ascore = gscore = tscore = spseries = pd.Series()
    df = pd.DataFrame({
        "ident": [],
        'stimes': [],
        'etimes': [],
        'maxes': [],
        'mins': []
    })
    seqdf = pd.DataFrame({'base': [], 'times': []})
    time1, intensity1 = image_plot.roiCurve1.getData()
    time2, intensity2 = image_plot.roiCurve2.getData()
    intensities = [intensity1, intensity2]
    # self.firingplotlist = []
    peak_series = pd.Series()
    min_series = pd.Series()
    firing_plot_list = []

    def threshold_data(data, threshold, noise, other_data):
        intensity_high_enough = data > thresh
        intensity_greater_than_other = data > other_data[idx - 1]
        firing_points = np.where(intensity_high_enough & intensity_greater_than_other)[0]
        boundary_points = np.diff(firing_points)
        boundary_idxs = list(np.concatenate(np.argwhere(boundary_points)))
        firing_bounds = list(zip(boundary_idxs[::2], boundary_idxs[1::2]))
        for start, end in firing_bounds:
            # sp = start_points[idx]
            # ep = end_points[idx] + 1
            # end += 1
            x_values = np.arange(start, end)
            intensities = intensity[start: end]
            curve_pen = pg.mkPen(colors[idx], width=2)
            curve = pg.PlotDataItem(x=x_values, y=intensities, pen=curve_pen)
            firing_plot_list.append(curve)
            # Plotting the curves highlighting each base
            image_plot.getRoiPlot().addItem(curve)

            try:
                peaks, mins = peakdet(v=intensities, delta=noise, x=x_values)
                if len(peaks) == 0 or len(mins) == 0:
                    peaks = np.NAN
                    substack = np.mean(stack[start: end], 0)
                    call = get_call(substack, czpro, idx)
                    seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [start]}), ignore_index=True)
                else:
                    # point = pg.PlotDataItem(maxes, pen = None, symbol = 'o', symbolBrush = 'g')
                    # self.p1.getRoiPlot().addItem(point)
                    # self.firingplotlist.append(point)
                    # point = pg.PlotDataItem(mins, pen = None, symbol = 'o', symbolBrush = 'r')
                    # self.p1.getRoiPlot().addItem(point)
                    # self.firingplotlist.append(point)
                    for idx, x in enumerate(peaks):
                        if idx == 0:
                            ssp = start
                            sep = int(mins[idx][0])
                        elif idx == len(peaks) - 1:
                            ssp = int(mins[idx - 1][0])
                            sep = end
                        else:
                            ssp = int(mins[idx - 1][0])
                            sep = int(mins[idx][0])
                        substack = np.mean(stack[ssp:sep + 1], 0)
                        call = get_call(substack, czpro, idx)
                        seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [ssp]}), ignore_index=True)
                peak_series = peak_series.append(pd.Series([peaks]))
                min_series = min_series.append(pd.Series([mins]))
            except Exception as e:
                raise ValueError

    for idx, intensity in enumerate(intensities):
        if idx == 0:
            thresh = redThresh
            noise = redPeakThresh
        elif idx == 1:
            thresh = blueThresh
            noise = bluePeakThresh
        else:
            raise IndexError("Intensity index out of range")
        intensity_high_enough = intensity > thresh
        intensity_greater_than_other = intensity > intensities[idx - 1]
        firing_points = np.where(intensity_high_enough & intensity_greater_than_other)[0]
        boundary_points = np.diff(firing_points)
        boundary_idxs = np.concatenate([[0], np.concatenate(np.argwhere(boundary_points > 1)), [len(boundary_points)]])
        # Skipping doesn't work for intervals of length 1
        firing_bounds = list(zip(firing_points[boundary_idxs[::2]], firing_points[boundary_idxs[1::2]]))

        df = df.append(pd.DataFrame({
            "ident": [idx] * len(firing_bounds),
            'stimes': [x for x, _ in firing_bounds],
            'etimes': [y for _, y in firing_bounds]}), ignore_index=True)

        # Check that length of start and end points is the same
        for start, end in firing_bounds:
            # sp = start_points[idx]
            # ep = end_points[idx] + 1
            # end += 1
            x_values = np.arange(start, end)
            intensities = intensity[start: end]
            curve_pen = pg.mkPen(colors[idx], width=2)
            curve = pg.PlotDataItem(x=x_values, y=intensities, pen=curve_pen)
            firing_plot_list.append(curve)
            # Plotting the curves highlighting each base
            image_plot.getRoiPlot().addItem(curve)

            try:
                maxes, mins = peakdet(v=intensities, delta=noise, x=x_values)
                if len(maxes) == 0 or len(mins) == 0:
                    maxes = np.NAN
                    substack = np.mean(stack[start: end], 0)
                    call = get_call(substack, czpro, idx)
                    seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [start]}), ignore_index=True)
                else:
                    # point = pg.PlotDataItem(maxes, pen = None, symbol = 'o', symbolBrush = 'g')
                    # self.p1.getRoiPlot().addItem(point)
                    # self.firing_plot_list.append(point)
                    # point = pg.PlotDataItem(mins, pen = None, symbol = 'o', symbolBrush = 'r')
                    # self.p1.getRoiPlot().addItem(point)
                    # self.firing_plot_list.append(point)
                    for idx, x in enumerate(maxes):
                        if idx == 0:
                            ssp = start
                            sep = int(mins[idx][0])
                        elif idx == len(maxes) - 1:
                            ssp = int(mins[idx - 1][0])
                            sep = end
                        else:
                            ssp = int(mins[idx - 1][0])
                            sep = int(mins[idx][0])
                        substack = np.mean(stack[ssp:sep + 1], 0)
                        call = get_call(substack, czpro, idx)
                        seqdf = seqdf.append(pd.DataFrame({'base': [call], 'times': [ssp]}), ignore_index=True)
                peak_series = peak_series.append(pd.Series([maxes]))
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
    seq_plot_list = []
    for idx, x in enumerate(seqdf.index):
        base = seqdf.base[idx]
        color = base_colors.get(base)
        if base in ['C', 'A']:
            intensity = intensities[0][int(seqdf.times[idx])]
        else:
            intensity = intensities[1][int(seqdf.times[idx])]
        text = pg.TextItem(base, color=color)
        seqplot = image_plot.getRoiPlot().addItem(text)
        if idx == 0:
            seq_plot_list = [text]
        else:
            seq_plot_list.append(text)
        text.setPos(seqdf.times[idx], intensity)

    print(seqdf.base.str.cat())
    return seqdf, firing_plot_list, seq_plot_list


def analyze(file, czpro, data_file_name):

    seq_plot_list = [[]]
    firing_plot_list = []

    dntp_names = ["dCTP", "dATP", "dGTP", "dTTP"]
    cdf = pd.DataFrame({'a': [], 't': [], 'g': [], 'c': []})
    dntps = [[], [], [], []]
    zpro = [[], [], [], []]

    # TODO: No reason to open another file; use data already loaded
    for idx, name in enumerate(dntp_names):
        zpro[idx] = np.array(file[name]).astype(float)
        zpro[idx] -= zpro[idx].mean()
        zpro[idx] /= zpro[idx].std()

    composite = zpro[0] + zpro[1] + zpro[2] + zpro[3]
    if len(composite.shape) == 2:
        # Single image
        ly, lx = composite.shape
    elif len(composite.shape) == 3:  # Should only be 3
        # Series of frames with time index being first
        lz, ly, lx = composite.shape
    else:
        raise ValueError("Composite dimensions are not 2D or 3D")

    shift_x, shift_y = calculate_xy_shift(lx, ly, composite, zpro)
    czpro = []  # Color zpro?
    for idx, name in enumerate(dntps):
        czpro.append(zpro[idx])
        czpro[idx] = ird.transform_img(zpro[idx], tvec=[shift_y, shift_x])

    # ------------
    seqdf, _, _ = peak_detection(blue_thresh, red_thresh, blue_peak_thresh, red_peak_thresh)
    predicted_seq = seqdf.base.str.cat()

    fn = data_file_name[:-3] + '_seq.fasta'
    predicted_seq = Seq.Seq(predicted_seq, generic_dna)
    predicted_seq = SeqRecord(predicted_seq, id=os.path.split(fn)[1])
    # TODO save as DATA_FILE_NAME + _analyzed.txt
    analyzed_name = data_file_name[:-3] + '_analyzed.txt'
    with open(analyzed_name, 'w') as data_file:
        data_file.write(str(predicted_seq))
    # SeqIO.write(predicted_seq, fn, "fasta")

    return zpro, czpro


def load():
    pass


def gauss(x, *p):
    A, mu1, sigma = p
    return A * np.exp(-(x - mu1) ** 2 / (2. * sigma ** 2))


def dubgauss(x, *p):
    A, mu1, mu2, sigma = p
    return A * np.exp(-(x - mu1) ** 2 / (2. * sigma ** 2)) + A * np.exp(-(x - mu2) ** 2 / (2. * sigma ** 2))


def z_project(image_plot, zpro):
    # Show all frames in one go; currently through stdev of each pixel across time
    # Transpose is here since numpy order is (y, x), qt imageView is (time, x, y)
    image_plot.setImage(np.transpose(256 * (1 + zpro), (1, 0)))


def view_stack(stack, image_plot):
    image_plot.setImage(np.transpose((1 + stack), (0, 2, 1)))
    image_plot.play(rate=60)


def filter_stack(stack, xyz_filter, filter_type, image_plot):
    x_filter, y_filter, z_filter = xyz_filter
    if filter_type == 0:
        stack = filters.median_filter(stack, size=(x_filter, y_filter, z_filter))
    elif filter_type == 1:
        stack = filters.gaussian_filter(stack, sigma=(x_filter, y_filter, z_filter))
    else:
        print('Incorrect filter type: ', filter_type)
    zpro = np.std(stack, 0)
    view_stack(stack, image_plot)
    return stack, zpro


def calculate_xy_shift(lx, ly, composite, zpro):
    p0 = [np.max(np.mean(composite, 1)), (ly / 2) - 4, (ly / 2) + 4, 1]
    coeff1, var_matrix1 = curve_fit(dubgauss, np.linspace(0, ly - 1, ly), np.mean(composite, 1), p0=p0)
    p0 = [np.max(np.mean(zpro, 1)), (ly / 2) - 4, (ly / 2) + 4, 1]
    coeff2, var_matrix2 = curve_fit(dubgauss, np.linspace(0, ly - 1, ly), np.mean(zpro, 1), p0=p0)
    shifty = np.mean((coeff2[1], coeff2[2])) - np.mean((coeff1[1], coeff1[2]))

    p0 = [np.max(np.mean(composite, 1)), lx / 2, 1]
    coeff1, var_matrix1 = curve_fit(gauss, np.linspace(0, lx - 1, lx), np.mean(composite, 0), p0=p0)
    p0 = [np.max(np.mean(zpro, 1)), lx / 2, 1]
    coeff2, var_matrix2 = curve_fit(gauss, np.linspace(0, lx - 1, lx), np.mean(zpro, 0), p0=p0)
    shiftx = coeff2[1] - coeff1[1]
    return shiftx, shifty


def get_call(substack, czpro, idx):
    if idx == 0:
        c_score = np.dot(substack.ravel(), czpro[0].ravel())
        a_score = np.dot(substack.ravel(), czpro[1].ravel())
        if c_score > a_score:
            call = 'C'
        else:
            call = 'A'
    elif idx == 1:
        g_score = np.dot(substack.ravel(), czpro[2].ravel())
        t_score = np.dot(substack.ravel(), czpro[3].ravel())
        if g_score > t_score:
            call = 'G'
        else:
            call = 'T'
    else:
        return
    return call
