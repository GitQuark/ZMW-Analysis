import numpy as np
from scipy.ndimage import filters
from scipy.optimize import curve_fit


def subtract_background(stack, background):
    if background:
        avg = stack.mean((1, 2))
        background = stack[avg < avg.mean() - .5 * avg.std()]
        background = np.median(background, 0)
        stack -= background
    return stack


def load():
    pass


def dubgauss(x, *p):
    A, mu1, mu2, sigma = p
    return A * np.exp(-(x - mu1) ** 2 / (2. * sigma ** 2)) + A * np.exp(-(x - mu2) ** 2 / (2. * sigma ** 2))


def gauss(x, *p):
    A, mu1, sigma = p
    return A * np.exp(-(x - mu1) ** 2 / (2. * sigma ** 2))


def z_project(image_plot, zpro):
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
    p0 = [max(np.mean(composite, 1)), (ly / 2) - 4, (ly / 2) + 4, 1]
    coeff1, var_matrix1 = curve_fit(dubgauss, np.linspace(0, ly - 1, ly), np.mean(composite, 1), p0=p0)
    p0 = [max(np.mean(zpro, 1)), (ly / 2) - 4, (ly / 2) + 4, 1]
    coeff2, var_matrix2 = curve_fit(dubgauss, np.linspace(0, ly - 1, ly), np.mean(zpro, 1), p0=p0)
    shifty = np.mean((coeff2[1], coeff2[2])) - np.mean((coeff1[1], coeff1[2]))

    p0 = [max(np.mean(composite, 1)), lx / 2, 1]
    coeff1, var_matrix1 = curve_fit(gauss, np.linspace(0, lx - 1, lx), np.mean(composite, 0), p0=p0)
    p0 = [max(np.mean(zpro, 1)), lx / 2, 1]
    coeff2, var_matrix2 = curve_fit(gauss, np.linspace(0, lx - 1, lx), np.mean(zpro, 0), p0=p0)
    shiftx = coeff2[1] - coeff1[1]
    return shiftx, shifty


def get_call(substack, czpro, idx):
    if idx == 0:
        cscore = np.dot(substack.ravel(), czpro[0].ravel())
        ascore = np.dot(substack.ravel(), czpro[1].ravel())
        if cscore > ascore:
            call = 'C'
        else:
            call = 'A'
    elif idx == 1:
        gscore = np.dot(substack.ravel(), czpro[2].ravel())
        tscore = np.dot(substack.ravel(), czpro[3].ravel())
        if gscore > tscore:
            call = 'G'
        else:
            call = 'T'
    else:
        return
    return call
