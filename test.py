import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pyqtgraph as pg
import numpy as np

imdat1 = np.random.random_integers(2,size = (100,100))
imdat2 = np.random.random_integers(2,size = (100,100))

#fig = plt.figure()
#
#cmap = plt.cm.Reds
#my_cmap = cmap(np.arange(cmap.N))
#my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
#my_cmap = ListedColormap(my_cmap)
#im1 = plt.imshow(imdat1, cmap = my_cmap)
#
#cmap = plt.cm.Blues
#my_cmap = cmap(np.arange(cmap.N))
#my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
#my_cmap = ListedColormap(my_cmap)
#im2 = plt.imshow(imdat2, cmap = my_cmap)
p1 = pg.plot()

#pos = np.array([0.0, 0.5, 1.0])
#color = np.array([[0,0,0,255], [0,0,0,255], [0,255,0,255]], dtype=np.ubyte)
#cmap = pg.ColorMap(pos, color)
#lut = cmap.getLookupTable(0.0, 1.0, 256)

pg.GradientWidget()

it1 = pg.ImageItem(imdat1)
it1.setLookupTable(lut)
p1.addItem(it1)