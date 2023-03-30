"""
A script to create a top BC condition field for a PFlotran calculation.

Assumes a regular (nx, ny) grid.
The BC field name: `surface_top`
"""
from h5py import *
import matplotlib.pyplot as plt
import numpy as np
import math

#==========1D pole - skusobne===============
nx = 50
ny = 50

x_range = np.linspace(0, 5000, 50)
x_norm = np.linspace(1, 3, 19)      #pocet buniek/3 pre vypocet sklonu
x_new = []
rarray = np.zeros((nx, ny), '=f8')

for i in range(0, 17):
    x = 1020
    x_new.append(x)
for i, y in zip(range(17, 36), x_norm):
    x = 10*math.sin(y*math.pi/2)+1010
    x_new.append(x)
for i in range(36, 50):
    x = 1000
    x_new.append(x)

fig, ax = plt.subplots()
ax.plot(x_range, x_new)
plt.show()

#======2D pole pre hdf5 dataset===========
L = 100   # rozmery bunky
dx = 5000  # rozmery domeny
dy = 2000   
dz = 1000
nx = int(dx/L) # pocet buniek
ny = int(dy/L)
nz = dz/L

p_0 = 101325.  # referencny tlak v 1000 m (Pa)


#x_left = np.arange(0., 1700.+L, L)
#xl = len(x_left)
x_middle = np.arange(1700., 3700.+L, L)
xm = len(x_middle)
#x_right = np.arange(3700., 5000.+L, L)
#xr = len(x_right)

x_norm = np.linspace(1, 3, xm)

rarray = np.zeros((nx, ny), '=f8')

for i in range(0, 17):
    for j in range(ny):
        rarray[i][j] = p_0+9800*(1020 - 1000)

for i, y in zip(range(17, 36), x_norm):
    for j in range(ny):
        x = 10*math.sin(y*math.pi/2)+1010   
        #y = x
        rarray[i][j] = p_0+9800*(x-1000)

for i in range(36, 50):
    for j in range(ny):
        rarray[i][j] = p_0  

a = rarray
plt.imshow(a, cmap='Blues', interpolation='nearest')
heatmap = plt.pcolor(a)
plt.colorbar(heatmap)
plt.show()


# ===============hdf5=============================
filename = 'pressure_boundary.h5'
h5file = File(filename, mode='w')

h5grp = h5file.create_group('surface_top')
h5grp.attrs['Dimension'] = np.string_('XY')
h5grp.attrs['Discretization'] = [1., 1.]
h5grp.attrs['Origin'] = [0., 0.]
h5grp.attrs['Transient'] = False

h5dset = h5grp.create_dataset('Data', data=rarray)
h5file.close()
