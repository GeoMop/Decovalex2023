'''
Usage:
    python mapdfn2pflotran.py <workdir>

Uses mapdfn.py to take output of dfnWorks-Version2.0, create
equivalent continuous porous medium representation, and write parameters
(permeability, porosity, tortuosity) to files for use with PFLOTRAN.

   Usage: Edit values for origin, nx, ny, nz, d, k_background, bulk_por,
          tortuosity factor, and h5origin.
          Paths and filenames are hardwired and may also need to be checked.
          As written, they assume script is being called from a subdirectory.
          Then: python mapdfn2pflotran.py
'''

# Access to modules in the same directory.
# Wihtout making the whole thing a Python package.
import os
import sys
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)

import numpy as np
from mapdfn import *
from h5py import File



# Edit these values
nx = 20
ny = 20
nz = 20
grid_step = 50
k_background = 1.e-24
bulk_por = 0.005
tortuosity_factor = 0.001
origin = [-500, -500, -500]
# origin of mapping in cpm domain (0,0,0 is origin of CPM, though it doesn't have to be)
dfn_origin = [-500, -500, -500]
# origin of area to map in DFN domain coordinates (0,0,0 is center of DFN)

# Call mapdfn functions
print('Mapping DFN to grid')
ellipses = readEllipse('../radii_Final.dat','../normal_vectors.dat','../translations.dat')
fracture = map_dfn(ellipses, dfn_origin, nx, ny, nz, [grid_step, grid_step, grid_step])
print('Calculating effective k')
T = findT('../aperture.dat','../perm.dat')
k_iso = permIso(fracture, T, grid_step, k_background)
k_aniso = permAniso(fracture, ellipses, T, grid_step, k_background)
print('Calculating fracture permeability')
por = porosity(fracture, grid_step, bulk_por, '../aperture.dat')

# Fill arrays that will go into PFLOTRAN h5 files.
# Also write mapELLIPSES.txt for later use.
# Potentially large file. If you have no use for it, don't write it.
with open('mapELLIPSES.txt','w') as fout:
    fout.write('#x, y, z, number of fractures, fracture numbers\n')

    #arrays for making PFLOTRAN hf file
    x = origin[0] + grid_step * np.linspace(0, 1, nx)
    y = origin[0] + grid_step * np.linspace(0, 1, ny)
    z = origin[0] + grid_step * np.linspace(0, 1, nz)

    n_cells = nx * ny * nz
    a = np.zeros((nx,ny,nz),'=f8')
    khdf5 = np.zeros((nx,ny,nz),'=f8')
    kx = np.zeros((nx,ny,nz),'=f8')
    ky = np.zeros((nx,ny,nz),'=f8')
    kz = np.zeros((nx,ny,nz),'=f8')
    phdf5 = np.zeros((nx,ny,nz),'=f8')

    print ('Writing text file')
    #write mapELLIPSES.txt and fill arrays
    for k in range(nz):
      for j in range(ny):
        for i in range(nx):
          index = i+nx*j+nx*ny*k
          khdf5[i][j][k] = k_iso[index]
          kx[i][j][k] = k_aniso[index][0]
          ky[i][j][k] = k_aniso[index][1]
          kz[i][j][k] = k_aniso[index][2]
          phdf5[i][j][k] = por[index]
          fout.write('%e  %e  %e  %i %e '
                     % (dfn_origin[0] + i * grid_step + grid_step / 2.,
                        dfn_origin[1] + j * grid_step + grid_step / 2.,
                        dfn_origin[2] + k * grid_step + grid_step / 2.,
                        fracture[index][0], k_iso[index]))
          if (fracture[index][0]) != 0:
                a[i][j][k] = fracture[index][1] #color by the first fracture number in the list
                for c in range(1,fracture[index][0]+1):
                    fout.write(' '+str(fracture[index][c]))
          else:
                a[i][j][k] = 0 #color it zero
                fout.write(' '+str(fracture[index][1])) #?
    fout.write('\n')


  def field_file(fname):
    """
    Open a new HDF5 file and return its handle.
    Usage:

    with field_file('abc.h5') as f:
        add_field(...)
    this way file is closed automatically.
    """
    print(f"Creating {fname} data file.")
    return File(fname, 'w')

# Write same information to mapELLIPSES.h5. This file can be opened in Paraview
# by chosing "PFLOTRAN file" as the format.
with field_file('mapELLIPSES.h5') as ff:
  ff.create_dataset('Coordinates/X [m]', data=x)
  ff.create_dataset('Coordinates/Y [m]', data=y)
  ff.create_dataset('Coordinates/Z [m]', data=z)
  ff.create_dataset('Time:  0.00000E+00 y/Perm', data=khdf5)
  ff.create_dataset('Time:  0.00000E+00 y/Fracture', data=a)
  ff.create_dataset('Time:  0.00000E+00 y/PermX', data=kx)
  ff.create_dataset('Time:  0.00000E+00 y/PermY', data=ky)
  ff.create_dataset('Time:  0.00000E+00 y/PermZ', data=kz)
  ff.create_dataset('Time:  0.00000E+00 y/Porosity', data=phdf5)




def add_field(f, name, data_array):
    # 3d uniform grid
    h5grp = f.create_group(name)
    # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
    h5grp.attrs['Dimension'] = np.string_('XYZ')
    # based on Dimension, specify the uniform grid spacing
    h5grp.attrs['Discretization'] = [grid_step, grid_step, grid_step]
    # again, depends on Dimension
    h5grp.attrs['Origin'] = origin
    # leave this line out if not cell centered.  If set to False, it will still
    # be true (issue with HDF5 and Fortran)
    h5grp.attrs['Cell Centered'] = [True]
    h5grp.attrs['Interpolation Method'] = np.string_('Step')
    h5grp.create_dataset('Data', data=data_array)

with field_file('isotropic_k.h5') as f:
  add_field(f, 'Permeability', khdf5)

with field_file('porosity.h5') as f:
  add_field(f, 'Porosity', phdf5)

with field_file('tortuosity.h5') as f:
    add_field(f, 'Tortuosity', tortuosity_factor/phdf5)

with field_file('anisotropic_k.h5') as f:
    add_field(f, 'PermeabilityX', kx)
    add_field(f, 'PermeabilityX', ky)
    add_field(f, 'PermeabilityX', kz)

with field_file('materials.h5') as f:
    group = f.create_group()
    iarray = np.zeros((nx*ny*nz),'=i4')
    marray = np.zeros((nx*ny*nz),'=i4')
    for i in range(nx*ny*nz):
      iarray[i] = i+1
      if k_iso[i] == k_background:
        marray[i] = 0
      else:
        marray[i] = 1
    group.create_dataset('Cell Ids', data=iarray)
    group.create_dataset('Material Ids', data=marray)


