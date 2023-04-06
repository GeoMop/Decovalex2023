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

from pathlib import Path
import numpy as np
import mapdfn
from h5py import File






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





def add_field(grid, f, name, data_array):
    # 3d uniform grid
    h5grp = f.create_group(name)
    # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
    h5grp.attrs['Dimension'] = np.string_('XYZ')
    # based on Dimension, specify the uniform grid spacing
    h5grp.attrs['Discretization'] = grid.step
    # again, depends on Dimension
    h5grp.attrs['Origin'] = grid.origin
    # leave this line out if not cell centered.  If set to False, it will still
    # be true (issue with HDF5 and Fortran)
    h5grp.attrs['Cell Centered'] = [True]
    h5grp.attrs['Interpolation Method'] = np.string_('Step')
    h5grp.create_dataset('Data', data=data_array)


class DFN:
    def __init__(self, workdir):
        self.config()
        self.create_dfn(workdir)

    def config(self):
        # Edit these values
        self.grid = mapdfn.Grid(
            origin=np.array([-500, -500, -500]),
            step=np.array([50, 50, 50]),
            dimensions=np.array([20, 20, 20], dtype=int)
        )
        self.k_background = 1.e-24
        self.bulk_por = 0.005
        self.tortuosity_factor = 0.001

        # origin of mapping in cpm domain (0,0,0 is origin of CPM, though it doesn't have to be)
        self.dfn_origin = [-500, -500, -500]
        # origin of area to map in DFN domain coordinates (0,0,0 is center of DFN)

    def create_dfn(self, workdir):
        # Call mapdfn functions
        print('Mapping DFN to grid')
        self.ellipses = mapdfn.readEllipse(workdir)
        self.fractures = mapdfn.map_dfn(self.grid, self.ellipses)
        print('Calculating effective k')
        transmissivity, appertre = mapdfn.fr_transmissivity_apperture(workdir)
        self.porosity = mapdfn.porosity(self.grid, self.fractures, appertre, self.bulk_por).reshape(self.grid.dimensions)
        self.k_iso = mapdfn.permIso(self.grid, self.fractures, transmissivity, self.k_background).reshape(self.grid.dimensions)
        # k_aniso = mapdfn.permAniso(fracture, ellipses, transmissivity, self.grid_step, self.k_background)
        print('Calculating fracture permeability')

    def main_output(self):

        self.xyz = [o + np.linspace(0, 1, n) for o, n in zip(self.grid.origin, self.grid.dimensions)]

        # first fracture index per cell
        self.fracture_idx = np.full(self.grid.dimensions.prod(), -1, dtype=float)
        for i, fr in enumerate(self.fractures):
            self.fracture_idx[fr.cells] = i

        # Write same information to mapELLIPSES.h5. This file can be opened in Paraview
        # by chosing "PFLOTRAN file" as the format.
        with field_file('mapELLIPSES.h5') as ff:
            ff.create_dataset('Coordinates/X [m]', data=self.xyz[0])
            ff.create_dataset('Coordinates/Y [m]', data=self.xyz[1])
            ff.create_dataset('Coordinates/Z [m]', data=self.xyz[2])
            ff.create_dataset('Time:  0.00000E+00 y/Perm', data=self.k_iso)
            ff.create_dataset('Time:  0.00000E+00 y/Fracture', data=self.fracture_idx)
            #ff.create_dataset('Time:  0.00000E+00 y/PermX', data=kx)
            #ff.create_dataset('Time:  0.00000E+00 y/PermY', data=ky)
            #ff.create_dataset('Time:  0.00000E+00 y/PermZ', data=kz)
            ff.create_dataset('Time:  0.00000E+00 y/Porosity', data=self.porosity)

        with field_file('isotropic_k.h5') as f:
            add_field(self.grid, f, 'Permeability', self.k_iso)

        with field_file('porosity.h5') as f:
            add_field(self.grid, f, 'Porosity', self.porosity)

        with field_file('tortuosity.h5') as f:
            add_field(self.grid, f, 'Tortuosity', self.tortuosity_factor / self.porosity)

        # with field_file('anisotropic_k.h5') as f:
        #     add_field(f, 'PermeabilityX', kx)
        #     add_field(f, 'PermeabilityX', ky)
        #     add_field(f, 'PermeabilityX', kz)

        iarray = np.arange(1, self.grid.dimensions.prod() + 1, dtype=int)
        marray = np.zeros(self.grid.dimensions.prod(), dtype=int)
        marray[self.k_iso.flatten() == self.k_background] = 1
        with field_file('materials.h5') as f:
            group = f.create_group("Materials")
            group.create_dataset('Cell Ids', data=iarray)
            group.create_dataset('Material Ids', data=marray)




def main(workdir):
    dfn = DFN(workdir)
    dfn.main_output()


if __name__ == "__main__":
    # get working directory
    if len(sys.argv) > 1:
        workdir = Path(sys.argv[1])
    else:
        workdir = Path.cwd()
    main(workdir)