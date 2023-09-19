
# Access to modules in the same directory.
# Wihtout making the whole thing a Python package.
# import os

import logging
import sys

logging.basicConfig(level=logging.INFO) #, filename='endorse_mlmc.log')
#logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


from pathlib import Path

import attrs
import matplotlib.pyplot as plt
import numpy as np
import yaml
from h5py import File

import mapdfn


# script_dir = os.path.dirname(os.path.realpath(__file__))
# sys.path.append(script_dir)


def field_file(fname):
    """
    Open a new HDF5 file and return its handle.
    Usage:

    with field_file('abc.h5') as f:
        add_field(...)
    this way file is closed automatically.
    """
    #print(f"Creating {fname} data file.")
    return File(fname, 'w')






@attrs.define
class RockMass:
    permeability: float
    porosity: float
    tortuosity_factor: float

@attrs.define
class RepositoryBlock:
    origin: np.array = attrs.field(converter=mapdfn.float_array)
    drift_direction: np.array = attrs.field(converter=mapdfn.int_array)
    drift_lengh: float
    drift_step: np.array = attrs.field(converter=mapdfn.float_array)
    n_drifts: int
    hole_spacing: float

    @staticmethod
    def _axis_unit(x):
        assert np.sum(x) in {1, -1}
        assert np.linalg.norm(x) == 1

    def validate(self):
        self._axis_unit(self.drift_direction)
        self._axis_unit(self.drift_step / np.linalg.norm(self.drift_step))

@attrs.define
class SurfaceBC:
    x0: float
    x1: float
    z0: float
    z1: float
    pressure_at_1000: float
    gravity_acceleration: float = 9.89
    water_density: float = 1000

    def surface_field(self, grid: mapdfn.Grid):
        """
        sin(x) transition from ZO to Z1 on the interval (X0, X1)
        """
        elevation = np.empty(grid.cell_dimensions[0:2], dtype=float)
        ix0 = grid.cell_coord([self.x0, 0,0])[0]
        ix1 = grid.cell_coord([self.x1, 0, 0])[0]
        elevation[0:ix0, :] = self.z0
        transition = 0.5 - 0.5 * np.sin(np.linspace(np.pi/2, np.pi*3/2, ix1 - ix0))
        # transition from 0 to 1 on the range ix0:ix1
        elevation[ix0:ix1, :] = transition[:, None] * (self.z1 - self.z0) + self.z0
        elevation[ix1:, :] = self.z1
        return self.pressure_at_1000 + (elevation - 1000) * self.water_density * self.gravity_acceleration


class DFN:
    @staticmethod
    def read_main_yaml(workdir: Path):
        with open(workdir / "main.yaml") as f:
            cfg = yaml.load(f, Loader=yaml.SafeLoader)
        return cfg

    def __init__(self, workdir: Path):
        self.workdir = workdir
        # input directory
        self.output_fields = {}
        # dictionary of the filds for the Paraview output

        # Read main input file
        cfg = self.read_main_yaml(workdir)
        self.grid = mapdfn.Grid.make_grid(**cfg['grid'])
        # configuration of the strutured grid
        self.rock_mass = RockMass(**cfg['rock'])
        # Bulk material parameters
        self.dfn_origin = cfg['dfn_origin']
        # Origion of the DFN domain, not used yet
        repo_cfg = cfg.get('repository', [])
        self.repository = [RepositoryBlock(**b) for b in repo_cfg]
        # Repository configuration (optional)
        bc_cfg = cfg.get("surface_bc", None)
        self.surface_bc = SurfaceBC(**bc_cfg) if bc_cfg is not None else None
                
        if bool(cfg.get("mean_eq_porosity", False)):
            self.porosity_fn = mapdfn.porosity_mean
        else:
            self.porosity_fn = mapdfn.porosity_min    
        self.ellipses = []
        # List of DFN fractures.
        self.fractures = []
        # Fractures with related cells.


    def create_dfn(self):
        # Call mapdfn functions
        logging.log(logging.INFO, 'Mapping DFN to grid...')
        self.ellipses = mapdfn.readEllipse(self.workdir / "input_dfn")
        logging.log(logging.INFO, f'loaded {len(self.ellipses)} fracture ellipses')
        self.fractures = mapdfn.map_dfn(self.grid, self.ellipses)


    def mark_line(self, mask, x, vec, length):
        """
        Mark cells intersected by line.
        DEal with lines out of the grid.
        """
        ia = self.grid.cell_coord(self.grid.trim(x))
        ia = np.maximum(np.zeros_like(ia), ia)
        ib = self.grid.cell_coord(self.grid.trim(x + length * vec))
        ib = np.minimum(self.grid.cell_dimensions - 1, ib)
        slices = tuple( [ a if s==0 else np.s_[a:b:s] for a,b,s in zip(ia,ib,vec) ])
        mask[slices] = 1

    def repository_fields(self):  # bude pouzito i pro generovani souboru pro Reaction Sandbox

        repo_cells = np.zeros(self.grid.cell_dimensions, dtype=int)

        for block in self.repository:
            block.validate()
            for i_drift in range(block.n_drifts):
                drift_origin = block.origin + i_drift * block.drift_step
                self.mark_line(repo_cells, drift_origin, block.drift_direction, block.drift_lengh)
        return repo_cells




    
        
####################################################################################################################################š
        """ 
        Vcházel jsem z funkce repository_fields, která bere hodnoty
        z yaml souboru, ale parametry virtualniho mineralu a difuzni 												
        konstanty jsem zatim ve zkusebni verzi zadaval naprimo do skriptu											
        """																									 
    def init_instant_fields(self):																					
        repo_cells = np.zeros(self.grid.cell_dimensions, dtype=int)													
        b_conc = 1e-25                                               # SW vyžaduje nenulovou hodnotu na pozadí		
        for block in self.repository:																				
            block.validate()																						
            for i_drift in range(block.n_drifts):																
                drift_origin = block.origin + i_drift * block.drift_step										
                self.mark_line(repo_cells, drift_origin, block.drift_direction, block.drift_lengh)		
        repo_cells = repo_cells * 2.49e-10                            # hodnota v miste UOS						
        repo_cells[repo_cells == 0] = b_conc																		
        return repo_cells 																							
        
    def init_fractional_fields(self):
        repo_cells = np.zeros(self.grid.cell_dimensions, dtype=int)
        b_conc = 1e-25
        for block in self.repository:
            block.validate()
            for i_drift in range(block.n_drifts):
                drift_origin = block.origin + i_drift * block.drift_step
                self.mark_line(repo_cells, drift_origin, block.drift_direction, block.drift_lengh)
        repo_cells = repo_cells * 2.24e-12
        repo_cells[repo_cells == 0] = b_conc
        return repo_cells 
        
    def diffusion_rate_fields(self):
        repo_cells = np.zeros(self.grid.cell_dimensions, dtype=int)
        for block in self.repository:
            block.validate()
            for i_drift in range(block.n_drifts):
                drift_origin = block.origin + i_drift * block.drift_step
                self.mark_line(repo_cells, drift_origin, block.drift_direction, block.drift_lengh)
        return repo_cells * 1e-9
        
    def decay_fields(self):
        repo_cells = np.zeros(self.grid.cell_dimensions, dtype=int)
        for block in self.repository:
            block.validate()
            for i_drift in range(block.n_drifts):
                drift_origin = block.origin + i_drift * block.drift_step
                self.mark_line(repo_cells, drift_origin, block.drift_direction, block.drift_lengh)
        return repo_cells * 1e-7
################################################################################################################# 
        
        



    def surface_plot(self, field_xy):
        fig, ax = plt.subplots()
        im = ax.imshow(field_xy, cmap='viridis', interpolation='nearest')
        cbar = ax.figure.colorbar(im, ax=ax)
        cbar.ax.set_ylabel('Color bar', rotation=-90, va="bottom")
        fig.savefig('pressure_top.pdf')

    def add_field(self, fname, name, data_array, axes=(0,1,2)):
        self.output_fields[name] = data_array

        with field_file(fname) as f:
            # 3d uniform grid
            h5grp = f.create_group(name)
            # 3D will always be XYZ where as 2D can be XY, XZ, etc. and 1D can be X, Y or Z
            h5grp.attrs['Dimension'] = np.string_(''.join(['XYZ'[ax] for ax in axes]))
            # based on Dimension, specify the uniform grid spacing
            h5grp.attrs['Discretization'] = [self.grid.step[ax] for ax in axes]
            # again, depends on Dimension
            h5grp.attrs['Origin'] = [self.grid.origin[ax] for ax in axes]
            # leave this line out if not cell centered.  If set to False, it will still
            # be true (issue with HDF5 and Fortran)
            h5grp.attrs['Cell Centered'] = [True]
            h5grp.attrs['Transient'] = False
            h5grp.attrs['Interpolation Method'] = np.string_('Step')
            h5grp.create_dataset('Data', data=data_array)

    def crate_fields(self):
        transmissivity, appertre = mapdfn.fr_transmissivity_apperture(self.workdir / "input_dfn")        
        porosity = self.porosity_fn(self.grid, self.fractures, appertre, self.rock_mass.porosity)
        k_iso = mapdfn.permIso(self.grid, self.fractures, transmissivity, self.rock_mass.permeability)
        # k_aniso = mapdfn.permAniso(fracture, ellipses, transmissivity, self.grid_step, self.k_background)
        self.add_field('porosity.h5', 'Porosity', porosity)
        self.add_field('isotropic_k.h5', 'Permeability', k_iso)
        self.add_field('tortuosity.h5', 'Tortuosity', self.rock_mass.tortuosity_factor / porosity)

        if len(self.repository) != 0:
            repository_mask = self.repository_fields()
            self.add_field('repository.h5', 'Repository', repository_mask)

        if self.surface_bc:
            bc_top_pressure = self.surface_bc.surface_field(self.grid)
            self.surface_plot(bc_top_pressure)
            self.add_field('pressure_top.h5', 'bc_top_presuure', bc_top_pressure, axes=(0,1))
            
            

########################################################################################################################
        """											  
		V pripade toho typu datasetu (v dokumentaci nazvany Gridded dataset)
		to podle vizualizace v Paraview vypada zatim funkcne
        """
        decay = self.decay_fields()	     
        init_instant = self.init_instant_fields()										
        init_fractional = self.init_fractional_fields()									
        diffusion_rate = self.diffusion_rate_fields()									
        self.add_field('init_instant.h5', 'Init_instant', init_instant)					
        self.add_field('init_fractional.h5', 'Init_fractional', init_fractional)     
        self.add_field('decay.h5', 'DecayRate', decay)	     
        self.add_field('diffusion.h5', 'DiffRate', diffusion_rate)
#############################################################################################################################																				




    def main_output(self):

        self.xyz = [o + s * np.arange(0, n + 1) for o, n, s in zip(self.grid.origin, self.grid.cell_dimensions, self.grid.step)]

        # first fracture index per cell
        self.fracture_idx = np.full(self.grid.cell_dimensions.prod(), -1, dtype=float)
        for i, fr in enumerate(self.fractures):
            self.fracture_idx[fr.cells] = i

        # Write same information to mapELLIPSES.h5. This file can be opened in Paraview
        # by chosing "PFLOTRAN file" as the format.



        
#############################################################################################
        """
		Tento postup asi nespravne (nebo alespon neocekavane) prirazuje souradnice, protoze ve 
		vyslednem rozlozeni UOS je evidentne nejaky system, ale neni konzistantny z HU v ramci
		datasetu vyse
        """

        decay = self.decay_fields()
        init_instant = self.init_instant_fields() 
        init_fractional = self.init_fractional_fields()
        diffusion_rate = self.diffusion_rate_fields()          															
																							
        with File('input_fields.h5', 'w') as ff:											
            ff.create_dataset('Cell Ids', data=np.arange(1, self.grid.cell_dimensions.prod() + 1, dtype=int))													
            ff.create_dataset('Coordinates/X [m]', data=self.xyz[0])							
            ff.create_dataset('Coordinates/Y [m]', data=self.xyz[1])						
            ff.create_dataset('Coordinates/Z [m]', data=self.xyz[2])        				
            ff.create_dataset('DecayRate', data=decay.flatten())							
            ff.create_dataset("InitInstant", data=init_instant.flatten())
            ff.create_dataset("InitFractional", data=init_fractional.flatten())
            ff.create_dataset("DiffusionRate", data=diffusion_rate.flatten())
#############################################################################################   




        
        with field_file('mapELLIPSES.h5') as ff:
            ff.create_dataset('Coordinates/Z [m]', data=self.xyz[2])
            ff.create_dataset('Coordinates/Y [m]', data=self.xyz[1])
            ff.create_dataset('Coordinates/X [m]', data=self.xyz[0])
            for name, field in self.output_fields.items():
                ff.create_dataset(f'Time:  0.00000E+00 y/{name}', data=field)
            #ff.create_dataset('Time:  0.00000E+00 y/Fracture', data=self.fracture_idx)
            #ff.create_dataset('Time:  0.00000E+00 y/PermX', data=kx)
            #ff.create_dataset('Time:  0.00000E+00 y/PermY', data=ky)
            #ff.create_dataset('Time:  0.00000E+00 y/PermZ', data=kz)


        # with field_file('anisotropic_k.h5') as f:
        #     add_field(f, 'PermeabilityX', kx)
        #     add_field(f, 'PermeabilityX', ky)
        #     add_field(f, 'PermeabilityX', kz)

        # iarray = np.arange(1, self.grid.cell_dimensions.prod() + 1, dtype=int)
        # marray = np.zeros(self.grid.cell_dimensions.prod(), dtype=int)
        # marray[self.k_iso.flatten() == self.rock_mass.permeability] = 1
        # with field_file('materials.h5') as f:
        #     group = f.create_group("Materials")
        #     group.create_dataset('Cell Ids', data=iarray)
        #     group.create_dataset('Material Ids', data=marray)


def process_dir(workdir):
    dfn = DFN(workdir)
    dfn.create_dfn()
    dfn.crate_fields()
    dfn.main_output()

def main():
    # get working directory
    if len(sys.argv) > 1:
        workdir = Path(sys.argv[1])
    else:
        workdir = workdir = Path.cwd()
    process_dir(workdir)




if __name__ == "__main__":
    main()
