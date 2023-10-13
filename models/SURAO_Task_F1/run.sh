#!/bin/bash

# needs input_dfn with DFN input files

#python3 main.py 

# homogenization of 10k fractures to  about 600k cells takes about: 4.4m
# 4 - 3000 cell intersection, time per fracture significantly depends on the number of intersections
#  => efficient intersection heuristic
# first intersection 3600 cells tooks about 15s of 4.4m

# output h5 files: 

# decay.h5
# diffusion.h5
# init_fractional.h5
# init_instant.h5
# input_fields.h5
# isotropic_k.h5
# mapELLIPSES.h5
# porosity.h5
# repository.h5
# tortuosity.h5


# Generats hdf5 input files for pflotran
docker pull flow123d/pflotran-gnu-rel
docker run -it  -v $(pwd):$(pwd) -w $(pwd) flow123d/pflotran-gnu-rel 
