#!/bin/bash

# needs input_dfn with DFN input files

python3 main.py

# Generats hdf5 input files for pflotran

docker run -it  -v $(pwd):$(pwd) -w $(pwd) flow123d/pflotran-gnu-rel 
