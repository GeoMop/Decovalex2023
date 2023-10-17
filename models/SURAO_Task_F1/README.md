# SURAO Model for Task F1 of the Decovalex 2023 Project

Ondrej Mikláš
Jan Březina

## Model description

Model is a proof of concept for a homogenization approach to the repository scale transport model.
A set of fractures generated as a sample of a DFN model is homogenized by a Python script to a structured grid. 
The homogenization is done for effective isotropic permeability and porosity. 
The Python script is also used for preparation of the spatial fields describing the homogenized contamination sources.
The the convection transport model of the conservative tracer is performed in PFlotran. The diffusion through the bentonite
buffer is described by a dual field model. On the cells containing the *Distributed Waste Package*, we consider 
mobile and immobile concentration of the tracer, where the immobile concentration represent the concentration in the container. 
The linear diffusion coeeficeint is determined from the geometrical and diffusion parameters of the container buffer. 
The immobile concentration is same for all containers, given as a time series. This time series is calculated by a simple 
model of initial fast degradation followid by the slow longterm degradation. 
See the project report for the details


## Model evaluation process
Whole model could be executed by the `run.sh` Bash script.
The main steps are:
1. Input files `input_dfn/*.dat`. 
2. Execution of `../../decodfn/main.py` (the symbolic link `main.py` is available for convenience)
   with `main.yaml` input file. 
3. *.h5 files with heterogeneous fields are produced
4. pflotran.in and lambda.dat PFlotran input files
4. Execution of the PFlotran model using the `flow123d/pflotran-gnu-rel` Docker image with our build of PFlotran 
   with custom reaction sandbox DWP model. The image is pulled from the public docker hub cloud repository.

   
with the fracture set
The generated sample of a DFN model is given by the `input_dfn` folder containing parameters of the fracture set.
The fracture set is then homogenized into a structured grid by the `../../src/main.py` script with the `main.yaml` input file.

## Repository Files (scripts and inputs)
lambda.dat - PFlotran chemistry database 
regions_2.txt - external file with region definitions, included from pflotran.in
pressure_boundary.h5 - prescribed field on the boundary (resource ??)
pflotran.in - PFlotran man input
stoch_fractures_WP.in - Some other ?? PFlotran input
main.py - symbolic link to the main script for homogenization and preparation of the PFlotran inputs 


## JB dev notes

- napodařilo se spustit, potřebuje speciální lambda.dat s Tracer_1
- asi rozumnější použít Iod jak tracer, bez nutnosti do tabulky zasahovat
- zahrnutí výpočtu zdrojového členu v update_vypoctu_zdroje;
  nutno mít test a zahrnout do repozitáře s pflotranem
  
  vypadá to, že výpočet waste_rate kopíruje dacay_rate, ale je tam asi chybně výpočet Jacobiánu a
  navíc se pomocí UpdateKinetics dále zmenšuje objem zdrojového členu, no úplně jasné mi to není, 
  protože to dělá vlastně ten dacay
