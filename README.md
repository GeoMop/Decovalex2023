# Decovalex2023

Tool for upsaling of a Discrete Fracture Network and repository sources into eqivalent 
heterogeneous continuum on a structed PFlotran grid.

Based on original SANDIA scripts (see below).


## Instalation
First download the [project package](https://github.com/GeoMop/Decovalex2023/archive/refs/heads/master.zip) to the `decovalex` directory. 
For development you may prefer to git clone the repository:
    
    ```
    git clone https://github.com/GeoMop/Decovalex2023.git
    ```

For installation to user space (no root / administrator rights required) change to the `decovalex` directory
and execute:

    ```
    python3 -m pip install ./decovalex
    ```

Posisbly the direct installation from Github may work as well:
    ```
    python3 -m pip install "decovalex @ git+https://github.com/GeoMop/Decovalex2023
    ``` 
    
# Usage
```
    decodfn <workdir>
```

`workdir` is the direcotry with the `main.yaml` configuration file and possible subdirectory `dfn_input` with fracture properties.
The output files are placed in current system directory.


# PFlotran with auxiliary reaction sandbox
On Windows:
1. install [Git for Windows](https://gitforwindows.org/) that provides Git BASH.
2. install [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)

On Linux just [install Docker](https://docs.docker.com/engine/install/ubuntu/) for Debian/Ubuntu bases distributions. 
For RedHat based distributions se e.g. [installation on CentOS](https://docs.docker.com/engine/install/centos/).

Having working terminal, installed Docker, and cloned [decovalex repository](https://github.com/GeoMop/Decovalex2023)
in a local `decovalex` directory one can run e.g. the `dwp_test` simulation as follows:
```
    cd decovalex/pflotran/dwp_test
    ../fterm.sh 
```

## Repository structure
- `decodfn` 
    - `main.py` : the main script, the `main()` called by the command `decodfn`
- `pflotran` : Docker file and test for building a pflotran with DWP diffusion reaction sandbox model
    - `pflotran_JB` : submodule with pflotran fork 
    - `dwp_test` : the test input for the Distributed Waste Package diffusion model
    - `build-image` : Dockerfile for building the pflotran image (use flow123D build images with compatible PETSC version)
    - `makefile` : targets for development and release builds 
    - `fterm.sh` : simple script for starting development Docker container
- `tests` : test calculations 
    - `dfn_253` : test sample of 253 fractures
    - `repo_multiscale` : basic multiscale repository model
    

## Diffusion Waste PAckage model
The modified Pflotran sources are in the [forked repository](https://github.com/flow123d/pflotran.git).
To build the release image use:
```
    cd ./pflotran && make rel_image
    docker push flow123d/pflotran-gnu-rel
```

### How to use the DWP model
We describe steps to execute the DWP test:

1. [Install Docker](https://docs.docker.com/engine/install/)
2. Install Git (versioning system)
    - Linux, Ubuntu: e.g. `sudo apt install git`
    - Windows: use [Git for Windows](https://gitforwindows.org/) which also provides the Git BASH terminal
3. Running the test. In terminal (native on Linux, Git Bash on Windows):
     
     ```
        cd ./pflotran/dwp_test                      # entry the test directory
        docker pull flow123d/pflotran-gnu-rel       # retrive the pflotran image with EOS diffusion model
        docker 
     ```
    
## Original scripts authors

Authors: Emily Stein (ergiamb@sandia.gov)
        Kris Kuhlman (klkuhl@sandia.gov)
        Applied Systems Analysis and Research, 8844
        Sandia National Laboratories

Date: 07/13/18
SAND Number: SAND2018-7604 O
