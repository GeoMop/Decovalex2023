#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
set -e    # exit on error

# needs input_dfn with DFN input files
src_dir="${SCRIPT_DIR}/../../decodfn"
#python3 "${src_dir}/main.py"

# homogenization of 10k fractures to  about 600k cells takes about: 4.4m
# 4 - 3000 cell intersection, time per fracture significantly depends on the number of intersections
#  => efficient intersection heuristic
# first intersection 3600 cells tooks about 15s of 4.4m

# output: input_fields.h5


#exit

# Generats hdf5 input files for pflotran
image=flow123d/pflotran-gnu-rel:OM_decay_source
docker pull ${image}

pf="docker run -it  -v $(pwd):$(pwd) -w $(pwd) ${image}"

#${pf} -pflotranin stoch_fractures_WP_flow.in
# 12s

#${pf} -pflotranin repository.in
# 31min

python3 "${src_dir}/pft_to_csv.py" repository-obs-0.pft

# move output files
mkdir -f repository_out
mv repository.out repository-int.dat repository-mas.dat repository-obs-0.pft repository-obs-0.csv repository.h5 repository_out

