#!/bin/bash
set -x
# Run the eos_test

# expect PFLOTRAN_CMD variable
# Example configurations:
# # Use a build in a development image.
#
# PFLOTRAN_CMD=../pflotran_JB/src/pflotran/pflotran
#
# # Use docker container.
# PFLOTRAN_CMD=docker run --rm -it -v `pwd`:/home/WORK -w /home/WORK/out flow123d/pflotran-gnu-rel 


prefix=$1
shift

python3 prepare_fields.py
# No obvious way how to separate output from input, output is generated at the same directory as the input file.
# So we rather copy the input files to the output to keep them separated.
rm -rf out
mkdir -p out
cp * out

#
$PFLOTRAN_CMD -input_prefix ${prefix} $@ 
# Redirection of STDOUT seems not necessary as the *.out file contains the same with more detailed information about initial and boundary conditions.
# | tee out/stdout
