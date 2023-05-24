#!/bin/bash
set -x

# Script for entering a build docker image to be used with the `pflotran_JB` submodule.
# Only meant for development of the specific sources: reaction sandbox files. 
# Final release image should be build via `release-image`
ABS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Simple docker 
gid=$(id -g)
uid=$(id -u)
uname="user"
IMG_WORKDIR="/home/${uname}/workdir"

envargs="-euid=$uid -egid=$gid -ewho=$uname -ehome=${IMG_WORKDIR} -ePETSC_DIR=/usr/local/petsc_v3.18.6/"
mountargs="-w ${IMG_WORKDIR} -v ${ABS_SCRIPT_DIR}:${IMG_WORKDIR}"
docker run --rm -it $envargs $mountargs flow123d/flow-dev-gnu-rel:4.0.3
