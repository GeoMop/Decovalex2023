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
HOST_ROOT_MOUNT="/$ABS_SCRIPT_DIR"
IMG_ROOT_MOUNT="/home/${uname}/workdir/"
HOST_PWD="/`pwd`"
IMG_PWD="$IMG_ROOT_MOUNT/${HOST_PWD##$HOST_ROOT_MOUNT/}"

envargs="-euid=$uid -egid=$gid -ewho=$uname -ehome=${IMG_ROOT_MOUNT} -ePETSC_DIR=//usr/local/petsc_v3.18.6/"
mountargs="-w /${IMG_PWD} -v ${HOST_ROOT_MOUNT}:${IMG_ROOT_MOUNT}"
interactive="--rm -it"
#image=flow123d/flow-dev-gnu-rel:4.0.3
image=${PFLOTRAN_IMAGE:-flow123d/pflotran-gnu-rel}
docker_call="docker run $interactive  $envargs $mountargs $image $cmd"

#docker_call=docker run $interactive  $envargs $mountargs flow123d/flow-dev-gnu-rel:4.0.3 $cmd

if [ -n "$1" ]
then        
    # command
    $docker_call $@
else
    # shell
    $docker_call
fi
