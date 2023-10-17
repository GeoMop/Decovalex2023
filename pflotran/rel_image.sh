#!/bin/bash

#repo=https://bitbucket.org/pflotran/pflotran
#repo=git@github.com:flow123d/pflotran.git
repo=https://github.com/flow123d/pflotran

# commit or branch to build
#commit=083857c
#commit=main
commit=OM_decay_source

docker build \
--build-arg source_repository=$repo \
--build-arg commit=$commit \
-t flow123d/pflotran-gnu-rel:OM_decay_source ./rel-image

