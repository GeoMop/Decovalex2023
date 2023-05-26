#!/bin/bash

#repo=https://bitbucket.org/pflotran/pflotran
repo=git@github.com:flow123d/pflotran.git

#commit=083857c
commit=main

docker build \
--build-arg source_repository=$repo \
--build-arg commit=$commit \
-t flow123d/pflotran-gnu-rel ./rel-image

