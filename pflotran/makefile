pflotran_build:
	./fterm.sh make -C ./pflotran_JB/src/pflotran pflotran 

#repo=git@github.com:flow123d/pflotran.git
repo=https://github.com/flow123d/pflotran.git
#commit=083857c
commit=main
NO_CACHE=--no-cache
build_image:
	docker build $(NO_CACHE) --build-arg source_repository=$(repo) --build-arg commit=$(commit) \
	-t flow123d/pflotran-gnu-rel ./rel-image

	
gid:=$(shell id -g)
uid:=$(shell id -u)	
uname=user
curr_dir:=$(shell pwd)
IMG_WORKDIR:=/home/$(uname)/workdir
envargs=-euid=$(uid) -egid=$(gid) -ewho=$(uname) -ehome=${IMG_WORKDIR}
mountargs=-v $(curr_dir):$(IMG_WORKDIR)
docker_call=

eos_test_img:
	export PFLOTRAN_CMD="docker run --rm -it $(envargs) $(mountargs) -w $(IMG_WORKDIR)/eos_test/out flow123d/pflotran-gnu-rel"; \
	cd eos_test && ./run.sh eos_test
