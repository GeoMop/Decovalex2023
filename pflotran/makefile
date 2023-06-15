
# ======================= DEveloplment build
pflotran_build:
	export PFLOTRAN_IMAGE="flow123d/flow-dev-gnu-rel:4.0.3"; ./fterm.sh make -j5 -C ./pflotran_JB/src/pflotran pflotran 

	

# ======================= EOS diffusion model release image
#repo=git@github.com:flow123d/pflotran.git
repo=https://github.com/flow123d/pflotran.git
#commit=083857c
commit=main
NO_CACHE=--no-cache

rel_image:
	docker build $(NO_CACHE) --build-arg source_repository=$(repo) --build-arg commit=$(commit) \
	-t flow123d/pflotran-gnu-rel ./rel-image

# ======================= EOS test	

dwp_run:
	cd dwp_test && ../fterm.sh python3 prepare_fields.py
	cd dwp_test; rm -rf out; mkdir -p out; cp * out; cd out
	cd dwp_test/out && ../fterm.sh  


dwp_dev:
	cd dwp_test; \
	export PFLOTRAN_IMAGE="flow123d/flow-dev-gnu-rel:4.0.3"; \
	python3 prepare_fields.py; \
	rm -rf out; mkdir -p out; cp * out; cd out; \
	../../fterm.sh ../../pflotran_JB/src/pflotran/pflotran  
