
# ======================= DEveloplment build
pflotran_build:
	./fterm.sh make -C ./pflotran_JB/src/pflotran pflotran 

	

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


eos_test_img:
	cd eos_test && ../run.sh -c python3 prepare_fields.py
	cd eos_test && ../run.sh 
