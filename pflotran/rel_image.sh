docker build \
--build-arg source_repository=https://bitbucket.org/pflotran/pflotran \
--build-arg commit=083857c \
-t flow123d/pflotran-gnu-rel ./rel-image

