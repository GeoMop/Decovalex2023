ARG source_image
ARG base_image
FROM flow123d/flow-dev-gnu-rel:4.0.3 as build_cont

ENV PETSC_DIR=/usr/local/petsc_v3.18.6/
WORKDIR ${HOME}/workspace/
RUN git clone https://bitbucket.org/pflotran/pflotran
# Last merge before update ro PETSC 3.19
WORKDIR ${HOME}/workspace/pflotran
RUN git checkout 083857c
WORKDIR ${HOME}/workspace/pflotran/src/pflotran
COPY reaction_sandbox_*.F90 pflotran_*.txt .
RUN make  pflotran
ENTRYPOINT ["/workspace/pflotran/src/pflotran/pflotran"]


# ARG base_image
# FROM ${base_image}
# 
# ARG base_image
# ARG source_location
# ARG flow_version
# ARG flow_install_location=/opt/flow123d
# ARG git_hash
# ARG build_date
# 
# # Build-time metadata as defined at http://label-schema.org
# LABEL org.label-schema.build-date="$build_date" \
#       org.label-schema.name="Flow123d" \
#       org.label-schema.description="Transport Processes in Fractured Media" \
#       org.label-schema.url="http://flow123d.github.io/" \
#       org.label-schema.vcs-ref="$git_hash" \
#       org.label-schema.vcs-url="https://github.com/flow123d/flow123d" \
#       org.label-schema.vendor="Technical University of Liberec" \
#       org.label-schema.version="$flow_version" \
#       org.label-schema.schema-version="1.0"
#       
#       
# COPY --from=build_cont ${source_location} ${flow_install_location} 
# 
# # update PATH
# ENV PATH="${flow_install_location}/bin:${PATH}"
# 
# 
# # save the image version
# #RUN echo ${flow_version} > /.dockerversion
