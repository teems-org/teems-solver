FROM ubuntu:20.04

LABEL org.label-schema.name="TEEMS solver base layer Dockerfile"
LABEL org.label-schema.version="0.91"
LABEL org.label-schema.description="Base Layer Docker build for TEEMS solver"
LABEL org.label-schema.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.label-schema.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.label-schema.source="https://github.com/matthewcantele/teems-solver"
LABEL org.label-schema.licenses="GPLv3.0"

ARG BUILD_DIR=/opt/teems-solver
ARG OPENMPI_VERSION=5.0.7
ARG OPENMPI_BASE_VERSION=5.0
ARG MPICH_VERSION=4.3.0
ARG PETSC_ARCH=native
ARG USERNAME=teems
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG TEEMS_TMP=/tmp/teems

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"

COPY ./src ${BUILD_DIR}/src
COPY ./teems-parser /opt/teems-parser

RUN cd ${TEEMS_TMP} \
    && tar xfz ma48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma48 \
    && cp ${BUILD_DIR}/lib/ma48/LICENCE ${BUILD_DIR}/LICENSES/hsl \
    && tar xfz ${TEEMS_TMP}/ma51.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/ma51 \
    && tar xfz ${TEEMS_TMP}/mc66.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mc66 \
    && tar xfz ${TEEMS_TMP}/mp48.tar.gz --strip-components=1 -C ${BUILD_DIR}/lib/mp48 \
    && cd ${BUILD_DIR}/lib/ma48 \
    && autoreconf -fi \
    && ./configure \
    && make \
    && ln -s ${BUILD_DIR}/lib/ma48/src/.libs/libma48.so.0.0.0 /usr/local/lib/libma48.so \
    && ln -s ${BUILD_DIR}/lib/ma48/src/.libs/libma48.so.0.0.0 /usr/local/lib/libma48.so.0 \
    && cd ${BUILD_DIR}/lib/ma51 \
    && autoreconf -fi \
    && ./configure \
    && make \
    && ln -s ${BUILD_DIR}/lib/ma51/src/.libs/libma51.so.0.0.0 /usr/local/lib/libma51.so \
    && ln -s ${BUILD_DIR}/lib/ma51/src/.libs/libma51.so.0.0.0 /usr/local/lib/libma51.so.0 \
    && cd ${BUILD_DIR}/lib/mp48 \
    && autoreconf -fi \
    && ./configure \
    && make \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp48d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/ddeps.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp01.mod ${BUILD_DIR}/src \
    && cd ${BUILD_DIR}/src \
    && chmod +x mp48_mod.sh \
    && ./mp48_mod.sh \
    && cp ${BUILD_DIR}/lib/mc66/src/hsl_mc66d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mc66/src/ddeps90.f90 ${BUILD_DIR}/src \
    && ldconfig \
    # final build
    && cd ${BUILD_DIR}/src \
    && make \
    && mv ${BUILD_DIR}/src ${BUILD_DIR}/solver \
    && rm -rf ${TEEMS_TMP} \
    # non-root user
    && groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # permissions for the bin_parser directory
    && chown -R 1000:1000 /opt/teems-parser \
    && chmod -R 755 /opt/teems-parser 

USER $USERNAME
WORKDIR /opt/teems
CMD ["/bin/bash"]
