FROM ubuntu:20.04

LABEL org.label-schema.name="TEEMS solver Dockerfile"
LABEL org.label-schema.version="0.9"
LABEL org.label-schema.description="Full Docker build for TEEMS solver"
LABEL org.label-schema.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.label-schema.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.label-schema.source="https://github.com/matthewcantele/teems-solver"
LABEL org.label-schema.licenses="GPLv3.0"

ARG PATH_HSL_MA48
ARG PATH_HSL_MA51
ARG PATH_HSL_MC66
ARG PATH_HSL_MP48
ARG BUILD_DIR=/home/teems-solver
ARG OPENMPI_VERSION=5.0.7
ARG OPENMPI_BASE_VERSION=5.0
ARG MPICH_VERSION=4.3.0
ARG PETSC_ARCH=native
ARG USERNAME=teems
ARG USER_UID=1000
ARG USER_GID=$USER_UID

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="${BUILD_DIR}/lib/mpi/bin:$PATH"

COPY ${PATH_HSL_MA48} /tmp/hsl/ma48.tar.gz
COPY ${PATH_HSL_MA51} /tmp/hsl/ma51.tar.gz
COPY ${PATH_HSL_MC66} /tmp/hsl/mc66.tar.gz
COPY ${PATH_HSL_MP48} /tmp/hsl/mp48.tar.gz
COPY ./src ${BUILD_DIR}/src
COPY ./bin_parser /home/bin_parser

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker \
    && echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker \
    && apt-get update \
    && apt-get install -y \
        automake \
        build-essential \
        gfortran \
        wget \
        libopenblas-dev \
        git \
        python3 \
        python3-pip \
        python3-distutils \
        python-setuptools \
        tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # Python libraries for bin_parser
    && pip3 install --no-cache-dir pandas toolz numpy Pyarrow \
    # libraries directory
    && mkdir -p ${BUILD_DIR}/tmp \
    ${BUILD_DIR}/lib/ucx \
    ${BUILD_DIR}/lib/mpi \
    ${BUILD_DIR}/lib/ma48 \
    ${BUILD_DIR}/lib/ma51 \
    ${BUILD_DIR}/LICENSES/openmpi \
    ${BUILD_DIR}/LICENSES/petsc \
    ${BUILD_DIR}/LICENSES/hsl \
    # unzip HSL libraries
    && mkdir -p ${BUILD_DIR}/tmp/hsl/ma48 \
    && tar xfz /tmp/hsl/ma48.tar.gz --strip-components=1 -C ${BUILD_DIR}/tmp/hsl/ma48 \
    && mkdir -p ${BUILD_DIR}/tmp/hsl/ma51 \
    && tar xfz /tmp/hsl/ma51.tar.gz --strip-components=1 -C ${BUILD_DIR}/tmp/hsl/ma51 \
    && mkdir -p ${BUILD_DIR}/tmp/hsl/mc66 \
    && tar xfz /tmp/hsl/mc66.tar.gz --strip-components=1 -C ${BUILD_DIR}/tmp/hsl/mc66 \
    && mkdir -p ${BUILD_DIR}/tmp/hsl/mp48 \
    && tar xfz /tmp/hsl/mp48.tar.gz --strip-components=1 -C ${BUILD_DIR}/tmp/hsl/mp48 \
    # OPENMPI
    && cd ${BUILD_DIR}/tmp \
    && mkdir openmpi \
    && cd openmpi \
    && wget https://download.open-mpi.org/release/open-mpi/v${OPENMPI_BASE_VERSION}/openmpi-${OPENMPI_VERSION}.tar.gz \
    && tar xfz openmpi-${OPENMPI_VERSION}.tar.gz --strip-components=1 \
    && cp ${BUILD_DIR}/tmp/openmpi/LICENSE ${BUILD_DIR}/LICENSES/openmpi \
    && ./configure --prefix=${BUILD_DIR}/lib/mpi \
    && make all install \
    # MPICH
    && cd ${BUILD_DIR}/tmp \
    && mkdir mpich \
    && cd mpich \
    && wget https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar xfz mpich-${MPICH_VERSION}.tar.gz --strip-components=1 \
    && ./configure --prefix=${BUILD_DIR}/lib/mpi \
    && make \
    && make install \
    # PETSC
    && cd ${BUILD_DIR}/lib \
    && git clone -b release https://gitlab.com/petsc/petsc.git petsc \
    && cd petsc \
    && cp ${BUILD_DIR}/lib/petsc/LICENSE ${BUILD_DIR}/LICENSES/petsc \
    && ldconfig \
    && pip3 install --no-cache-dir pandas toolz \
    && python3 ./config/configure.py PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=$PWD --with-debugging=0 --with-mpi-dir=${BUILD_DIR}/lib/mpi --with-afterimage COPTFLAGS="-Ofast -march=native -mtune=native" CXXOPTFLAGS="-Ofast -march=native -mtune=native" FOPTFLAGS="-Ofast -march=native -mtune=native" --with-shared-libraries=0 --with-openmp --with-threadsafety --with-log=0 --with-x=0 --force \
    && make PETSC_DIR=${PWD} PETSC_ARCH=${PETSC_ARCH} all \
    # HSL
    && cp -r ${BUILD_DIR}/tmp/hsl/ma48 ${BUILD_DIR}/lib \
    && cd ${BUILD_DIR}/lib/ma48 \
    && cp ${BUILD_DIR}/lib/ma48/LICENCE ${BUILD_DIR}/LICENSES/hsl \
    && ./configure \
    && make \
    && cp -r ${BUILD_DIR}/tmp/hsl/ma51 ${BUILD_DIR}/lib \
    && cd ${BUILD_DIR}/lib/ma51 \
    && ./configure \
    && make \
    && cp -r ${BUILD_DIR}/tmp/hsl/mp48 ${BUILD_DIR}/lib \
    && cd ${BUILD_DIR}/lib/mp48 \
    && ./configure \
    && make \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp48d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/ddeps.f ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/lib/mp48/src/hsl_mp01.mod ${BUILD_DIR}/src \
    && cd ${BUILD_DIR}/src \
    && chmod +x mp48_mod.sh \
    && ./mp48_mod.sh \
    && cp ${BUILD_DIR}/tmp/hsl/mc66/src/hsl_mc66d.f90 ${BUILD_DIR}/src \
    && cp ${BUILD_DIR}/tmp/hsl/mc66/src/ddeps90.f90 ${BUILD_DIR}/src \
    # final build
    && cd ${BUILD_DIR}/src \
    && make \
    && mv ${BUILD_DIR}/src ${BUILD_DIR}/solver \
    && rm -rf ${BUILD_DIR}/tmp \
    # non-root user
    && groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # fix permissions for the bin_parser directory
    && chown -R 1000:1000 /home/bin_parser && \
    chmod -R 755 /home/bin_parser

USER $USERNAME
WORKDIR /home/launchpad
CMD ["/bin/bash"]
