FROM ubuntu:20.04

LABEL org.opencontainers.image.name="TEEMS solver base layer Dockerfile"
LABEL org.opencontainers.image.version="0.91"
LABEL org.opencontainers.image.description="Base Layer Docker build for TEEMS solver"
LABEL org.opencontainers.image.author="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.maintainer="Matthew Cantele <matthew.cantele@protonmail.com>"
LABEL org.opencontainers.image.source="https://github.com/matthewcantele/teems-solver"
LABEL org.opencontainers.image.licenses="GPLv3.0"

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

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker \
    && echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker \
    && apt-get update \
    && apt-get install -y \
        automake \
        autoconf \
        libtool \
        build-essential \
        gfortran \
        wget \
        libopenblas-dev \
        git \
        python3 \
        python3-pip \
        python3-distutils \
        python3-setuptools \
        tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # Python libraries for bin_parser
    && pip3 install --no-cache-dir pandas toolz numpy Pyarrow \
    # OPENMPI
    && mkdir ${TEEMS_TMP} \
    && cd ${TEEMS_TMP} \
    && mkdir openmpi \
    && cd openmpi \
    && wget https://download.open-mpi.org/release/open-mpi/v${OPENMPI_BASE_VERSION}/openmpi-${OPENMPI_VERSION}.tar.gz \
    && tar xfz openmpi-${OPENMPI_VERSION}.tar.gz --strip-components=1 \
    && mkdir -p ${BUILD_DIR}/LICENSES/openmpi \
    && cp ${TEEMS_TMP}/openmpi/LICENSE ${BUILD_DIR}/LICENSES/openmpi \
    && mkdir -p ${BUILD_DIR}/lib/mpi \
    && ./configure --prefix=${BUILD_DIR}/lib/mpi \
    && make all install \
    # MPICH
    && cd ${TEEMS_TMP} \
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
    && mkdir ${BUILD_DIR}/LICENSES/petsc \ 
    && cp ${BUILD_DIR}/lib/petsc/LICENSE ${BUILD_DIR}/LICENSES/petsc \
    && ldconfig \
    && python3 ./config/configure.py PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=$PWD --with-debugging=0 --with-mpi-dir=${BUILD_DIR}/lib/mpi --with-afterimage COPTFLAGS="-Ofast -march=native -mtune=native" CXXOPTFLAGS="-Ofast -march=native -mtune=native" FOPTFLAGS="-Ofast -march=native -mtune=native" --with-shared-libraries=0 --with-openmp --with-threadsafety --with-log=0 --with-x=0 --force \
    && make PETSC_DIR=${PWD} PETSC_ARCH=${PETSC_ARCH} all \
    # HSL
    && mkdir -p ${BUILD_DIR}/LICENSES/hsl \
    && mkdir ${BUILD_DIR}/lib/ma48 \
    && mkdir ${BUILD_DIR}/lib/ma51 \
    && mkdir ${BUILD_DIR}/lib/mc66 \
    && mkdir ${BUILD_DIR}/lib/mp48 \
    # non-root user
    && groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # permissions for the bin_parser directory
    && chown -R 1000:1000 /opt/teems-parser \
    && chmod -R 755 /opt/teems-parser \
    # cleanup
    && rm -rf ${TEEMS_TMP} \
    && apt-get autoremove -y \
    && apt-get autoclean

USER $USERNAME
WORKDIR /opt/teems
CMD ["/bin/bash"]
