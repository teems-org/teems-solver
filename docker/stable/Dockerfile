# syntax=docker/dockerfile:1-labs

FROM ubuntu:22.04
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update \
  && apt-get install -y python3 \
  && rm -rf /var/lib/apt/lists/*

# linux packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    autoconf \
    automake \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config \
    libglib2.0-dev \
    libfuse3-dev \
    libtool \
    cryptsetup-bin \
    wget \
    libxml2 \
    libglpk-dev \
    acl \
    python3-pip \
    tzdata

ENV PATH="$PATH:/usr/local/go/bin"

# Go installation (Singularity dependency)
RUN cd /tmp \
    && wget https://go.dev/dl/go1.21.6.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go1.21.6.linux-amd64.tar.gz \
    && rm go1.21.6.linux-amd64.tar.gz \
    && echo 'export PATH=$PATH:/usr/local/go/bin' >> /etc/profile \
    && . /etc/profile

# Singularity installation
RUN export SINGULARITY_VERSION=4.0.0 \
    && wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce-$SINGULARITY_VERSION.tar.gz \
    && tar -xzf singularity-ce-$SINGULARITY_VERSION.tar.gz -C /opt \
    && rm singularity-ce-$SINGULARITY_VERSION.tar.gz \
    && cd /opt/singularity-ce-$SINGULARITY_VERSION \
    && ./mconfig \
    && make -C ./builddir \
    && make -C ./builddir install

# Python libraries for bin_parser
RUN pip3 install --no-cache-dir pandas toolz numpy Pyarrow

# Inner layer: Singularity image build
COPY ./singularity/latest/Singularity.intel /tmp
COPY ./solver /tmp/solver
COPY ./tarballs /tmp/tarballs
COPY ./bin_parser /home/bin_parser

RUN --security=insecure singularity build /home/teems_solver.sif /tmp/Singularity.intel \
    && mkdir -p /home/launchpad

WORKDIR /home/launchpad
CMD ["/init"]

# build with:
#docker buildx build --load --allow security.insecure -t teems:beta -f ./docker/latest/Dockerfile .
