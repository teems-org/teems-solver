# syntax=docker/dockerfile:1-labs

FROM ubuntu:22.04
ARG PATH_HSL_MA48
ARG PATH_HSL_MA51
ARG PATH_HSL_MC66
ARG PATH_HSL_MP48

ENV DEBIAN_FRONTEND=noninteractive \
    PATH="$PATH:/usr/local/go/bin" \
    SINGULARITY_VERSION=4.0.0 \
    GO_VERSION=1.24.2 \
    RELEASE="dev" \
    PATH="$PATH:/usr/local/go/bin"

RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker \
    && echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker \
    && apt-get update \
    && apt-get install -y \
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
        python3 \
        python3-pip \
        tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Go installation (Singularity dependency)
RUN cd /tmp \
    && wget https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go${GO_VERSION}.linux-amd64.tar.gz \
    && rm go${GO_VERSION}.linux-amd64.tar.gz \
    && echo 'export PATH=$PATH:/usr/local/go/bin' >> /etc/profile \
    && . /etc/profile

# Singularity installation
RUN wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce-${SINGULARITY_VERSION}.tar.gz \
    && tar -xzf singularity-ce-${SINGULARITY_VERSION}.tar.gz -C /opt \
    && rm singularity-ce-${SINGULARITY_VERSION}.tar.gz \
    && cd /opt/singularity-ce-${SINGULARITY_VERSION} \
    && ./mconfig \
    && make -C ./builddir \
    && make -C ./builddir install

# Python libraries for bin_parser
RUN pip3 install --no-cache-dir pandas toolz numpy Pyarrow

# HSL libraries
COPY ${PATH_HSL_MA48} /tmp/hsl/ma48.tar.gz
COPY ${PATH_HSL_MA51} /tmp/hsl/ma51.tar.gz
COPY ${PATH_HSL_MC66} /tmp/hsl/mc66.tar.gz
COPY ${PATH_HSL_MP48} /tmp/hsl/mp48.tar.gz

RUN mkdir -p /tmp/hsl/ma48 \
    && tar xfz /tmp/hsl/ma48.tar.gz --strip-components=1 -C /tmp/hsl/ma48 \
    && mkdir -p /tmp/hsl/ma51 \
    && tar xfz /tmp/hsl/ma51.tar.gz --strip-components=1 -C /tmp/hsl/ma51 \
    && mkdir -p /tmp/hsl/mc66 \
    && tar xfz /tmp/hsl/mc66.tar.gz --strip-components=1 -C /tmp/hsl/mc66 \
    && mkdir -p /tmp/hsl/mp48 \
    && tar xfz /tmp/hsl/mp48.tar.gz --strip-components=1 -C /tmp/hsl/mp48

# Inner layer: Singularity image build (copy from local machine into docker container)
COPY ./singularity/${RELEASE}/Singularity.intel /tmp
COPY ./solver /tmp/solver
COPY ./tarballs /tmp/tarballs
COPY ./bin_parser /home/bin_parser

RUN --security=insecure singularity build /home/teems_solver.sif /tmp/Singularity.intel \
    && mkdir -p /home/launchpad

WORKDIR /home/launchpad
CMD ["/init"]

# build with:
#docker buildx build --build-arg PATH_HSL_MA48="tarballs/hsl/ma48-2.2.0.tar.gz" --build-arg PATH_HSL_MA51="tarballs/hsl/ma51-1.0.0.tar.gz" --load --allow security.insecure -t teems:beta -f ./docker/development/Dockerfile .
