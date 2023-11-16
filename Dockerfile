# use the ubuntu base image
FROM ubuntu:22.04

MAINTAINER Tobias Rausch rausch@embl.de

# install packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install tracy
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/tracy.git \
    && cd /opt/tracy/ \
    && make STATIC=1 all \
    && make install

# Multi-stage build
FROM alpine:latest
RUN apk add --no-cache bash
RUN mkdir -p /opt/tracy/bin
WORKDIR /opt/tracy/bin
COPY --from=0 /opt/tracy/bin/tracy .

# Workdir
WORKDIR /home

# Add Tracy to PATH
ENV PATH="/opt/tracy/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/bash"]
