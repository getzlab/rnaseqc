# Dockerfile for RNASeQC
FROM ubuntu:16.04
MAINTAINER Aaron Graubert

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && \
    cmake .. && make && make install && make clean && \
    cp /usr/local/lib/bamtools/libbamtools.a /usr/local/lib/
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/lib/bamtools

#RNASeQC
COPY src /opt/rnaseqc/src
COPY Makefile /opt/rnaseqc
COPY args.hxx /opt/rnaseqc
RUN cd /opt/rnaseqc && make && mv rnaseqc /usr/local/bin && make clean

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
