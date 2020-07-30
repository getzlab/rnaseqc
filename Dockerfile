# Dockerfile for RNASeQC
FROM ubuntu:20.04
MAINTAINER Aaron Graubert

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        python3 \
        python3-pip \
        libboost-filesystem-dev \
        libboost-regex-dev \
        libboost-system-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libpthread-stubs0-dev \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
    
# Python
RUN python3 -m pip install --upgrade pip setuptools pyarrow jupyter

# SeqLib
COPY Makefile /opt/rnaseqc/Makefile
RUN cd /opt/rnaseqc && git clone --recursive https://github.com/walaj/SeqLib.git && \
    cd SeqLib && git checkout 7e1f98267b5057f9505dbff119308137a0e006db && cd .. && \
    make SeqLib/lib/libseqlib.a

# python
RUN cd /opt && git clone https://github.com/francois-a/rnaseq-utils rnaseq && cd rnaseq && \
  git checkout f1c6a5677bbca465ea1edd06c2293a5d1078a18b && python3 -m pip install --upgrade pip setuptools && \
  python3 -m pip install numpy && python3 -m pip install pandas matplotlib scipy pyBigWig bx-python \
  agutil nbformat seaborn sklearn qtl && mkdir -p /root/.config/matplotlib && echo "backend	:	Agg" > /root/.config/matplotlib/matplotlibrc
ENV PYTHONPATH $PYTHONPATH:/opt/

#RNASeQC
COPY src /opt/rnaseqc/src
COPY python /opt/rnaseqc/python
COPY args.hxx /opt/rnaseqc
COPY bioio.hpp /opt/rnaseqc
RUN cd /opt/rnaseqc && make && ln -s /opt/rnaseqc/rnaseqc /usr/local/bin/rnaseqc && make clean && python3 -m pip install -e /opt/rnaseqc/python

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
