# Dockerfile for RNASeQC
FROM ubuntu:16.04
MAINTAINER Aaron Graubert

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        git \
        python3 \
        python3-pip \
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

# python
RUN cd /opt && git clone https://github.com/francois-a/rnaseq-utils rnaseq && cd rnaseq && \
  git checkout f1c6a5677bbca465ea1edd06c2293a5d1078a18b && python3 -m pip install --upgrade pip setuptools && \
  python3 -m pip install numpy && python3 -m pip install pandas matplotlib scipy pyBigWig bx-python \
  agutil nbformat seaborn sklearn && mkdir -p /root/.config/matplotlib && echo "backend	:	Agg" > /root/.config/matplotlib/matplotlibrc
ENV PYTHONPATH $PYTHONPATH:/opt/

#RNASeQC
COPY src /opt/rnaseqc/src
COPY python /scripts
COPY Makefile /opt/rnaseqc
COPY args.hxx /opt/rnaseqc
COPY bioio.hpp /opt/rnaseqc
RUN cd /opt/rnaseqc && make && ln -s /opt/rnaseqc/rnaseqc /usr/local/bin/rnaseqc && make clean

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
