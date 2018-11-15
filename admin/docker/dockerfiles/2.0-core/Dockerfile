# pull ubuntu 18.04 from Docker repo
FROM ubuntu:18.04

MAINTAINER "Zhicheng Geng <zhichenggeng@utexas.com>"

# install packages for madagascar
RUN apt-get update && apt-get install -y \
        git \
        python2.7 \
        python-pip \
        openssh-client \
        tar \
        gzip \
        wget \
        vim \
        emacs \
        make \
        man \
 && apt-get install -y \
        libblas-dev \
        liblapack-dev \
        swig \
        libxaw7-dev \
        freeglut3-dev \
        libnetpbm10-dev \
        libtiff5-dev \
        libgd-dev \
        libplplot-dev \
        libavcodec-dev \
        libcairo2-dev \
        libjpeg-dev \
        libopenmpi-dev \
        libfftw3-dev \
        libsuitesparse-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# install python packages
RUN pip install numpy scipy

# get code from sourceforge
RUN  wget -O madagascar-core-2.0.tar.gz https://sourceforge.net/projects/rsf/files/madagascar/madagascar-2.0/madagascar-core-2.0.tar.gz/download \
 && tar -xzvf madagascar-core-2.0.tar.gz -C $HOME \
 && rm madagascar-core-2.0.tar.gz

# set environment variable for installing madagascar
ENV RSFROOT /root/RSFROOT

# install madagascar
RUN cd ~/madagascar-core-2.0 \
 && ./configure \
 && make install

RUN echo 'export RSFROOT="$HOME/RSFROOT"'                   >> $HOME/.bashrc \
 && echo 'source $RSFROOT/share/madagascar/etc/env.sh'      >> $HOME/.bashrc

WORKDIR /root

CMD ["/bin/bash"]
