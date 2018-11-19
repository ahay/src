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

ENV LANG C.UTF-8

# install python packages
RUN pip install numpy scipy jupyter

# get code from github
RUN git clone https://github.com/ahay/src.git $HOME/RSFSRC

# set environment variable for installing madagascar
ENV RSFROOT /root/RSFROOT

# install madagascar
RUN cd ~/RSFSRC \
 && ./configure \
 && make install

# install latex
RUN apt-get update && apt-get install -y \
        texlive-latex-recommended \
		texlive-latex-extra \
		texlive-fonts-recommended \
		texlive-bibtex-extra \
		texlive-lang-english \
		texlive-generic-extra \
		biber \
        --no-install-recommends \
 && rm -rf /var/lib/apt/lists/*

# install segtex
RUN git clone https://github.com/SEGTeX/texmf $HOME/texmf \
 && texhash

RUN echo 'export RSFROOT="$HOME/RSFROOT"'                   >> $HOME/.bashrc \
 && echo 'source $RSFROOT/share/madagascar/etc/env.sh'      >> $HOME/.bashrc 

WORKDIR /root
EXPOSE 8888

CMD ["bash", "-c" ,"export SHELL=/bin/bash && source $RSFROOT/share/madagascar/etc/env.sh && jupyter notebook --notebook-dir=/root --ip 0.0.0.0 --no-browser --allow-root"]
