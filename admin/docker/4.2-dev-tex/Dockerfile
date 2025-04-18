# pull ubuntu 24.04 from Docker repo
FROM ubuntu:24.04

MAINTAINER "Sergey Fomel <sergey.fomel@gmail.com>"

# install packages for madagascar
RUN apt-get update && apt-get install -y \
        git \
        python3-full \
        python-is-python3 \
        python3-pip \
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
        python-dev-is-python3 \
        python3-numpy \
        libxaw7-dev \
        freeglut3-dev \
        libnetpbm11-dev \
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

RUN echo 'export RSFROOT="/root/RSFROOT"'                   >> $HOME/.bashrc \
 && echo 'source $RSFROOT/share/madagascar/etc/env.sh'      >> $HOME/.bashrc

WORKDIR /root

CMD ["/bin/bash"]
