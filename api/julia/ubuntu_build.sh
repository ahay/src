#!/bin/bash

# Check if we are on ubuntu
lsb_release -a | grep "Ubuntu" || exit 1
UBUNTU_VERSION=$(lsb_release -r | awk '{print $2}')

# Get the git hash
GIT_HASH=$(git rev-parse HEAD)

cd ../..

PREFIX=$(pwd)/__PRECOMPILED_BUILD_UBUNTU_${UBUNTU_VERSION}__${GIT_HASH}__
INCLUDE_PATH=$PREFIX/include
LIB_PATH=$PREFIX/lib
./configure \
    --prefix=$PREFIX \
    CC=$(which gcc) \
    CXX=$(which g++) \
    CAIROPATH=$INCLUDE_PATH/cairo \
    FFMPEGPATH=$INCLUDE_PATH/ffmpeg \
    AR=$(which ar) \
    CPPPATH=$INCLUDE_PATH \
    LIBPATH=$LIB_PATH \
    BLAS=openblas

# Build
NUM_CPUS=$(lscpu | grep "^CPU(s):" | grep -v "NUMA" | awk '{print $2}')
make -j$NUM_CPUS

# Install
make install