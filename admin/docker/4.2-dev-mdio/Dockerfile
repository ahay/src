# Madagascar 4.2 with MDIO (mdio-cpp / TGS) support
# Base: ubuntu:26.04
#
# Build from the Madagascar source root (docker branch):
#   docker build \
#     -f admin/docker/4.2-dev-mdio/Dockerfile \
#     -t madagascar/m8r:4.2-dev-mdio \
#     .
#
# The local source tree is COPYed (not cloned from GitHub) because the MDIO
# programs (sfmdioread, sfmdiowrite) live on the 'docker' branch and are not
# yet merged upstream.
#
# NOTE: mdio-cpp compilation (tensorstore + absl + ~100 other deps) takes
# 20-30 minutes on a machine with 8+ cores.

FROM ubuntu:26.04

LABEL maintainer="Sergey Fomel <sergey.fomel@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Install packages for Madagascar + mdio-cpp build dependencies.
# Notes vs. the 24.04 reference Dockerfile:
#   - libtiff-dev replaces libtiff5-dev (libtiff5-dev is transitional on 26.04)
#   - nasm is required by tensorstore (used in mdio-cpp)
#   - cmake ships as 3.30+ on 26.04, satisfying the mdio-cpp 3.24 minimum
#   - build-essential provides g++ needed by mdio-cpp C++ sources
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
        build-essential \
        nasm \
        cmake \
 && apt-get install -y \
        libblas-dev \
        liblapack-dev \
        swig \
        python-dev-is-python3 \
        python3-numpy \
        libxaw7-dev \
        freeglut3-dev \
        libnetpbm11-dev \
        libtiff-dev \
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
# setuptools is required because Python 3.12+ removed the 'distutils' module
# from the standard library; Madagascar's SCons scripts import distutils.sysconfig
RUN pip install --break-system-packages numpy scipy setuptools

# Build mdio-cpp (header-only library, but tensorstore static archives must
# be compiled before Madagascar can link against them).
RUN git clone https://github.com/TGSAI/mdio-cpp.git /root/mdio-cpp

# Configure: download FetchContent deps and generate build files.
# BUILD_TESTING=OFF suppresses test registration but specific targets are
# still built below to compile the tensorstore/absl archives.
# GCC 15 on Ubuntu 26.04 needs two workarounds for older vendored code:
#   CMAKE_C_STANDARD=17 — C23 makes 'bool' a keyword, breaking blosc's
#                          'typedef _Bool bool;'
#   CMAKE_CXX_FLAGS="-include cstdint" — GCC 15's -Wtemplate-body now checks
#                          template definitions; abseil lts_20240116 uses
#                          uintptr_t without including <cstdint>
RUN cmake -B /root/mdio-cpp/build \
      -S /root/mdio-cpp \
      -DCMAKE_INSTALL_PREFIX=/root/mdio-cpp-install \
      -DBUILD_TESTING=OFF \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_C_STANDARD=17 \
      -DCMAKE_CXX_FLAGS="-include cstdint"

# Compile the tensorstore / absl / S3 static archives (20-30 min).
# mdio_dataset_test: local filesystem kvstore drivers (zarr, file).
# mdio_s3_test:      S3 kvstore — needed for sfmdioread mdio=s3://... URIs.
RUN cmake --build /root/mdio-cpp/build \
      --target mdio_dataset_test \
      -j"$(nproc)"

RUN cmake --build /root/mdio-cpp/build \
      --target mdio_s3_test \
      -j"$(nproc)"

# Build the JSON schema validator archive installed by cmake --install.
RUN cmake --build /root/mdio-cpp/build \
      --target nlohmann_json_schema_validator \
      -j"$(nproc)"

# Install mdio headers and libnlohmann_json_schema_validator.a.
# The tensorstore .a archives remain in the build tree; configure.py reads
# their paths directly from link.txt (requires MDIO_BUILD_DIR at configure time).
RUN cmake --install /root/mdio-cpp/build

# Copy the local Madagascar source (docker branch with sfmdioread/sfmdiowrite).
COPY . /root/RSFSRC

# set environment variable for installing Madagascar
ENV RSFROOT=/root/RSFROOT

# GCC 15 treats -Wincompatible-pointer-types as an error by default (C23
# constraint). Some Madagascar C files have pre-existing int**/float**
# mismatches that are harmless at runtime.  Append the downgrade flag to
# every gcc flag string that configure.py tests.
RUN sed -i "s/-Wall -pedantic/-Wall -pedantic -Wno-error=incompatible-pointer-types/g" \
    /root/RSFSRC/framework/configure.py

# Configure and build Madagascar with MDIO support.
# MDIO_HOME:      install prefix (headers + libnlohmann_json_schema_validator.a)
# MDIO_BUILD_DIR: cmake build tree (configure.py reads link.txt from here)
WORKDIR /root/RSFSRC
RUN rm -f config.py \
 && ./configure \
      --prefix="${RSFROOT}" \
      MDIO_HOME=/root/mdio-cpp-install \
      MDIO_BUILD_DIR=/root/mdio-cpp/build \
 && make install

RUN echo 'export RSFROOT="/root/RSFROOT"'                   >> "$HOME/.bashrc" \
 && echo 'source $RSFROOT/share/madagascar/etc/env.sh'      >> "$HOME/.bashrc"

WORKDIR /root

CMD ["/bin/bash"]
