#!/bin/bash

# Converts "old"-style SEPlib complex files (esize=8, data_format=xdr_float) to
# RSF (esize=8, data_format=xdr_complex).
#
# Usage: sepcmplx2rsf.sh directory
# 
# Processes correctly files in which the binary is merged with the header, and
# does not change the timestamp of any file

for y in `find $1 -name "*.H"`;
  do 
    sfsepcmplx2rsf file=$y verb=y;
  done
