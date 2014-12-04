these directories are for processng a 2d extraction of the 3D SEAM project

SEAM stands for SEG Advanced Modeling project.

SEAM Phase1 is a 3d synthetic for benchmarking subsalt techniques.  I mostly 
know about the weve equation seismic synthetic data, but there is also
synthetic gravity.

The 2D extraction is in the [rpcess of release Summer of 2014.  I wonder
if formal release will be at the 2014 SEG convention in October.

This data was used at the 2d Madagascar Working Workshop Summer 2014.
Our objective was (is) to compare 2d finite difference wave equation programs
in Madagascar.  Statements about relative runtime and accuracy would be great.
 
Directories in this test are:

fetch - download the 2D extraction of velocity and 3D seismic
sfmshots - compute 2D synthetic using Pengliang's 4th order space 3rd order time
sfawefd3d - compute 2D synthetic using Sava's 4th order space 3rd order time
sfmshots1 - run 151 shots using pscons parallel on 151 cores, 10 nodes
sfawefd2d_split - run 151 shots using split and pscons on 151 cores, 10 nodes
ntg
