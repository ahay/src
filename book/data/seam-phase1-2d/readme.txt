these directories are for processng a 2d extraction of the 3D SEAM project

SEAM stands for SEG Advanced Modeling project.

SEAM Phase1 is a 3d synthetic for benchmarking subsalt techniques.  I mostly 
know about the weve equation seismic synthetic data, but there is also
synthetic gravity.

The 2D extraction was release Summer of 2014.

This data was used at the 2d Madagascar Working Workshop Summer 2014.
Our objective was to compare 2d finite difference wave equation programs
in Madagascar.  Statements about relative runtime and accuracy would be great.
 
Directories in this test are:

fetch - download and dipslay the pdf's that describe the data.  Fetch and 
        display the seam 3D elastic shot records, vp, vs, density and wavelets.
sfmshots - compute 2D synthetic using Pengliang's 4th order space 3rd order time
sfawefd3d - compute 2D synthetic using Sava's 4th order space 3rd order time
sfmshots1 - run 151 shots using pscons parallel on 151 cores, 10 nodes
sfawefd2d_split - run 151 shots using split and pscons on 151 cores, 10 nodes

These directories download, display and apply basic processing to the 2D 
subset of the 3D elastic modeling:
fetch - the directory described above must be run to down load data for this 
      processing.
ntg - (near trace gather) print some headers, compute additiona trace headers 
      (offsets, cmps).  Sort to shotx,groupx order, create and display fold 
      plot.  make and plot near trace gather.  Plot every 10th trace
cvs - (constant velocity stack)  plot gathers with nmo applied.  Make and plot 
      constant velocity stacks for range of velocities.  Compute v(t) stack.
vscan - (velocity scan)  compute and plot velocity semblance.


