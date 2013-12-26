Author: Pengliang Yang, Xi'an Jiatong Universtiy

=================================================================
Preamble: This readme files is devoted to explain my programs. 
Some of them has been tested, and the rest are under construction.
Be careful when you try to use them!
=================================================================

Under this directory, I implemented the algorithms:

1) POCS (projection onto convex sets), FFTW required

Main:		Mpocs3d.c, Mfpocs3d.c

Test file: SConstruct.fpocs3d (put it in another directory 
and rename it as SConstruct, run 'scons view')

2) DLCT (discrete linear chirp transform), FFTW required

Main:		Mdlct.c, 
Depends on: 	dlct.c

Test file: SConstruct.dlct (put it in another directory 
and rename it as SConstruct, run 'scons view')

Note: To make the adjoint of DLCT same as inverse, I normalized 
	the forward and inverse DLCT with a factor. 

3) 2D and 3D FD for forward modelling

Main:		MTestfd2d.c, MTestfd3d.c
Depends on: 	fd3dutil.c 
(fd3dutil.c is modified from fdutil.c in Madagascar.)

Test file:	SConstruct.Testfd2d SConstruct.Testfd3d

===================================================================
The following codes are under construction. Be careful!
===================================================================

4) MWNI (minimum weighted norm interpolation), FFTW requred

Main:		Mmwni2d.c Mmwni3d.c

Test file: SConstruct.mwni2d

Note: I use conjugate gradient algorithm here. Although the testing
seems nice, I found the residual of my implementation not converged
well. Be careful! It is under modification!

5) LSRTM (least squares reverse time migration), under construction

Main: 		Mlsrtm2d.c

Note:  Be careful! It is under construction!

==================================================================
Miscellaneous:
_segy.h: segy format definations
triangle.c,ctriange.c: taking from others for future usage.
ft3d.c(fftw required): 3D FFT is coded as a linear operator in 
Madagascar linear operator standard.
	oper(adj, add, nm, nd, mod, dat);
maskft3d.c: a mask applied on 3d fft to constrain the bandlimitedness.


