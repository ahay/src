Author: Pengliang Yang, Xi'an Jiatong Universtiy

=================================================================
Preamble: This readme files is devoted to explain my programs. 
Some of them has been tested, and the rest are under construction.
Be careful when you try to use them!
=================================================================

Under this directory, I implemented the algorithms:

1) POCS (projection onto convex sets), FFTW required

Main:		Mpocs3d.c, Mfpocs3d.c

Test file: 	/book/pyang/test/fpocs3d/SConstruct
		/book/pyang/test/fpocs2d/SConstruct

Note: fpocs is a two-step version of POCS. 

2) DLCT (discrete linear chirp transform), FFTW required

Main:		Mdlct.c, 
Depends on: 	dlct.c

Test file: 	/book/pyang/test/dlct/SConstruct

Note: To make the adjoint of DLCT same as inverse, I normalized 
	the forward and inverse DLCT with a factor. 

3) 2D and 3D FD for forward modelling

Main:		MTestfd2d.c, MTestfd3d.c
Depends on: 	fd3dutil.c 
(fd3dutil.c is modified from fdutil.c in Madagascar.)

Test file:	/book/pyang/test/Testfd3d/SConstruct
		/book/pyang/test/Testfd2d/SConstruct

4) RTM and LSRTM (2-D zero-offset least squares reverse time migration)

Main: 		Mrtm2d.c Mlsrtm2d.c
Depends on: 	rtm2d.c

Test file:	/book/pyang/test/lsrtm2d/hyper/SConstruct
		/book/pyang/test/lsrtm2d/sigsbee/SConstruct

Note: rtm2d.c is coded following the linear operator standard in 
	Madagascar:	oper(adj, add, nm, nd, mod, dat)

===================================================================
The following codes are under construction. Be careful!
===================================================================

5) MWNI (minimum weighted norm interpolation), FFTW requred

Main:		Mmwni2d.c Mmwni3d.c

Test file: 	/book/pyang/test/mwni2d/SConstruct

Note: I use conjugate gradient algorithm here. Although the testing
seems nice, I found the residual of my implementation not converged
well. Be careful! It is under modification!



==================================================================
Miscellaneous:
_segy.h: segy format definations
triangle.c,ctriange.c: taking from others for future usage.
ft3d.c(fftw required): 3D FFT is coded as a linear operator in 
Madagascar linear operator standard.
	oper(adj, add, nm, nd, mod, dat);
maskft3d.c: a mask applied on 3d fft to constrain the bandlimitedness.


