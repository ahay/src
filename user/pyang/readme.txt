Author: Pengliang Yang, Xi'an Jiatong Universtiy
Email: 	ypl.2100@gmail.com

=================================================================
Preamble: This readme files is devoted to explain my programs. 
Some of them has been tested, and the rest are under construction.
Be careful when you try to use them!
=================================================================

Under this directory, I implemented the algorithms:

1) POCS (projection onto convex sets), FFTW required
Main:		Mpocs3d.c, Mfpocs3d.c
Depends on:	ft3d.c
Test file: 	/book/xjtu/test/fpocs3d/SConstruct
		/book/xjtu/test/fpocs2d/SConstruct
Note: fpocs is a two-step version of POCS. 

2) DLCT (discrete linear chirp transform), FFTW required
Main:		Mdlct.c, Mdlct2.c
Depends on: 	dlct.c
Test file: 	/book/xjtu/test/dlct/SConstruct
Note: To make the adjoint of DLCT same as inverse, I normalized 
	the forward and inverse DLCT with a factor. 

3) 2D and 3D FD for forward modelling
Main:		MTestfd2d.c, MTestfd3d.c
Test file:	/book/xjtu/test/Testfd3d/SConstruct
		/book/xjtu/test/Testfd2d/SConstruct

4) RTM and LSRTM (2-D zero-offset least squares RTM)
Main: 		Mrtm2d.c Mlsrtm2d.c
Depends on: 	rtm2d.c
Test file:	/book/xjtu/rtm2d/hyper/SConstruct
		/book/xjtu/rtm2d/marmousi/SConstruct
		/book/xjtu/rtm2d/sigsbee/SConstruct
Note: rtm2d.c is coded following the linear operator standard in 
	Madagascar:	oper(adj, add, nm, nd, mod, dat)


5) Prestack RTM using GPU
Main: 		staggered_fdcoeff.m, MTesteb.c, Mgpurtm.c
Depends on:	cuda_kernels.cu
Test file: 	/book/xjtu/test/Testeb/SConstruct
Note: 	(a)staggered_fdcoeff.m is a matlab script to find the finite 
	difference coefficients with order-NJ(NJ=2N);
	(b) MTesteb.c is a file to test the validity of the proposed
	effective boundary saving strategy! 
	(c) Most of the detail explaination for GPU-based RTM can be
	found in the codes.

6) Seislet-based POCS, IST and MCA algorithm (for 2D validation)
Main:		Mpocsseislet.c, Mistseislet.c, Mmcaseislet
Test file:	/book/xjtu/mcaseislet/deblend/SConstruct
		/book/xjtu/mcaseislet/interp/SConstruct
		/book/xjtu/mcaseislet/sep1/SConstruct
		/book/xjtu/mcaseislet/sep2/SConstruct

7) 2-D forward modeling to generate shot records
Main: 		Mmodeling2d.c
Test file: 	/book/xjtu/test/modeling2d/SConstruct


8) Elastic and anisotropic modeling
Main: 		MTestelastic2d.c, MTestaniso.c
Test file:	/book/xjtu/test/testelastic2d/SConstruct
		/book/xjtu/test/testaniso/SConstruct

9) Radon transforms

Main: 		Mmyradon1.c Mmyradon2.c
Test file: 	/book/xjtu/test/myradon2/SConstruct


===================================================================
The following codes are under construction. Be careful!
===================================================================
10) 3D coherence calculation
Main: 		Mcohn.c
Depends on: 	svd.c

11) MWNI (minimum weighted norm interpolation), FFTW requred
Main:		Mmwni2d.c Mmwni3d.c
Test file: 	/book/xjtu/test/mwni2d/SConstruct
Note: I use conjugate gradient algorithm here. Although the testing
seems nice, I found the residual of my implementation not converged
well. Be careful! It is under modification!

==================================================================
Miscellaneous:
maskft3d.c: a mask applied on 3d fft


