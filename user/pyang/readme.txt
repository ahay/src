Author: Pengliang Yang, Xi'an Jiatong Universtiy, UT Austin
Email: 	ypl.2100@gmail.com

=================================================================
Preamble: This readme files is committed to explain my programs. 
Some of them have been tested. The rest are under construction.
Be careful when you try to use them! Please feel free to contact
me if you find errors or have suggestions!
=================================================================

In this directory, I implemented the algorithms:

1) POCS (projection onto convex sets), FFTW required
Main:		Mpocs3d.c, Mfpocs2d.c, Mfpocs3d.c, Mpocs.c,
Depends on:	fftn.c
Test file: 	/book/xjtu/test/fpocs3d/SConstruct
		/book/xjtu/test/fpocs2d/SConstruct
		/book/xjtu/test/pocs5d/SConstruct
Note: fpocs is a two-step version of POCS. You are able to test 
    sfpocs3d and sffpocs3d by changing the name of the
    program in /book/xjtu/test/fpocs3d/SConstruct. 

    POCS implemented in frequency domain (see sfpocs) data gives 
    much faster speed. It saves the storage and computational cost 
    of FFT. sfpocs can handle any dimensional data interpolation.
    Thanks ot the data generation effort done by Sergey in the sript
    of /book/xjtu/test/pocs5d/SConstruct as the validation.

For more information on FPOCS, check the paper:
 Yang Pengliang, Gao Jinghuai, Chen Wenchao, "On analysis-based 
 two-step interpolation methods for randomly sampled seismic data"
 Computers & Geosciences 51 449-461 2013


2) DLCT (discrete linear chirp transform), FFTW required
Main:		Mdlct.c, Mdlct2.c
Depends on: 	dlct.c
Test file: 	/book/xjtu/test/dlct/SConstruct
		/book/user/pyang/test_dlct1.m
		/book/user/pyang/test_dlct2.m
Note: To make the adjoint of DLCT same as inverse, I normalized 
	the forward and inverse DLCT with a factor. 
The matlab scripts can run directly within MATLAB software environment.

3) 2-D and 3-D FD for forward modelling
Main:		MTestfd2d.c, MTestfd3d.c
Test file:	/book/xjtu/test/Testfd3d/SConstruct
		/book/xjtu/test/Testfd2d/SConstruct

4) 2-D Least-squares reverse time migration (LSRTM)
Main:		Mrtm2d.c, Mlsrtm2d.c, Mlsprtm2d.c
Depends on:	rtm2d.c, prtm2d.c
Test file:	/book/xjtu/test/zortm2d/marmousi/SConstruct
		/book/xjtu/test/zortm2d/sigsbee/SConstruct
		(zero-offset RTM)
Test file: 	/book/xjtu/test/lsprtm2d/SConstruct (under construction)
		(prestack RTM)

5) Prestack RTM using GPU with staggered grid
Main: 		staggered_fdcoeff.m, MTesteb.c, Mgpurtm.c
Depends on:	cuda_kernels.cu
Test file: 	/book/xjtu/gpurtm/marmousi/SConstruct
		/book/xjtu/gpurtm/sigsbee/SConstruct
		/book/xjtu/test/Testeb/SConstruct
Note: 	(a)staggered_fdcoeff.m is a matlab script to find the finite 
	difference coefficients with order-NJ (NJ=2N);
	(b) MTesteb.c is a file to test the validity of the proposed
	effective boundary saving strategy! 
	(c) Most of the detail explaination for GPU-based RTM can be
	found in the codes.
For more information, check the paper:
 Yang, Pengliang, Jinghuai Gao, and Baoli Wang. "RTM using effective 
 boundary saving: A staggered grid GPU implementation." Computers & 
 Geosciences (2014).


6) Seislet-based POCS, IST and MCA algorithm (for 2D validation)
Main:		Mpocsseislet.c, Mistseislet.c, Mmcaseislet
Test file:	/book/xjtu/mcaseislet/deblend/SConstruct
		/book/xjtu/mcaseislet/interp/SConstruct
		/book/xjtu/mcaseislet/sep1/SConstruct
		/book/xjtu/mcaseislet/sep2/SConstruct
		/book/xjtu/test/istpocsseislet/SConstruct

7) 2-D forward modeling to generate shot records
Main: 		Mmodeling2d.c
Test file: 	/book/xjtu/test/modeling2d/SConstruct


8) Elastic and anisotropic modeling
Main: 		MTestelastic2d.c, MTestaniso.c
Test file:	/book/xjtu/test/testelastic2d/SConstruct
		/book/xjtu/test/testaniso/SConstruct

9) Radon transforms
Main: 		Mmyradon1.c Mmyradon2.c
Depends on:	ctoeplitz_reg.c, myradon2.c
Test file: 	/book/xjtu/test/myradon2/SConstruct
NB: 	Myradon1.c is not tested yet!

10) Angle gather (ADCIG) extraction using Poynting vector
Main: 		Mrtmadcig.c
Test file: 	/book/xjtu/primer/rtmadcig/SConstruct
NB: ADCIG computation is much more expensive than RTM imaging.
An MPI version of this program is in preparation!

11) 2D GPU-based full waveform inversion (FWI)
Main:		Mgenshots.cu, Mfbrec.cu,Mgpufwi.cu
Test file:	/book/xjtu/gpufwi/fbrec/SConstruct
		/book/xjtu/gpufwi/syntest/SConstruct
		/book/xjtu/gpufwi/marmtest/SConstruct
Note: Mgenshots.cu is used to generate shots by forward modeling using
the exact velocity model. We can use a starting model to do FWI by
invoking Mgpufwi.cu. We are using boundary saving strategy in FWI. To 
demonstrate that the modeled wavefield can be precisely reconstructed,
we design the code Mfbrec.cu befere going to FWI. It is important to
note that the top boundary is free surface boundary condition (no ABC 
applied here!). 

12) 3D coherence calculation
Main: 		Mcohn.c
Test file:	/book/xjtu/test/coherence/SConstruct


13) MWNI (minimum weighted norm interpolation), FFTW requred
Main:		Mmwni2d.c 
Test file: 	/book/xjtu/test/mwni2d/SConstruct
Note: I use conjugate gradient algorithm here. Although the testing
seems nice, I found the residual of my implementation not converged
well. 

14) 3D FD using GPU
Main: 		Mgpufd3d.cu, Mgpufbrec3d.cu
Test file:	/book/xjtu/test/gpufd3d/SConstruct
NB: Mgpufbrec3d.cu is performing backward reconstruction for the forward 
modeled wavefield in 3D with GPU. It is implemented 2nd order FD, and 
prepared for 3D GPU-based RTM

===================================================================
I try my best to make my code self-contained. I believe it brings 
convenience and readability, because for readers much effort will
be saved on understanding how to invoke complicated functions which may 
not be written by the author. I strongly discourage that kind of style!
===================================================================



