#ifndef __SEAM_CUBIC__
#define __SEAM_CUBIC__

#include "cstd.h"

/** \page cubic Cubic Spline Interpolation and Adjoint Interpolation

Documentation in cubic.h.

SEAM 1D cubic spline interpolation and adjoint interpolation. Based on
TRIP Fortran code written by Ann Campbell as an undergraduate summer
project in 1991. See source code for original documentation.

The cubic interpolation routine implements a perfectly ordinary spline algorithm from Forsythe, Malcolm, and Moler's book. This package's main claim to fame is the adjoint interpolation implementation; I have no idea where else you can find such a thing.

WWS, July 08
*/

/** Returns necessary workspace for interpolation */
int cubic_getworksize(int nt);

/** Cubic spline interpolation.

@param[in] nin (int) number of input gridpoints;
@param[in] din (float) step for input grid;
@param[in] oin (float) coordinate of first input sample;
@param[in] vin (float *) input data array;
@param[in] nout (int) number of output gridpoints;
@param[in] dout (float) step for output grid;
@param[in] oout (float) coordinate of first output sample;
@param[out] vout (float *) output data array;
@param[in] iend (int) endpoint code:  
<ol>
<li> iend=1: linear ends, s(1)=s(n)=0 </li>
<li> iend=2: parabolic ends, s(1)=s(2),s(n)=s(n-1) </li>
<li> iend=3 cubic ends, s(1),s(n) are extrapolated </li>
</ol>
@param[in] work (float *) workspace (memory managed externally);
@param[in] wl (int) size of workspace.

Code translated from Fortran by f2c, so all args passed by pointer,
whether they are logically const or not.

The helper function cubic_getworksize returns the necessary workspace length.

\retval 0 normal return
\retval 11 incorrect end condition given (must lie between 1 and 3)
\retval 12 not enough input data points, must be > 3 
\retval 10 not enough workspace, needs 4*ni-3 words 

*/ 

int cubic_(float const *oin,  
	   float const *din,  
	   float const *vin,  
	   int const *nin, 
	   float const *oout, 
	   float const *dout, 
	   float *vout, 
	   int const *nout,
	   int const *iend, 
	   float *work, 
	   int const *wl);

/** Returns workspace size for adjoint interpolation. */

int cubicadj_getworksize(int nout, int nin);

/** Cubic spline adjoint interpoation.

@param[in] nin (int) number of input gridpoints;
@param[in] din (float) step for input grid;
@param[in] oin (float) coordinate of first input sample;
@param[in] vin (float *) input data array;
@param[in] nout (int) number of output gridpoints;
@param[in] dout (float) step for output grid;
@param[in] oout (float) coordinate of first output sample;
@param[out] vout (float *) output data array;
@param[in] iend (int) endpoint code:  
<ol>
<li> iend=1: linear ends, s(1)=s(n)=0 </li>
<li> iend=2: parabolic ends, s(1)=s(2),s(n)=s(n-1) </li>
<li> iend=3 cubic ends, s(1),s(n) are extrapolated </li>
</ol>
@param[in] work (float *) workspace (memory managed externally);
@param[in] wl (int) size of workspace.

Code translated from Fortran by f2c, so all args passed by pointer,
whether they are logically const or not.

The helper function cubicadg_getworksize returns the necessary workspace length.

\retval 0 normal return
\retval 11 incorrect end condition given (must lie between 1 and 3)
\retval 12 not enough input data points, must be > 3 
\retval 10 not enough workspace, needs nin+7*nout-8 words.
*/

int cubicadj_(float const *oin,  
	      float const *din,  
	      float const *vin,  
	      int const *nin,
	      float const *oout, 
	      float const *dout, 
	      float *vout, 
	      int const *nout, 
	      int const *iend, 
	      float *work, 
	      int const *wl);

#endif

