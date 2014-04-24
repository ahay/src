/* Linear/parabolic radon operator in frequency domain
Note: I borrowed a lot from /system/seismic/radon+Mradon.c. 
The distinction:
	I am using FFTW because I am inexperienced in invoking kiss_fft. 
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  References: 
	1) Kostov C., 1990. "Toeplitz structure in Slant-Stack inversion": 
	SEG Extended Abstracts, 1647-1650.
	2) Sacchi, Mauricio D., and Milton Porsani. "Fast high resolution 
	parabolic Radon transform." Society of Exploration Geophysicists 
	69th Annual International Meeting, SPRO P. Vol. 1. No. 1. 1999.
*/

#include <rsf.h>
#include <time.h>
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "myradon2.h"
#include "ctoeplitz_reg.h"

static int np, nx;
static float w, dp, *p, *xx;

void myradon2_init(int np_, int nx_, float dp_, float *p_, float *xx_)
/*< initialization >*/
{
	np=np_;
	nx=nx_;
	dp=dp_;
	p=p_;
	xx=xx_;
}

void myradon2_set(float w_)
/*< set up frequency w >*/
{
	w=w_;
}

void myradon2_lop(bool adj, bool add, int nm, int nd, sf_complex *mm, sf_complex *dd)
/*< radon linear operator >*/
{
	int ix, ip;
	sf_complex sumc;

    	if (nm != np || nd != nx) sf_error("%s: mismatched data sizes",__FILE__);
	
	sf_cadjnull(adj, add, nm, nd, mm, dd);

	if(adj){// mm(p,w)=sum_{ix=0}^{nx} dd(xx[ix],w)*exp(i*w*p*xx[ix])
		for(ip=0; ip<np; ip++) // loop over slopes
		{
			sumc=sf_cmplx(0,0);
			for(ix=0; ix<nx; ix++) 
				sumc+=cexpf(I*w*p[ip]*xx[ix])*dd[ix];
			mm[ip]=sumc;
		}
	}else{// dd(xx,w)=sum_{ip=0}^{np} mm(p[ip],w)*exp(-i*w*p[ip]*xx)
		for(ix=0; ix<nx; ix++) 
		{
			sumc=sf_cmplx(0,0);
			for(ip=0; ip<np; ip++)
				sumc+=cexpf(-I*w*p[ip]*xx[ix])*mm[ip];
			dd[ix]=sumc;
		}
	}
}

void myradon2_inv(sf_complex *mm, sf_complex *mm_adj, float eps)
/*< fast Toeplitz matrix inversion for radon transform 
mm: model to be inverted
mm_adj: adjoint of radon
eps: regularization parameter
>*/
{
	int ip, ix;
	sf_complex *c, sumc;
	c=(sf_complex *)malloc(np*sizeof(sf_complex));
	
	for(ip=0; ip<np; ip++) 
	{
		sumc=sf_cmplx(0,0);
		for(ix=0; ix<nx; ix++)
			sumc+=cexpf(I*w*ip*dp*xx[ix]);
		c[ip]=sumc;
	}
	
	ctoeplitz_inv(np, eps, c, mm, mm_adj);
}
