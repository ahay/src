/* 2-D zero-offset reverse-time migration linear operator
Note: 	1) Migration should be the adjoint of modeling. If you pass the 
    dot-product test, you can use least-squares method with this 
    rtm operator! 
	2) Here, the sponge absorbing boundary condition is applied!
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
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "rtm2d.h"

static int n0, n1, n2, nb, nt, nm, nd;
static float c0, c11, c21, c12, c22;
static float *bndr, *mod, *dat;
static float **u0, **u1, **u2, **vv, **ptr=NULL;

void step_forward(float **u0, float **u1, float **u2, float **vv)
{
	int i1, i2;
/*
// This part is designed for the adjointness, every points in the 
// grid has been included, but slow running.
	for (i2=0; i2<n2; i2++) 
	for (i1=0; i1<n1; i1++) 
	{
		float u = c0*u1[i2][i1];
		if(i1 >= 1) u += c11*u1[i2][i1-1];
		if(i1 >= 2) u += c12*u1[i2][i1-2];
		if(i1 < n1-1) u += c11*u1[i2][i1+1];
		if(i1 < n1-2) u += c12*u1[i2][i1+2];
		if(i2 >= 1) u += c21*u1[i2-1][i1];
		if(i2 >= 2) u += c22*u1[i2-2][i1];
		if(i2 < n2-1) u += c21*u1[i2+1][i1];
		if(i2 < n2-2) u += c22*u1[i2+2][i1];
		u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+vv[i2][i1]*u;

		//float 	u = vv[i2][i1  ]*c0 *u1[i2][i1  ];
		//if(i1 >= 1)   u+= vv[i2][i1-1]*c11*u1[i2][i1-1];
		//if(i1 >= 2)   u+= vv[i2][i1-2]*c12*u1[i2][i1-2];
		//if(i1 < n1-1) u+= vv[i2][i1+1]*c11*u1[i2][i1+1];
		//if(i1 < n1-2) u+= vv[i2][i1+2]*c12*u1[i2][i1+2];
		//if(i2 >= 1)   u+= vv[i2-1][i1]*c21*u1[i2-1][i1];
		//if(i2 >= 2)   u+= vv[i2-2][i1]*c22*u1[i2-2][i1];
		//if(i2 < n2-1) u+= vv[i2+1][i1]*c21*u1[i2+1][i1];
		//if(i2 < n2-2) u+= vv[i2+2][i1]*c22*u1[i2+2][i1];
		//u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+u;
	}
*/

// This part is devoted to efficiency reasons, if you do not use 
// absorbing boundary, it is hard to be adjoint. If you use ABC, 
// you can have approximate adjoint.
#ifdef _OPENMP
#pragma omp parallel for		\
    private(i2,i1)		    	\
    shared(n1,n2,u1,u0,u2,c0,c11,c12,c21,c22)
#endif
	for(i2=2; i2<n2-2; i2++)
	for(i1=2; i1<n1-2; i1++)
	{
		u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+vv[i2][i1]*(	c0*u1[i2][i1]
				+c11*(u1[i2][i1-1]+u1[i2][i1+1])
				+c12*(u1[i2][i1-2]+u1[i2][i1+2])
				+c21*(u1[i2-1][i1]+u1[i2+1][i1])
				+c22*(u1[i2-2][i1]+u1[i2+2][i1]));
	}
}

void apply_sponge(float**p0, float **p1, float *bndr)
/* apply sponge absorbing boundary condition */
{
	int ix,iz;

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,p0,p1)
#endif
	for(ix=0; ix<n2; ix++)
	{
		for(iz=0;iz<nb;iz++){	// top ABC			
			p0[ix][iz]=bndr[iz]*p0[ix][iz];
			p1[ix][iz]=bndr[iz]*p1[ix][iz];
		}
		for(iz=n1-nb;iz<n1;iz++){// bottom ABC			
			p0[ix][iz]=bndr[n1-iz-1]*p0[ix][iz];
			p1[ix][iz]=bndr[n1-iz-1]*p1[ix][iz];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,p0,p1)
#endif
	for(iz=0; iz<n1; iz++)
	{
		for(ix=0;ix<nb;ix++){	// left ABC			
			p0[ix][iz]=bndr[ix]*p0[ix][iz];
			p1[ix][iz]=bndr[ix]*p1[ix][iz];
		}	
		for(ix=n2-nb;ix<n2;ix++){// right ABC			
			p0[ix][iz]=bndr[n2-ix-1]*p0[ix][iz];
			p1[ix][iz]=bndr[n2-ix-1]*p1[ix][iz];
		}	
	}
}

void rtm2d_init(float dz_, float dx_, float dt_, int n0_, int n1_, 
int n2_, int nb_, int nt_, float **vv_, float *mod_, float *dat_)
/*< allocate variables and initialize parameters >*/
{  
   	/* initialize OpenMP support */ 
#ifdef _OPENMP
   	omp_init();
#endif
	float t;
	t = 1.0/(dz_*dz_);
	c11 = 4.0*t/3.0;
	c12= -t/12.0;

	t = 1.0/(dx_*dx_);
	c21 = 4.0*t/3.0;
	c22= -t/12.0;
	c0=-2.0*(c11+c12+c21+c22);

	n0=n0_;
	n1=n1_;
	n2=n2_;
	nb=nb_;
	nt=nt_;
	nm=n1*n2;
	nd=nt*n2;

    	/* allocate temporary arrays */
    	bndr=sf_floatalloc(nb);
    	u0=sf_floatalloc2(n1,n2);
    	u1=sf_floatalloc2(n1,n2);
    	u2=sf_floatalloc2(n1,n2);
	/* initialized sponge ABC coefficients */
	for(int ib=0;ib<nb;ib++){
		t=(nb-ib)/(sqrt(2.0)*4.0);
		bndr[ib]=expf(-t*t);
		//t=expf(-t*t);
		//bndr[ib]=powf(t,10.0);
	}
	vv=vv_;
	mod=mod_;
	dat=dat_;
    	float dt2 = dt_*dt_;
    	for (int i2=0; i2<n2; i2++) 
	for (int i1=0; i1<n1; i1++) 
	{
	    vv[i2][i1] *= vv[i2][i1]*dt2;
	}
}

void rtm2d_close()
/*< free allocated variables >*/
{
	free(bndr);
	free(*u0); free(u0);
	free(*u1); free(u1);
	free(*u2); free(u2);
}

void rtm2d_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< rtm2d linear operator >*/
{
	int i1,i2,it;
	sf_adjnull (adj, add, nm, nd, mod, dat);
 
	memset(u0[0], 0, n1*n2*sizeof(float));
	memset(u1[0], 0, n1*n2*sizeof(float));
	memset(u2[0], 0, n1*n2*sizeof(float));

    	if(adj){// migration
	    	for (it=nt-1; it >-1; it--) {
			sf_warning("%d;",it);

			apply_sponge(u0, u1, bndr); 
			step_forward(u0,u1,u2,vv);
			ptr=u0; u0=u1; u1=u2; u2=ptr;

#ifdef _OPENMP
#pragma omp parallel for	\
    private(i2)		    	\
    shared(it,u1,dat,n0)
#endif
			for (i2=0; i2<n2; i2++)	u1[i2][n0] += dat[i2*nt+it];
	    	}
		/* output image (mod is image) */
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
		{
			mod[i1+n1*i2]=u1[i2][i1];
		}
    	}else{ // modeling
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
		{
			u1[i2][i1]=mod[i1+n1*i2];
		}
    		for (it=0; it <nt; it++) {
			sf_warning("%d;",it);
	  
		/* record data */
#ifdef _OPENMP
#pragma omp parallel for	\
    private(i2)		    	\
    shared(it,u1,dat,n0)
#endif
			for (i2=0; i2<n2; i2++) dat[i2*nt+it]=u1[i2][n0];
	
			step_forward(u0,u1,u2,vv);
			apply_sponge(u1, u2, bndr);
			ptr=u0; u0=u1; u1=u2; u2=ptr;
    		}
    	}
}

