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

Reference: Ji, Jun. "An exact adjoint operation pair in time extrapolation 
and its application in least-squares reverse-time migration." Geophysics 
74.5 (2009): H27-H33.
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>

#if defined __INTEL_COMPILER && __INTEL_COMPILER_BUILD_DATE >= 20130607 && _OPENMP < 201107
#undef  _OPENMP
#define _OPENMP 201107
#endif

/* collapse is a feature from OpenMP 3 (2008) */
#if _OPENMP < 200805
    #define collapse(x) 
#endif

#endif

#include "rtm2d.h"

static int n0, nzpad, nxpad, nb, nz, nx, nt, nm, nd;
static float c0, c11, c21, c12, c22;
static float *bndr, *mod, *dat;
static float **u0, **u1, **u2, **vv, **ptr=NULL;

void step_forward(float **u0, float **u1, float **u2, float **vv, bool adj)
{
	int i1, i2;
	float u;

	if(adj){
#ifdef _OPENMP
#if _OPENMP >= 201107
#pragma omp parallel for collapse(2) default(none)	\
    private(i2,i1,u)		    			\
    shared(nzpad,nxpad,u1,vv,u0,u2,c0,c11,c12,c21,c22)
#else
#pragma omp parallel for default(none)	\
    private(i2,i1,u)		    			\
    shared(nzpad,nxpad,u1,vv,u0,u2,c0,c11,c12,c21,c22)
#endif
#endif
		for (i2=0; i2<nxpad; i2++) 
		for (i1=0; i1<nzpad; i1++) 
		{
			u 		= vv[i2][i1  ]*c0 *u1[i2][i1  ];
			if(i1 >= 1)   u+= vv[i2][i1-1]*c11*u1[i2][i1-1];
			if(i1 >= 2)   u+= vv[i2][i1-2]*c12*u1[i2][i1-2];
			if(i1 < nzpad-1) u+= vv[i2][i1+1]*c11*u1[i2][i1+1];
			if(i1 < nzpad-2) u+= vv[i2][i1+2]*c12*u1[i2][i1+2];
			if(i2 >= 1)   u+= vv[i2-1][i1]*c21*u1[i2-1][i1];
			if(i2 >= 2)   u+= vv[i2-2][i1]*c22*u1[i2-2][i1];
			if(i2 < nxpad-1) u+= vv[i2+1][i1]*c21*u1[i2+1][i1];
			if(i2 < nxpad-2) u+= vv[i2+2][i1]*c22*u1[i2+2][i1];
			u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+u;
		}
	}else{
#ifdef _OPENMP
#if _OPENMP >= 201107
#pragma omp parallel for collapse(2) default(none)	\
    private(i2,i1,u)		    			\
    shared(nzpad,nxpad,u1,vv,u0,u2,c0,c11,c12,c21,c22)
#else
#pragma omp parallel for default(none)	\
    private(i2,i1,u)		    			\
    shared(nzpad,nxpad,u1,vv,u0,u2,c0,c11,c12,c21,c22)
#endif
#endif
		for (i2=0; i2<nxpad; i2++) 
		for (i1=0; i1<nzpad; i1++) 
		{
			u 		= c0*u1[i2][i1];
			if(i1 >= 1)   u+= c11*u1[i2][i1-1];
			if(i1 >= 2)   u+= c12*u1[i2][i1-2];
			if(i1 < nzpad-1) u+= c11*u1[i2][i1+1];
			if(i1 < nzpad-2) u+= c12*u1[i2][i1+2];
			if(i2 >= 1)   u+= c21*u1[i2-1][i1];
			if(i2 >= 2)   u+= c22*u1[i2-2][i1];
			if(i2 < nxpad-1) u+= c21*u1[i2+1][i1];
			if(i2 < nxpad-2) u+= c22*u1[i2+2][i1];
			u2[i2][i1]=2*u1[i2][i1]-u0[i2][i1]+vv[i2][i1]*u;
		}
	}
}

void apply_sponge(float **p0, float *bndr)
/* apply sponge absorbing boundary condition */
{
	int ix,iz,ib,ibx,ibz;
	float w;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(p0,bndr,nzpad,nxpad,nb)
#endif
    for(ib=0; ib<nb; ib++) {
	w = bndr[ib];

	ibz = nzpad-ib-1;
	for(ix=0; ix<nxpad; ix++) {
	    p0[ix][ib ] *= w; /*    top sponge */
	    p0[ix][ibz] *= w; /* bottom sponge */
	}

	ibx = nxpad-ib-1;
	for(iz=0; iz<nzpad; iz++) {
	    p0[ib ][iz] *= w; /*   left sponge */
	    p0[ibx][iz] *= w; /*  right sponge */
	}
    }
}


void rtm2d_init(float dz_, float dx_, float dt_, int n0_, int nz_, 
int nx_, int nb_, int nt_, float **v0, float *mod_, float *dat_)
/*< allocate variables and initialize parameters >*/
{  
   	/* initialize OpenMP support */ 
#ifdef _OPENMP
   	omp_init();
#endif
	float t;
	int ib,i1,i2;
	t = 1.0/(dz_*dz_);
	c11 = 4.0*t/3.0;
	c12= -t/12.0;

	t = 1.0/(dx_*dx_);
	c21 = 4.0*t/3.0;
	c22= -t/12.0;
	c0=-2.0*(c11+c12+c21+c22);

	n0=n0_;
	nb=nb_;
	nz=nz_;
	nx=nx_;
	nt=nt_;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nm=nzpad*nxpad;
	nd=nt*nx;

    	/* allocate temporary arrays */
    	bndr=sf_floatalloc(nb);
    	u0=sf_floatalloc2(nzpad,nxpad);
    	u1=sf_floatalloc2(nzpad,nxpad);
    	u2=sf_floatalloc2(nzpad,nxpad);
    	vv=sf_floatalloc2(nzpad,nxpad);
	/* initialized sponge ABC coefficients */
	for(ib=0;ib<nb;ib++){
		t=0.015*(nb-1-ib);
		bndr[ib]=expf(-t*t);
	}
	mod=mod_;
	dat=dat_;
    	for (i2=0; i2<nx; i2++) 
	for (i1=0; i1<nz; i1++) 
	{
	    t=v0[i2][i1]*dt_;
	    vv[i2+nb][i1+nb] = t*t;
	}	
    	for (i2=0; i2<nxpad; i2++)
	for (i1=0; i1<nb; i1++) 
	{
	    vv[i2][   i1  ] =vv[i2][   nb  ];
	    vv[i2][nzpad-i1-1] =vv[i2][nzpad-nb-1];
	}
    	for (i2=0; i2<nb; i2++)
	for (i1=0; i1<nzpad; i1++) 
	{
	    vv[   i2  ][i1] =vv[   nb  ][i1];
	    vv[nxpad-i2-1][i1] =vv[nxpad-nb-1][i1];
	}
}

void rtm2d_close()
/*< free allocated variables >*/
{
	free(bndr);
	free(*u0); free(u0);
	free(*u1); free(u1);
	free(*u2); free(u2);
	free(*vv); free(vv);
}

void rtm2d_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< rtm2d linear operator >*/
{
	int i1,i2,it;
	sf_adjnull (adj, add, nm, nd, mod, dat);
 
	memset(u0[0], 0, nzpad*nxpad*sizeof(float));
	memset(u1[0], 0, nzpad*nxpad*sizeof(float));
	memset(u2[0], 0, nzpad*nxpad*sizeof(float));

    	if(adj){/* migration */
	    	for (it=nt-1; it >-1; it--) {
			sf_warning("%d;",it);

			apply_sponge(u0, bndr); 
			apply_sponge(u1, bndr); 
			step_forward(u0, u1, u2, vv, true);
			ptr=u0; u0=u1; u1=u2; u2=ptr;

#ifdef _OPENMP
#pragma omp parallel for	\
    private(i2)		    	\
    shared(it,u1,dat,n0,nb)
#endif
			for (i2=0; i2<nx; i2++)	u1[i2+nb][n0+nb] += dat[i2*nt+it];
	    	}
		/* output image (mod is image) */
		for(i2=0; i2<nxpad; i2++)
		for(i1=0; i1<nzpad; i1++)
		{
			mod[i1+nzpad*i2]+=u1[i2][i1];
		}
    	}else{ /* modeling */
		for(i2=0; i2<nxpad; i2++)
		for(i1=0; i1<nzpad; i1++)
		{
			u1[i2][i1]+=mod[i1+nzpad*i2];
		}
    		for (it=0; it <nt; it++) {
			sf_warning("%d;",it);
	  
		/* record data */
#ifdef _OPENMP
#pragma omp parallel for	\
    private(i2)		    	\
    shared(it,u1,dat,n0,nb)
#endif
			for (i2=0; i2<nx; i2++) dat[i2*nt+it]+=u1[i2+nb][n0+nb];
	
			step_forward(u0, u1, u2, vv, false);
			apply_sponge(u1, bndr);
			apply_sponge(u2, bndr);
			ptr=u0; u0=u1; u1=u2; u2=ptr;
    		}
    	}
}

