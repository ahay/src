/* 2-D zero-offset reverse-time migration linear operator
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
#endif

#include "rtm2d.h"

static int n0, nzpad, nxpad, nbl, nbr, nbt, nbb, nz, nx, nt, nm, nd;
static float c0, c11, c21, c12, c22;
static float *mod, *dat;
static float **u0, **u1, **vv, **ptr=NULL;

void step_forward(float **u0, float **u1, float **vv, bool adj)
{
  int i1, i2;

  if(adj){
#ifdef _OPENMP
#pragma omp parallel for default(none)			\
  private(i2,i1)					\
  shared(nzpad,nxpad,u1,vv,u0,c0,c11,c12,c21,c22)
#endif
    for (i2=2; i2<nxpad-2; i2++) 
      for (i1=2; i1<nzpad-2; i1++) 
	{
	  u0[i2][i1]=2.*u1[i2][i1]-u0[i2][i1]+
	    c0 *vv[i2][i1  ]*u1[i2][i1  ]+
	    c11*(vv[i2][i1-1]*u1[i2][i1-1]+vv[i2][i1+1]*u1[i2][i1+1])+
	    c12*(vv[i2][i1-2]*u1[i2][i1-2]+vv[i2][i1+2]*u1[i2][i1+2])+
	    c21*(vv[i2-1][i1]*u1[i2-1][i1]+vv[i2+1][i1]*u1[i2+1][i1])+
	    c22*(vv[i2-2][i1]*u1[i2-2][i1]+vv[i2+2][i1]*u1[i2+2][i1]);
	}
  }else{
#ifdef _OPENMP
#pragma omp parallel for default(none)			\
  private(i2,i1)					\
  shared(nzpad,nxpad,u1,vv,u0,c0,c11,c12,c21,c22)
#endif
    for (i2=2; i2<nxpad-2; i2++) 
      for (i1=2; i1<nzpad-2; i1++) 
	{
	  u0[i2][i1]=2.*u1[i2][i1]-u0[i2][i1]+
	    vv[i2][i1]*(c0*u1[i2][i1]+
			c11*(u1[i2][i1-1]+u1[i2][i1+1])+
			c12*(u1[i2][i1-2]+u1[i2][i1+2])+
			c21*(u1[i2-1][i1]+u1[i2+1][i1])+
			c22*(u1[i2-2][i1]+u1[i2+2][i1]));
	}
  }
}


void rtm2d_init(float dz_, float dx_, float dt_, int n0_, int nz_, 
		int nx_, int nt_, float **v0, float *mod_, float *dat_)
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
  nz=nz_;
  nx=nx_;
  nt=nt_;
  nbl=nx/2;//number of left boundary layers
  nbr=nx-nbl;//number of right boundary layers
  nbt=nz/2;//number of top boundary layers
  nbb=nz-nbt;//number of bottom boundary layers
  nzpad=nbt+nz+nbb;
  nxpad=nbl+nx+nbr;
  nm=nzpad*nxpad;
  nd=nt*nx;

  /* allocate temporary arrays */
  u0=sf_floatalloc2(nzpad,nxpad);
  u1=sf_floatalloc2(nzpad,nxpad);
  vv=sf_floatalloc2(nzpad,nxpad);
  mod=mod_;
  dat=dat_;
  for (i2=0; i2<nx; i2++) 
    for (i1=0; i1<nz; i1++) 
      {
	t=v0[i2][i1]*dt_;
	vv[i2+nbl][i1+nbt] = t*t;
      }	
  for (i2=0; i2<nxpad; i2++)
    {
      for(i1=0; i1<nbt; i1++)//top boundary
	vv[i2][   i1  ] =vv[i2][   nbt  ];
      for(i1=0; i1<nbb; i1++)//bottom boundary
	vv[i2][nz+nbt+i1] =vv[i2][nz+nbt-1];
    }
  for(i1=0; i1<nzpad; i1++)
    {   
      for(i2=0; i2<nbl; i2++)//left boundary
	vv[   i2  ][i1] =vv[   nbl  ][i1];
      for(i2=0; i2<nbr; i2++)//right boundary
	vv[nx+nbl+i2][i1] =vv[nx+nbl-1][i1];
    }
}

void rtm2d_close()
/*< free allocated variables >*/
{
  free(*u0); free(u0);
  free(*u1); free(u1);
  free(*vv); free(vv);
}

void rtm2d_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< rtm2d linear operator >*/
{
  int i1,i2,it;

  sf_adjnull (adj, add, nm, nd, mod, dat); 
  memset(u0[0], 0, nzpad*nxpad*sizeof(float));
  memset(u1[0], 0, nzpad*nxpad*sizeof(float));

  if(adj){// migration
    for (it=nt-1; it >-1; it--) {
      sf_warning("%d;",it);

      step_forward(u0, u1, vv, true);
      ptr=u0; u0=u1; u1=ptr;

#ifdef _OPENMP
#pragma omp parallel for			\
  private(i2)					\
  shared(it,u1,dat,n0,nbl,nbt)
#endif
      for (i2=0; i2<nx; i2++)	u1[i2+nbl][n0+nbt] += dat[i2*nt+it];

    }
    /* output image (mod is image, i.e. truncation of u1) */
    for(i2=0; i2<nx; i2++)
      for(i1=0; i1<nz; i1++)
	mod[i1+nz*i2]+=u1[i2+nbl][i1+nbt];
  }else{ // modeling
    for(i2=0; i2<nx; i2++)
      for(i1=0; i1<nz; i1++)
	u1[i2+nbl][i1+nbt]+=mod[i1+nz*i2];// zero-paded u1

    for (it=0; it <nt; it++) {
      sf_warning("%d;",it);	
	  
      /* record data */
#ifdef _OPENMP
#pragma omp parallel for			\
  private(i2)					\
  shared(it,u1,dat,n0,nbl,nbt)
#endif
      for (i2=0; i2<nx; i2++) dat[i2*nt+it]+=u1[i2+nbl][n0+nbt];

      step_forward(u0, u1, vv, false);
      ptr=u0; u0=u1; u1=ptr;
    }
  }
}

