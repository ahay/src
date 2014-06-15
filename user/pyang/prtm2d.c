/* 2-D prestack LSRTM linear operator using wavefield reconstruction method
NB: Sponge ABC is applied!
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

#include "prtm2d.h"

static bool csdgather, verb;
static int sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz;
static int distx, distz, *sxz, *gxz;
static float *wlt, *bndr;

static int nzpad, nxpad, nb, nz, nx, nt, ns, ng, nm, nd;
static float c0, c11, c21, c12, c22;
static float *rwbndr, *mod, *dat;
static float **sp0, **sp1, **gp0, **gp1, **vv, **ptr=NULL;


void bndr_rw(bool read, float *rwbndr, float **p)
/*< read/write effective boundary >*/
{
	int i1, i2;

	if(read){// read bndr into wavefield	
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<2; i1++)
		{	
			p[i2+nb][i1+nb]=rwbndr[i1+4*i2];
			p[i2+nb][i1+nz+nb-2]=rwbndr[i1+2+4*i2];
		}
		for(i2=0; i2<2; i2++)
		for(i1=0; i1<nz; i1++)
		{
			p[i2+nb][i1+nb]=rwbndr[4*nz+i1+nz*i2];
			p[i2+nb+nx-2][i1+nb]=rwbndr[4*nz+i1+nz*(i2+2)];
		}
	}else{	// write effective boundary into bndr	
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<2; i1++)
		{	
			rwbndr[i1+4*i2]=p[i2+nb][i1+nb];
			rwbndr[i1+2+4*i2]=p[i2+nb][i1+nz+nb-2];
		}
		for(i2=0; i2<2; i2++)
		for(i1=0; i1<nz; i1++)
		{
			rwbndr[4*nz+i1+nz*i2]=p[i2+nb][i1+nb];
			rwbndr[4*nz+i1+nz*(i2+2)]=p[i2+nb+nx-2][i1+nb];
		}
	}
}



void step_forward(float **u0, float **u1, float **vv, bool adj)
/*< forward step for wave propagation >*/
{
	int i1, i2;

	if(adj){
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(i2,i1)		    	\
    shared(nzpad,nxpad,u1,vv,u0,c0,c11,c12,c21,c22)
#endif
		for (i2=2; i2<nxpad-2; i2++) 
		for (i1=2; i1<nzpad-2; i1++) 
		{
			u0[i2][i1]=2.*u1[i2][i1]-u0[i2][i1]+c0*vv[i2][i1  ]*u1[i2][i1  ]+
				c11*(vv[i2][i1-1]*u1[i2][i1-1]+vv[i2][i1+1]*u1[i2][i1+1])+
				c12*(vv[i2][i1-2]*u1[i2][i1-2]+vv[i2][i1+2]*u1[i2][i1+2])+
				c21*(vv[i2-1][i1]*u1[i2-1][i1]+vv[i2+1][i1]*u1[i2+1][i1])+
				c22*(vv[i2-2][i1]*u1[i2-2][i1]+vv[i2+2][i1]*u1[i2+2][i1]);
		}
	}else{
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(i2,i1)		    	\
    shared(nzpad,nxpad,u1,vv,u0,c0,c11,c12,c21,c22)
#endif
		for (i2=2; i2<nxpad-2; i2++) 
		for (i1=2; i1<nzpad-2; i1++) 
		{
			u0[i2][i1]=2.*u1[i2][i1]-u0[i2][i1]+vv[i2][i1]*(c0*u1[i2][i1]+
			c11*(u1[i2][i1-1]+u1[i2][i1+1])+c12*(u1[i2][i1-2]+u1[i2][i1+2])+
			c21*(u1[i2-1][i1]+u1[i2+1][i1])+c22*(u1[i2-2][i1]+u1[i2+2][i1]));
		}
	}
}


void apply_sponge(float **p0)
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

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns)
/*< shot/geophone position initialize
sxz/gxz; szbeg/gzbeg; sxbeg/gxbeg; jsz/jgz; jsx/jgx; ns/ng; >*/
{
	int is, sz, sx;

	for(is=0; is<ns; is++)
	{
		sz=szbeg+is*jsz;
		sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

void add_source(int *sxz, float **p, int ns, float *source, bool add)
/*< add seismic sources in grid >*/
{
	int is, sx, sz;

	if(add){/* add sources*/
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(is,sx,sz)		\
	shared(p,source,sxz,nb,ns,nz)
#endif
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx+nb][sz+nb]+=source[is];
		}
	}else{ /* subtract sources */
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(is,sx,sz)		\
	shared(p,source,sxz,nb,ns,nz)
#endif
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx+nb][sz+nb]-=source[is];
		}
	}
}


void prtm2d_init(bool verb_, bool csdgather_, float dz_, float dx_, float dt_, 
	int nz_, int nx_, int nb_, int nt_, int ns_, int ng_,
	int sxbeg_, int szbeg_, int jsx_, int jsz_, 
	int gxbeg_, int gzbeg_, int jgx_, int jgz_,
	float *wlt_, float **v0, float *mod_, float *dat_)
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

	verb=verb_;
	csdgather=csdgather_;
	ns=ns_;
	ng=ng_;
	nb=nb_;
	nz=nz_;
	nx=nx_;
	nt=nt_;
	sxbeg=sxbeg_;
	szbeg=szbeg_;
	jsx=jsx_;
	jsz=jsz_;
	gxbeg=gxbeg_;
	gzbeg=gzbeg_;
	jgx=jgx_;
	jgz=jgz_;
	wlt=wlt_;

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nm=nzpad*nxpad;
	nd=nt*nx;

    	/* allocate temporary arrays */
    	bndr=sf_floatalloc(nb);
    	sp0=sf_floatalloc2(nzpad,nxpad);
    	sp1=sf_floatalloc2(nzpad,nxpad);
    	gp0=sf_floatalloc2(nzpad,nxpad);
    	gp1=sf_floatalloc2(nzpad,nxpad);
    	vv=sf_floatalloc2(nzpad,nxpad);
	sxz=sf_intalloc(ns);
	gxz=sf_intalloc(ng);

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

	/* configuration of sources and geophones */
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_warning("sources exceeds the computing zone!"); exit(1);}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_warning("geophones exceeds the computing zone!"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ sf_warning("geophones exceeds the computing zone!"); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
}

void prtm2d_close()
/*< free allocated variables >*/
{
	free(bndr);
	free(*sp0); free(sp0);
	free(*sp1); free(sp1);
	free(*gp0); free(gp0);
	free(*gp1); free(gp1);
	free(*vv); free(vv);
	free(sxz);
	free(gxz);
}

void prtm2d_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< prtm2d linear operator: it may be parallized using MPI >*/
{
	int i1,i2,it,is,ig, gx, gz;
	if(nm!=nx*nz) sf_error("model size mismatch: %d!=%d",nm, nx*nz);
	if(nd!=nt*ng*ns) sf_error("data size mismatch: %d!=%d",nd,nt*ng*ns);
	sf_adjnull (adj, add, nm, nd, mod, dat); 

    	if(adj){// migration
		for(is=0; is<ns; is++){
			memset(sp0[0], 0, nzpad*nxpad*sizeof(float));
			memset(sp1[0], 0, nzpad*nxpad*sizeof(float));
			/* generate is-th source wavefield */
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
			}
			for(it=0; it<nt; it++){		
				bndr_rw(false, &rwbndr[it*4*(nx+nz)], sp1);	
				add_source(&sxz[is], sp1, 1, &wlt[it], true);
				step_forward(sp0, sp1, vv, false);
				apply_sponge(sp0);
				apply_sponge(sp1);
				ptr=sp0; sp0=sp1; sp1=ptr;	
			}

			memset(gp0[0], 0, nzpad*nxpad*sizeof(float));
			memset(gp1[0], 0, nzpad*nxpad*sizeof(float));
		    	for (it=nt-1; it >-1; it--) {// dat-->mod
				if(verb) sf_warning("%d;",it);
				
				/* reconstruct source wavefield */
				ptr=sp0; sp0=sp1; sp1=ptr;
				step_forward(sp0, sp1, vv, false);
				add_source(&sxz[is], sp1, 1, &wlt[it], false);
				bndr_rw(true, &rwbndr[it*4*(nx+nz)], sp1);
				
				/* backpropagate receiver wavefield */
				for(ig=0;ig<ng; ig++){
					gx=gxz[ig]/nz+nb;
					gz=gxz[ig]%nz+nb;
					gp1[gx][gz]+=dat[it+ig*nt+is*nt*ng];
				}
				step_forward(gp0, gp1, vv, false);
				apply_sponge(gp0); 
				apply_sponge(gp1); 
				ptr=gp0; gp0=gp1; gp1=ptr;

				for(i2=0; i2<nx; i2++)
				for(i1=0; i1<nz; i1++)
					mod[i1+nz*i2]+=sp1[i2+nb][i1+nb]*gp1[i2+nb][i1+nb];
		    	}
		}
    	}else{ // demigration, think of it as ajoint of migration 
		for(is=0; is<ns; is++){
			memset(sp0[0], 0, nzpad*nxpad*sizeof(float));
			memset(sp1[0], 0, nzpad*nxpad*sizeof(float));
			/* generate is-th source wavefield */
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
			}
			for(it=0; it<nt; it++){		
				bndr_rw(false, &rwbndr[it*4*(nx+nz)], sp1);
	
				add_source(&sxz[is], sp1, 1, &wlt[it], true);
				step_forward(sp0, sp1, vv, false);
				apply_sponge(sp0);
				apply_sponge(sp1);
				ptr=sp0; sp0=sp1; sp1=ptr;	
			}

			memset(gp0[0], 0, nzpad*nxpad*sizeof(float));
			memset(gp1[0], 0, nzpad*nxpad*sizeof(float));
		    	for (it=0; it<nt; it++) {//mod-->dat
				if(verb) sf_warning("%d;",it);
				
				/* reconstruct source wavefield */
				ptr=sp0; sp0=sp1; sp1=ptr;
				step_forward(sp0, sp1, vv, false);				
				add_source(&sxz[is], sp1, 1, &wlt[it], false);
				bndr_rw(true, &rwbndr[it*4*(nx+nz)], sp1);
					
				for(i2=0; i2<nx; i2++)
				for(i1=0; i1<nz; i1++)
					gp1[i2+nb][i1+nb]+=mod[i1+nz*i2];

				ptr=gp0; gp0=gp1; gp1=ptr;
				apply_sponge(gp0); 
				apply_sponge(gp1); 
				step_forward(gp0, gp1, vv, true);

				for(ig=0;ig<ng; ig++){
					gx=gxz[ig]/nz;
					gz=gxz[ig]%nz;
					dat[it+ig*nt+is*nt*ng]+=sp1[gx+nb][gz+nb]*gp1[gx+nb][gz+nb];
				}
		    	}
		}
    	}
}


