/* 2-D prestack LSRTM linear operator using wavefield reconstruction method
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

  Comments: 
  Operator modified to pass the adjoint. Added adjoint tests as well. 
  Daniel Trad, Uinversity of Calgary, 2017
*/
#include <rsf.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "prtm2d.h"
static float ***sp0array;
static bool fromBoundary;
static bool csdgather, verb;
static int nzpad, nxpad, nb, nz, nx, nt, ns, ng, sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz, distx, distz;
static int *sxz, *gxz;
static float c0, c11, c21, c12, c22;
static float *wlt, *bndr,*rwbndr, *mod, *dat;
static float **sp0, **sp1, **gp0, **gp1, **vv, **ptr=NULL;
int getTime(struct timeval t1, struct timeval  t2){
    // return time in milliseconds;
    return ((t2.tv_sec - t1.tv_sec)* 1000 + (t2.tv_usec - t1.tv_usec)/1000);
}
float compareWavefields(int it, float **u){
    ///*< compare current wavefield with true wavefield saved in sp0array >*/

    float sum =0;
    for (int ix=nb;ix<nxpad-nb;ix++){
	for (int iz=nb;iz<nzpad-nb;iz++){
	    if (gp1[ix][iz]) sum += fabs((u[ix][iz]-sp0array[it][ix][iz]));
	}
    }
    return sum;
}
void rw_snapshot(float** p, int it, bool read){
    // read/write snapshot completely 
    // if read=true read, else write
  
    if (!read){
	for (int ix=0;ix< nxpad;ix++){
	    for (int iz=0;iz< nzpad;iz++){
		sp0array[it][ix][iz]=p[ix][iz];
	    }
	}
    }
    else{
	for (int ix=0;ix< nxpad;ix++){
	    for (int iz=0;iz< nzpad;iz++){
		p[ix][iz]=sp0array[it][ix][iz];
	    }
	}
    }

}

void boundary_rw(float **p, float *spo, bool read)
/* read/write using effective boundary saving strategy: 
   if read=true, read the boundary out; else save/write the boundary*/
{
    int ix,iz;

    if (read){
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(p,spo,nx,nz,nb)  
#endif	
	for(ix=0; ix<nx; ix++)
	    for(iz=0; iz<2; iz++)
	    {
		p[ix+nb][iz+nb]=spo[iz+4*ix];
		p[ix+nb][iz-2+nz+nb]=spo[iz+2+4*ix];
	    }
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(p,spo,nx,nz,nb)  
#endif	
	for(iz=0; iz<nz; iz++)
	    for(ix=0; ix<2; ix++)
	    {
		p[ix+nb][iz+nb]=spo[4*nx+iz+nz*ix];
		p[ix-2+nx+nb][iz+nb]=spo[4*nx+iz+nz*(ix+2)];
	    }
    }else{
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(p,spo,nx,nz,nb)  
#endif	
	for(ix=0; ix<nx; ix++)
	    for(iz=0; iz<2; iz++)
	    {
		spo[iz+4*ix]=p[ix+nb][iz+nb];
		spo[iz+2+4*ix]=p[ix+nb][iz-2+nz+nb];
	    }
#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(p,spo,nx,nz,nb)  
#endif	
	for(iz=0; iz<nz; iz++)
	    for(ix=0; ix<2; ix++)
	    {
		spo[4*nx+iz+nz*ix]=p[ix+nb][iz+nb];
		spo[4*nx+iz+nz*(ix+2)]=p[ix-2+nx+nb][iz+nb];
	    }
    }
}

void step_forward(float **u0, float **u1, float **vv, bool adj)
/*< forward step for wave propagation >*/
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
		u0[i2][i1]=2.*u1[i2][i1]-u0[i2][i1]+c0*vv[i2][i1]*u1[i2][i1]+
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


void apply_sponge(float **p0)
/*< apply sponge (Gaussian taper) absorbing boundary condition
  L=Gaussian taper ABC; L=L*, L is self-adjoint operator. >*/
{
    int ix,iz,ib,ibx,ibz;
    float w;

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
#pragma omp parallel for default(none)		\
    private(is,sx,sz)				\
    shared(p,source,sxz,nb,ns,nz)
#endif
	for(is=0;is<ns; is++){
	    sx=sxz[is]/nz;
	    sz=sxz[is]%nz;
	    p[sx+nb][sz+nb]+=source[is];
	}
    }else{ /* subtract sources */
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
    private(is,sx,sz)				\
    shared(p,source,sxz,nb,ns,nz)
#endif
	for(is=0;is<ns; is++){
	    sx=sxz[is]/nz;
	    sz=sxz[is]%nz;
	    p[sx+nb][sz+nb]-=source[is];
	}
    }
}


void expand2d(float** b, float** a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
    private(ix,iz)				\
    shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    b[nb+ix][nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<nxpad; ix++) {
	for (iz=0; iz<nb;    iz++) {
	    b[ix][      iz  ] = b[ix][nb  ];
	    b[ix][nzpad-iz-1] = b[ix][nzpad-nb-1];
	}
    }

    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[ix 	 ][iz] = b[nb  		][iz];
	    b[nxpad-ix-1 ][iz] = b[nxpad-nb-1	][iz];
	}
    }
}


void window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
    private(ix,iz)				\
    shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}

void prtm2d_init(bool verb_, bool csdgather_, bool fromBoundary_, float dz_, float dx_, float dt_, 
		 float amp, float fm, 
		 int nz_, int nx_, int nb_, int nt_, int ns_, int ng_, 
		 int sxbeg_, int szbeg_, int jsx_, int jsz_, 
		 int gxbeg_, int gzbeg_, int jgx_, int jgz_,
		 float **v0, float *mod_, float *dat_)
/*< allocate variables and initialize parameters >*/
{  
#ifdef _OPENMP
    omp_init();  /* initialize OpenMP support */ 
#endif
    float t;
    int i1, i2, it,ib;
    t = 1.0/(dz_*dz_);
    c11 = 4.0*t/3.0;
    c12= -t/12.0;

    t = 1.0/(dx_*dx_);
    c21 = 4.0*t/3.0;
    c22= -t/12.0;
    c0=-2.0*(c11+c12+c21+c22);

    fromBoundary=fromBoundary_; 
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

    nzpad=nz+2*nb;
    nxpad=nx+2*nb;

    /* allocate temporary arrays */
    bndr=sf_floatalloc(nb);
    sp0=sf_floatalloc2(nzpad,nxpad);
    sp1=sf_floatalloc2(nzpad,nxpad);
    gp0=sf_floatalloc2(nzpad,nxpad);
    gp1=sf_floatalloc2(nzpad,nxpad);
    vv=sf_floatalloc2(nzpad,nxpad);
    wlt=sf_floatalloc(nt);
    sxz=sf_intalloc(ns);
    gxz=sf_intalloc(ng);
    rwbndr=sf_floatalloc(nt*4*(nx+nz));
    if (!fromBoundary) sp0array=sf_floatalloc3(nzpad,nxpad,nt);

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
    for(it=0; it<nt; it++){
	t=SF_PI*fm*(it*dt_-1.0/fm);t=t*t;
	wlt[it]=amp*(1.0-2.0*t)*expf(-t);
    }

    /* configuration of sources and geophones */
    if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
    { sf_warning("sources exceeds the computing zone!"); exit(1);}
    sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
    distx=sxbeg-gxbeg;
    distz=szbeg-gzbeg;
    if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
    { sf_warning("geophones exceeds the computing zone!"); exit(1);}
    if (csdgather && !( (sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  
			&& (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	{
	sf_warning("geophones exceeds the computing zone!"); exit(1);
    }
    sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
}

void prtm2d_close()
/*< free allocated variables >*/
{
    free(bndr);
    if (!fromBoundary){ free(**sp0array);free(*sp0array);free(sp0array);}
    free(*sp0); free(sp0);
    free(*sp1); free(sp1);
    free(*gp0); free(gp0);
    free(*gp1); free(gp1);
    free(*vv); free(vv);
    free(wlt);
    free(sxz);
    free(gxz);
}

void prtm2d_shotwav(int is)
/*< prtm2d shot wavefield generated at front >*/
{
    int it;

    /* initialize is-th source wavefield Ps[] */
    memset(sp0[0], 0, nzpad*nxpad*sizeof(float));
    memset(sp1[0], 0, nzpad*nxpad*sizeof(float));
    memset(gp0[0], 0, nzpad*nxpad*sizeof(float));
    memset(gp1[0], 0, nzpad*nxpad*sizeof(float));

    for(it=0; it<nt; it++){			
	add_source(&sxz[is], sp1, 1, &wlt[it], true);
	step_forward(sp0, sp1, vv, false);
	ptr=sp0; sp0=sp1; sp1=ptr;
	apply_sponge(sp0);
	apply_sponge(sp1);
	rw_snapshot(sp0,it,false);
    }
}


void prtm2d_lop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< prtm2d linear operator >*/
{
    int i1,i2,it,is,ig, gx, gz;

    if(nm!=nx*nz) sf_error("model size mismatch: %d!=%d",nm, nx*nz);
    if(nd!=nt*ng*ns) sf_error("data size mismatch: %d!=%d",nd,nt*ng*ns);
    sf_adjnull(adj, add, nm, nd, mod, dat); 

    for(is=0; is<ns; is++) {/* it may be parallized using MPI */
	/* initialize is-th source wavefield Ps[] */
	memset(sp0[0], 0, nzpad*nxpad*sizeof(float));
	memset(sp1[0], 0, nzpad*nxpad*sizeof(float));
	memset(gp0[0], 0, nzpad*nxpad*sizeof(float));
	memset(gp1[0], 0, nzpad*nxpad*sizeof(float));
	if(adj){/* migration: mm=Lt dd */
	    //==========================================================
	    for(it=0; it<nt; it++){			
		if (fromBoundary) boundary_rw(sp1, &rwbndr[it*4*(nx+nz)], false);
		else rw_snapshot(sp1,it,false);

		add_source(&sxz[is], sp1, 1, &wlt[it], true);
		step_forward(sp0, sp1, vv, false);
		ptr=sp0; sp0=sp1; sp1=ptr;
		apply_sponge(sp0);
		apply_sponge(sp1);
	    }
	    ptr=sp0; sp0=sp1; sp1=ptr;
	    for (it=nt-1; it >-1; it--) {
		/* reverse time order, Img[]+=Ps[]* Pg[]; */
		for(i2=0; i2<nx; i2++)
		    for(i1=0; i1<nz; i1++)
			mod[i1+nz*i2]+=sp1[i2+nb][i1+nb]*gp1[i2+nb][i1+nb];

		/* backpropagate receiver wavefield */
		for(ig=0;ig<ng; ig++){
		    gx=gxz[ig]/nz;
		    gz=gxz[ig]%nz;
		    gp1[gx+nb][gz+nb]+=dat[it+ig*nt+is*nt*ng];
		}
		step_forward(gp0, gp1, vv, false);
		ptr=gp0; gp0=gp1; gp1=ptr;
		apply_sponge(gp0); 
		apply_sponge(gp1); 

		if (fromBoundary){ 
		    boundary_rw(sp1, &rwbndr[(it)*4*(nx+nz)], true);
		    step_forward(sp0, sp1, vv, false);
		    add_source(&sxz[is], sp1, 1, &wlt[it], false);
		    ptr=sp0; sp0=sp1; sp1=ptr;
		}else rw_snapshot(sp0,it,true); 
	    }
	}else{/* Born modeling/demigration: dd=L mm */	
	    //==========================================================
	    for(it=0; it<nt; it++){	
		/* forward time order, Pg[]+=Ps[]* Img[]; */	
		for(i2=0; i2<nx; i2++)
		    for(i1=0; i1<nz; i1++)
			gp1[i2+nb][i1+nb]+=sp1[i2+nb][i1+nb]*mod[i1+nz*i2];

		step_forward(gp0, gp1, vv, true);	
		ptr=gp0; gp0=gp1; gp1=ptr;
		apply_sponge(gp0); 
		apply_sponge(gp1); 
		for(ig=0;ig<ng; ig++){//record data
		    gx=gxz[ig]/nz;
		    gz=gxz[ig]%nz;
		    dat[it+ig*nt+is*nt*ng]+=gp1[gx+nb][gz+nb];
		}

		add_source(&sxz[is], sp1, 1, &wlt[it], true);
		step_forward(sp0, sp1, vv, false);
		ptr=sp0; sp0=sp1; sp1=ptr;
		apply_sponge(sp0);
		apply_sponge(sp1);
	    }	
	}
    }	 
}

void prtm2d_adjtest()
/*< prtm2d adjoint test >*/	  
{
    unsigned int i;
    float dotdata,dotmodel;

    unsigned int datasize =ng*nt*ns;
    unsigned int modelsize=nx*nz;
    float* m1=sf_floatalloc(modelsize);
    float* m2=sf_floatalloc(modelsize);
    float* d1=sf_floatalloc(datasize);
    float* d2=sf_floatalloc(datasize);
  
    for (i=0;i<datasize;i++)  d1[i]=(drand48()-0.5);
    for (i=0;i<modelsize;i++) m1[i]=(drand48()-0.5);

    memset(m2,0,modelsize*sizeof(float));
    memset(d2,0,datasize*sizeof(float));

    prtm2d_lop(true,false,nx*nz,nt*ng*ns,m2,d1);
    prtm2d_lop(false,false,nx*nz,nt*ng*ns,m1,d2);

    dotdata=0;
    for (i=0;i<datasize;i++)  dotdata+=(d1[i]*d2[i]);
    dotmodel=0;
    for (i=0;i<modelsize;i++) dotmodel+=(m1[i]*m2[i]);
  
    fprintf(stderr,"dot=%g, dotdata=%g, dotmodel=%g\n",dotdata/dotmodel,dotdata,dotmodel);

    free(m1);
    free(m2);
    free(d1);
    free(d2);
}
