/* 2-D forward modeling to generate shot records 
Note: 	Here, the sponge absorbing boundary condition is applied!
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

static int nb, nz, nx, nt, nzpad, nxpad;
static float dz, dx, dt, fm, c0, c11, c12, c21, c22;
static float *bndr;
static float **vv, **p0, **p1, **ptr=NULL;

void expand2d(float** b, float** a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(ix,iz)			\
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
#pragma omp parallel for default(none)	\
	private(ix,iz)			\
	shared(b,a,nb,nz,nx)
#endif
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[ix][iz]=b[nb+ix][nb+iz] ;
	}
    }
}


void wavefield_init(float **p0, float**p1)
/*< initialize wavefield >*/
{
	memset(p0[0],0,nzpad*nxpad*sizeof(float));
	memset(p1[0],0,nzpad*nxpad*sizeof(float));
}

void step_forward(float **p0, float **p1)
/*< forward modeling step >*/
{
	int ix,iz;
	float tmp;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz,tmp)			\
    shared(p1,p0,vv,nzpad,nxpad,c0,c11,c12,c21,c22)  
#endif	
	for (ix=2; ix < nxpad-2; ix++) 
	for (iz=2; iz < nzpad-2; iz++) 
	{
		tmp =	c0*p1[ix][iz]+
			c11*(p1[ix][iz-1]+p1[ix][iz+1])+
			c12*(p1[ix][iz-2]+p1[ix][iz+2])+
			c21*(p1[ix-1][iz]+p1[ix+1][iz])+
			c22*(p1[ix-2][iz]+p1[ix+2][iz]);
		p0[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*tmp;
	}
}

void apply_sponge(float**p0)
/*< apply absorbing boundary condition >*/
{
	int ix,iz,ib,ibx,ibz;
	float w;

#ifdef _OPENMP
#pragma omp parallel for		\
    private(ib,iz,ix,ibz,ibx,w)		\
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



void boundary_rw(float **p, float *spo, bool read)
/*< read/write using effective boundary saving strategy: 
if read=true, read the boundary out; else save/write the boundary >*/
{
	int ix,iz;
	if (read){
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(p,spo,nx,nz,nb)  
#endif	
		for(ix=0; ix<nx; ix++)
		for(iz=0; iz<2; iz++)
		{
			p[ix+nb][iz-2+nb]=spo[iz+4*ix];
			p[ix+nb][iz+nz+nb]=spo[iz+2+4*ix];
		}
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(p,spo,nx,nz,nb)  
#endif	
		for(iz=0; iz<nz; iz++)
		for(ix=0; ix<2; ix++)
		{
			p[ix-2+nb][iz+nb]=spo[4*nx+iz+nz*ix];
			p[ix+nx+nb][iz+nb]=spo[4*nx+iz+nz*(ix+2)];
		}
	}else{
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(p,spo,nx,nz,nb)  
#endif	
		for(ix=0; ix<nx; ix++)
		for(iz=0; iz<2; iz++)
		{
			spo[iz+4*ix]=p[ix+nb][iz-2+nb];
			spo[iz+2+4*ix]=p[ix+nb][iz+nz+nb];
		}
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(p,spo,nx,nz,nb)  
#endif	
		for(iz=0; iz<nz; iz++)
		for(ix=0; ix<2; ix++)
		{
			spo[4*nx+iz+nz*ix]=p[ix-2+nb][iz+nb];
			spo[4*nx+iz+nz*(ix+2)]=p[ix+nx+nb][iz+nb];
		}
	}
}

void add_source(int *sxz, float **p, int ns, float *source, bool add)
/*< add source term >*/
{
	int is, sx, sz;
	if(add){
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]+=source[is];
		}
	}else{
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]-=source[is];
		}
	}
}

void record_seis(float *seis_it, int *gxz, float **p, int ng)
/*< record seismogram at time it into a vector length of ng >*/
{
	int ig, gx, gz;
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz+nb;
		gz=gxz[ig]%nz+nb;
		seis_it[ig]=p[gx][gz];
	}
}

void matrix_transpose(float *matrix, int n1, int n2)
/*< transpose a matrix n1xn2 into n2xn1 >*/
{
	float *tmp=(float*)malloc(n1*n2*sizeof(float));
	if (tmp==NULL) {printf("out of memory!"); exit(1);}
	for(int i2=0; i2<n2; i2++)
	for(int i1=0; i1<n1; i1++)
	{
		tmp[i2+n2*i1]=matrix[i1+n1*i2];
	}
	memcpy(matrix, tmp, n1*n2*sizeof(float));
	free(tmp);
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns)
/*< shot/geophone position initialize >*/
{
	int is, sz, sx;
	for(is=0; is<ns; is++)
	{
		sz=szbeg+is*jsz;
		sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

int main(int argc, char* argv[])
{
	int iz, ix, ib,is, it, ns, ng, jsx, jsz, jgx, jgz, sxbeg, szbeg, gxbeg, gzbeg;
	int *sxz, *gxz;
	float tmp;
	float *wlt, *spo, *dobs, *dcal;
	float **v0;
	sf_file Fv, Fs;

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif
	Fv = sf_input("in");/* veloctiy model */
	Fs = sf_output("out");/* shot records */

    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    	if (!sf_getint("nb",&nb)) nb=20; /* thickness of sponge ABC */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=15.0; /*dominant freq of Ricker wavelet */
	if (!sf_getint("ns",&ns)) ns=1;	/* number of shots */
	if (!sf_getint("ng",&ng)) ng=nx;/* number of receivers */
	if (ng>nx) sf_error("make sure ng<=nx!");

    	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0;/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0;/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");/* z-begining index of receivers, starting from 0 */

	sf_putint(Fs,"n1",nt);
	sf_putint(Fs,"n2",ng);
    	sf_putint(Fs,"n3",ns);
    	sf_putfloat(Fs,"d1",dt);

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;

	/*< initialize 4-th order FD coefficients >*/
	tmp = 1.0/(dz*dz);
	c11 = 4.0*tmp/3.0;
	c12= -tmp/12.0;
	tmp = 1.0/(dx*dx);
	c21 = 4.0*tmp/3.0;
	c22= -tmp/12.0;
	c0=-2.0*(c11+c12+c21 + c22);


	spo=(float*)malloc(nt*4*(nx+nz)*sizeof(float));
	dobs=(float*)malloc(ns*ng*nt*sizeof(float));
	dcal=(float*)malloc(ng*nt*sizeof(float));
	vv=sf_floatalloc2(nzpad, nxpad);
	p0=sf_floatalloc2(nzpad, nxpad);
	p1=sf_floatalloc2(nzpad, nxpad);
	bndr=sf_floatalloc(nb);
	wlt=sf_floatalloc(nt);
	sxz=sf_intalloc(ns);
	gxz=sf_intalloc(ng);

	memset(spo,0,nt*4*(nx+nz)*sizeof(float));
	memset(dobs,0,ns*ng*nt*sizeof(float));
	memset(dcal,0,ng*nt*sizeof(float));
	v0=sf_floatalloc2(nz,nx); 
	sf_floatread(v0[0],nz*nx,Fv);
	expand2d(vv, v0);
	for(ib=0;ib<nb;ib++){
		tmp=0.015*(nb-ib);
		bndr[ib]=expf(-tmp*tmp);
	}
	for(ix=0;ix<nxpad;ix++){
	    for(iz=0;iz<nzpad;iz++){
		tmp=vv[ix][iz]*dt;
		vv[ix][iz]=tmp*tmp;// vv=vv^2*dt^2
	    }
	}	
	for(it=0; it<nt; it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
	}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);

	for(is=0; is<ns; is++)
	{
		wavefield_init(p0, p1);
		for(it=0; it<nt; it++)
		{
			add_source(&sxz[is], p1, 1, &wlt[it], true);
			step_forward(p0, p1);
			apply_sponge(p1);
			apply_sponge(p0);
			ptr=p0; p0=p1; p1=ptr;

			record_seis(&dcal[it*ng], gxz, p0, ng);
		}
		matrix_transpose(dcal, ng, nt);
		memcpy(&dobs[is*ng*nt], dcal, ng*nt*sizeof(float));

	}
	sf_floatwrite(dobs, ns*ng*nt, Fs);

	free(sxz);
	free(gxz);
	free(spo);
	free(dobs);
	free(dcal);
	free(wlt);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*p0); free(p0);
	free(*p1); free(p1);
	free(bndr);

    	exit(0);
}

