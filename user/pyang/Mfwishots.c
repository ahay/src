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
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int nb, nz, nx, nt, ns, ng, nzpad, nxpad;
static int *sxz, *gxz;
static float dz, dx, dt, fm, c0, c11, c12, c21, c22;
static float *wlt, *bndr, *dobs;
static float **v0,**vv, **p0, **p1, **ptr=NULL;

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
	    b[nb+ix][iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<nxpad; ix++) 
	for (iz=0; iz<nb;    iz++) 
	    b[ix][nzpad-iz-1] = b[ix][nzpad-nb-1];/* bottom*/

    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[ix 	 ][iz] = b[nb  		][iz];/* left */
	    b[nxpad-ix-1 ][iz] = b[nxpad-nb-1	][iz];/* right */
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
	    a[ix][iz]=b[nb+ix][iz] ;
	}
    }
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
	    p0[ix][ibz] *= w; /* bottom sponge */
	}

	ibx = nxpad-ib-1;
	for(iz=0; iz<nzpad; iz++) {
	    p0[ib ][iz] *= w; /*   left sponge */
	    p0[ibx][iz] *= w; /*  right sponge */
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
			sz=sxz[is]%nz;
			p[sx][sz]+=source[is];
		}
	}else{
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz;
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
		gz=gxz[ig]%nz;
		seis_it[ig]=p[gx][gz];
	}
}

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
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
	bool csdgather;
	int iz, ix, ib,is, it, jsx, jsz, jgx, jgz, sxbeg, szbeg, gxbeg, gzbeg, distx, distz;
	float *trans, tmp, amp;
	clock_t start, end;
	sf_file vinit, shots;

    	/* initialize Madagascar */
    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

    	/*< set up I/O files >*/
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
    	shots=sf_output("out");  /* output image with correlation imaging condition */ 

    	/* get parameters for forward modeling */
    	if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");
    	if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

	if (!sf_getfloat("amp",&amp)) amp=1000;
	/* maximum amplitude of ricker */
    	if (!sf_getfloat("fm",&fm)) fm=10;	
	/* dominant freq of ricker */
    	if (!sf_getint("nb",&nb))   nb=30;
	/* thickness of sponge ABC  */
    	if (!sf_getfloat("dt",&dt)) sf_error("no dt");	
	/* time interval */
    	if (!sf_getint("nt",&nt))   sf_error("no nt");	
	/* total modeling time steps */
    	if (!sf_getint("ns",&ns))   sf_error("no ns");	
	/* total shots */
    	if (!sf_getint("ng",&ng))   sf_error("no ng");	
	/* total receivers in each shot */	
    	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
	/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0;
	/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;
	/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0;
	/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
	/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
	/* z-begining index of receivers, starting from 0 */
	if (!sf_getbool("csdgather",&csdgather)) csdgather=false;
	/* default, common shot-gather; if n, record at every point*/

	sf_putint(shots,"n1",nt);	
	sf_putint(shots,"n2",ng);
	sf_putint(shots,"n3",ns);
	sf_putfloat(shots,"d1",dt);
	sf_putfloat(shots,"d2",jgx*dx);
	sf_putfloat(shots,"o1",0);
	sf_putstring(shots,"label1","Time");
	sf_putstring(shots,"label2","Lateral");
	sf_putstring(shots,"label3","Shot");
	sf_putstring(shots,"unit1","sec");
	sf_putstring(shots,"unit2","m");
	sf_putfloat(shots,"amp",amp);
	sf_putfloat(shots,"fm",fm);
	sf_putint(shots,"ng",ng);
	sf_putint(shots,"szbeg",szbeg);
	sf_putint(shots,"sxbeg",sxbeg);
	sf_putint(shots,"gzbeg",gzbeg);
	sf_putint(shots,"gxbeg",gxbeg);
	sf_putint(shots,"jsx",jsx);
	sf_putint(shots,"jsz",jsz);
	sf_putint(shots,"jgx",jgx);
	sf_putint(shots,"jgz",jgz);
	sf_putint(shots,"csdgather",csdgather?1:0);

	nzpad=nz+nb;
	nxpad=nx+2*nb;

	/*< initialize 4-th order FD coefficients >*/
	tmp = 1.0/(dz*dz);
	c11 = 4.0*tmp/3.0;
	c12= -tmp/12.0;
	tmp = 1.0/(dx*dx);
	c21 = 4.0*tmp/3.0;
	c22= -tmp/12.0;
	c0=-2.0*(c11+c12+c21+c22);

	v0=sf_floatalloc2(nz,nx); 
	vv=sf_floatalloc2(nzpad, nxpad);
	p0=sf_floatalloc2(nzpad, nxpad);
	p1=sf_floatalloc2(nzpad, nxpad);
	dobs=(float*)malloc(ng*nt*sizeof(float));
	trans=(float*)malloc(ng*nt*sizeof(float));
	bndr=sf_floatalloc(nb);
	wlt=sf_floatalloc(nt);
	sxz=sf_intalloc(ns);
	gxz=sf_intalloc(ng);

	sf_floatread(v0[0],nz*nx,vinit);
	expand2d(vv, v0);
	memset(dobs,0,ng*nt*sizeof(float));
	memset(trans,0,ng*nt*sizeof(float));
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
		wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
	}
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);

	for(is=0; is<ns; is++)
	{
		start = clock();

		if (csdgather){
			gxbeg=sxbeg+is*jsx-distx;
			sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
		}
		memset(p0[0],0,nzpad*nxpad*sizeof(float));
		memset(p1[0],0,nzpad*nxpad*sizeof(float));
		for(it=0; it<nt; it++)
		{
			add_source(&sxz[is], p1, 1, &wlt[it], true);
			step_forward(p0, p1);
			apply_sponge(p1);
			apply_sponge(p0);
			ptr=p0; p0=p1; p1=ptr;

			record_seis(&dobs[it*ng], gxz, p0, ng);
		}
		matrix_transpose(dobs, trans, ng, nt);
		sf_floatwrite(trans, ng*nt, shots);
		
 		end = clock();
 		sf_warning("shot %d finished: %f", is+1,((float)(end-start))/CLOCKS_PER_SEC); 
	}

	free(sxz);
	free(gxz);
	free(dobs);
	free(trans);
	free(wlt);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*p0); free(p0);
	free(*p1); free(p1);
	free(bndr);

    	exit(0);
}

