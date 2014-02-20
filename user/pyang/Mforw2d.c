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

static int nz, nx, nt;
static float dz, dx, dt, fm, c10, c11, c12, c20, c21, c22;
static float *bndr;
static float **vv, **p0, **p1, **p2, **ptr=NULL;

void fd2d_init(float **v0)
/* allocate and initialize variables */
{
	int iz, ix;
	float tmp;

#ifdef _OPENMP
    	omp_init();
#endif
	
	vv=sf_floatalloc2(nz, nx);
	p0=sf_floatalloc2(nz, nx);
	p1=sf_floatalloc2(nz, nx);
	p2=sf_floatalloc2(nz, nx);

	for(ix=0;ix<nx;ix++){
	    for(iz=0;iz<nz;iz++){
		tmp=v0[ix][iz]*dt;
		vv[ix][iz]=tmp*tmp;// vv=vv^2*dt^2
	    }
	}

	/*< initialize 4-th order fd coefficients >*/
    	c11 = 4.0/3.0;
    	c12 =  -1.0/12.0;
    	c21 = 4.0/3.0;
    	c22 =  -1.0/12.0;
    	c10  = -2.0 * (c11+c12);
    	c20  = -2.0 * (c21+c22);
}


void wavefield_init(float **p0, float**p1, float **p2)
{
	memset(p0[0],0,nz*nx*sizeof(float));
	memset(p1[0],0,nz*nx*sizeof(float));
	memset(p2[0],0,nz*nx*sizeof(float));
}

void step_forward(float **p0, float **p1, float**p2)
{
	int ix,iz;
	float diff1,diff2;
	float _dz2=1.0/(dz*dz);
	float _dx2=1.0/(dx*dx);

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
    private(ix,iz,diff1,diff2)				\
    shared(vv,p2,p1,p0,_dz2,_dx2,c10,c11,c12,c20,c21,c22,nz,nx)  
#endif	
/*
	for (ix=2; ix < nx-2; ix++) 
	for (iz=2; iz < nz-2; iz++) 
	{
		diff1 =	c10*p1[ix][iz]+
			c11*(p1[ix][iz-1]+p1[ix][iz+1])+
			c12*(p1[ix][iz-2]+p1[ix][iz+2]);
		diff2 =	c20*p1[ix][iz]+
			c21*(p1[ix-1][iz]+p1[ix+1][iz])+
			c22*(p1[ix-2][iz]+p1[ix+2][iz]);
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]
			+vv[ix][iz]*(_dz2*diff1+_dx2*diff2);
	}
*/
	for (ix=1; ix < nx-1; ix++) 
	for (iz=1; iz < nz-1; iz++) 
	{
		diff1 =	p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
		diff2 =	p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]
			+vv[ix][iz]*(_dz2*diff1+_dx2*diff2);
	}
}


void apply_abc(float **p0, float **p1, float **p2)
{
	int ix,iz;
	float diff1,diff2;

    for (ix=2; ix < nx-2; ix++) { 
	//top boundary
/*
	for (iz=0; iz < 2; iz++){
	    	diff1=	(p1[ix][iz+1]-p1[ix][iz])-
			(p0[ix][iz+1]-p0[ix][iz]);
	    	diff2=	c21*(p1[ix-1][iz]+p1[ix+1][iz]) +
			c22*(p1[ix-2][iz]+p1[ix+2][iz]) +
			c20*p1[ix][iz];
		diff1*=sqrtf(vv[ix][iz])/dz;
 	    	diff2*=vv[ix][iz]/(2.0*dx*dx);
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	} 
*/
	p2[ix][0]=0;
	//bottom boundary
	for (iz=nz-2; iz < nz; iz++){
	    	diff1=-(p1[ix][iz]-p1[ix][iz-1])+(p0[ix][iz]-p0[ix][iz-1]);
	    	diff2=p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
		diff1*=sqrtf(vv[ix][iz])/dz;
 	    	diff2*=vv[ix][iz]/(2.0*dx*dx);
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
    }

    for (iz=2; iz <nz-2; iz++){ 
	//left boundary
    	for(ix=0;ix<2;ix++){
	    	diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	    	diff2=(p1[ix+1][iz]-p1[ix][iz])-(p0[ix+1][iz]-p0[ix][iz]);
 	    	diff1*=vv[ix][iz]/(2.0*dz*dz);
		diff2*=sqrtf(vv[ix][iz])/dx;
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
	//right boundary
    	for(ix=nx-2;ix<nx;ix++){
	    	diff1=	p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	    	diff2=	-(p1[ix][iz]-p1[ix-1][iz])+(p0[ix][iz]-p0[ix-1][iz]);
 	    	diff1*=vv[ix][iz]/(2.0*dz*dz);
		diff2*=sqrtf(vv[ix][iz])/dx;
		p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	}
    }  
}


void fd2d_close()
/* free the allocated variables */
{
	free(*vv); free(vv);
	free(*p0); free(p0);
	free(*p1); free(p1);
	free(*p2); free(p2);
}


void add_source(int *sxz, float **p, int ns, float *source, bool add)
{
	if(add){
		for(int is=0;is<ns; is++){
			int sx=sxz[is]/nz;
			int sz=sxz[is]%nz;
			p[sx][sz]+=source[is];
		}
	}else{
		for(int is=0;is<ns; is++){
			int sx=sxz[is]/nz;
			int sz=sxz[is]%nz;
			p[sx][sz]-=source[is];
		}
	}
}

void record_seis(float *seis_it, int *gxz, float **p, int ng)
/* record seismogram at time it into a vector length of ng */
{
	for(int ig=0;ig<ng; ig++)
	{
		int gx=gxz[ig]/nz;
		int gz=gxz[ig]%nz;
		seis_it[ig]=p[gx][gz];
	}
}


void matrix_transpose(float *matrix, int n1, int n2)
/* transpose a matrix n1xn2 into n2xn1 */
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
/* shot/geophone position initialize */
{
	for(int is=0; is<ns; is++)
	{
		int sz=szbeg+is*jsz;
		int sx=sxbeg+is*jsx;
		sxz[is]=sz+nz*sx;
	}
}

void rw_bndr(float **p, float *bndr, bool read)
/*< if read==true, read boundaries into variables;
 else write/save boundaries (for 2nd order FD) >*/
{
	if(read){
		for(int i=0; i<nz; i++)
		{
			p[0][i]=bndr[i];
			p[nx-1][i]=bndr[i+nz];
		}
		for(int i=0; i<nx; i++) p[i][nz-1]=bndr[i+2*nz];
	}else{
		for(int i=0; i<nz; i++)
		{
			bndr[i]=p[0][i];
			bndr[i+nz]=p[nx-1][i];
		}
		for(int i=0; i<nx; i++) bndr[i+2*nz]=p[i][nz-1];
	}
}


int main(int argc, char* argv[])
{
	int ns, ng, jsx, jsz, jgx, jgz, sxbeg, szbeg, gxbeg, gzbeg;
	float **v0,*dobs, *dcal;	
	sf_file Fv, Fs;

    	sf_init(argc,argv);

	Fv = sf_input("in");/* veloctiy model */
	Fs = sf_output("out");/* shot records */

    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
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

	bndr=(float*)malloc(nt*(2*nz+nx)*sizeof(float));
	dobs=(float*)malloc(ns*ng*nt*sizeof(float));
	dcal=(float*)malloc(ng*nt*sizeof(float));
	memset(bndr,0,nt*(2*nz+nx)*sizeof(float));
	memset(dobs,0,ns*ng*nt*sizeof(float));
	memset(dcal,0,ng*nt*sizeof(float));
	v0=sf_floatalloc2(nz,nx); 
	sf_floatread(v0[0],nz*nx,Fv);
	fd2d_init(v0);

	float *wlt=(float*)malloc(nt*sizeof(float));
	for(int it=0; it<nt; it++){
		float a=SF_PI*fm*(it*dt-1.0/fm);a=a*a;
		wlt[it]=(1.0-2.0*a)*expf(-a);
	}

	int *sxz, *gxz;
	sxz=(int*)malloc(ns*sizeof(int));
	gxz=(int*)malloc(ng*sizeof(int));
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);

	for(int is=0; is<ns; is++)
	{
		wavefield_init(p0, p1, p2);
		for(int it=0; it<nt; it++)
		{
			add_source(&sxz[is], p1, 1, &wlt[it], true);
			step_forward(p0, p1, p2);
			apply_abc(p0,p1,p2);
			ptr=p0; p0=p1; p1=p2; p2=ptr;
			record_seis(&dcal[it*ng], gxz, p0, ng);
			rw_bndr(p0, &bndr[it*(2*nz+nx)], false);

		}
		matrix_transpose(dcal, ng, nt);
		memcpy(&dobs[is*ng*nt], dcal, ng*nt*sizeof(float));
/*
		ptr=p0; p0=p1; p1=ptr;
		for(int it=nt-1; it>-1; it--)
		{
			rw_bndr(p1, &bndr[it*(2*nz+nx)], true);
			step_forward(p0, p1, p2);
			add_source(&sxz[is], p1, 1, &wlt[it], false);
			ptr=p0; p0=p1; p1=p2; p2=ptr;

			if (it==200) sf_floatwrite(p0[0],nz*nx,Fs);
		}
*/
	}
	sf_floatwrite(dobs, ns*ng*nt, Fs);



	free(sxz);
	free(gxz);
	free(bndr);
	free(dobs);
	free(dcal);
	free(wlt);
	free(*v0); free(v0);
	fd2d_close();

    	exit(0);
}

