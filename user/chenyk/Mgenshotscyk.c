/* Generate shots for FWI test 
 */
/*
  Copyright (C) 2014  The University of Texas at Austin

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

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
}


void step_forward(float **p0, float **p1, float **p2, float **vv, float dtz, float dtx, int nz, int nx)
/*< forward modeling step, Clayton-Enquist ABC incorporated >*/
{
    int ix,iz;
    float v1,v2,diff1,diff2;

    for (ix=0; ix < nx; ix++) 
    for (iz=0; iz < nz; iz++) 
    {
	    v1=vv[ix][iz]*dtz; v1=v1*v1;
	    v2=vv[ix][iz]*dtx; v2=v2*v2;
	    diff1=diff2=-2.0*p1[ix][iz];
	    diff1+=(iz-1>=0)?p1[ix][iz-1]:0.0;
	    diff1+=(iz+1<nz)?p1[ix][iz+1]:0.0;
	    diff2+=(ix-1>=0)?p1[ix-1][iz]:0.0;
	    diff2+=(ix+1<nx)?p1[ix+1][iz]:0.0;
	    diff1*=v1;
	    diff2*=v2;
	    p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	    //if(fabsf(p2[ix][iz])>0.0) sf_warning("p2[%d][%d]=%g",ix,iz,p2[ix][iz]);
	 //   sf_warning("p2[ix][iz]=%g",p2[ix][iz]);
    }

    for (ix=1; ix < nx-1; ix++) { 
	/* top boundary */
/*
	iz=0;
	diff1=	(p1[ix][iz+1]-p1[ix][iz])-
		(p0[ix][iz+1]-p0[ix][iz]);
	diff2=	c21*(p1[ix-1][iz]+p1[ix+1][iz]) +
		c22*(p1[ix-2][iz]+p1[ix+2][iz]) +
		c20*p1[ix][iz];
	diff1*=sqrtf(vv[ix][iz])/dz;
	diff2*=vv[ix][iz]/(2.0*dx*dx);
	p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
*/
	/* bottom boundary */
	iz=nz-1;
	v1=vv[ix][iz]*dtz; 
	v2=vv[ix][iz]*dtx;
	diff1=-(p1[ix][iz]-p1[ix][iz-1])+(p0[ix][iz]-p0[ix][iz-1]);
	diff2=p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
	diff1*=v1;
 	diff2*=0.5*v2*v2;
	p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }

    for (iz=1; iz <nz-1; iz++){ 
	/* left boundary */
	ix=0;
	v1=vv[ix][iz]*dtz; 
	v2=vv[ix][iz]*dtx;
	diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	diff2=(p1[ix+1][iz]-p1[ix][iz])-(p0[ix+1][iz]-p0[ix][iz]);
	diff1*=0.5*v1*v1;
	diff2*=v2;
	p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	/* right boundary */
	ix=nx-1;
	v1=vv[ix][iz]*dtz; 
	v2=vv[ix][iz]*dtx;
	diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
	diff2=-(p1[ix][iz]-p1[ix-1][iz])+(p0[ix][iz]-p0[ix-1][iz]);
 	diff1*=0.5*v1*v1;
	diff2*=v2;
	p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }  
}

void add_source(float **p, float *source, int *sxz, int ns, int nz, bool add)
/*< add/subtract seismic sources >*/
{
	int is, sx, sz;
	if(add){
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]+=source[is];
		}
	}else{
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz;
			sz=sxz[is]%nz;
			p[sx][sz]-=source[is];
		}
	}
}

void record_seis(float *seis_it, int *gxz, float **p, int ng, int nz)
/*< record seismogram at time it into a vector length of ng >*/
{
	int ig, gx, gz;
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz;
		gz=gxz[ig]%nz;
		seis_it[ig]=p[gx][gz];
	}
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz)
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



void rw_bndr(float *bndr, float **p, int nz, int nx, bool write)
/*< if write==true, write/save boundaries out of variables;
 else  read boundaries into variables (for 2nd order FD) >*/
{
	int i;
	if(write){
		for(i=0; i<nz; i++)
		{
			bndr[i]=p[0][i];
			bndr[i+nz]=p[nx-1][i];
		}
		for(i=0; i<nx; i++) bndr[i+2*nz]=p[i][nz-1];
	}else{
		for(i=0; i<nz; i++)
		{
			p[0][i]=bndr[i];
			p[nx-1][i]=bndr[i+nz];
		}
		for(i=0; i<nx; i++) p[i][nz-1]=bndr[i+2*nz];
	}
}

int main(int argc, char *argv[])
{
	/* variables on host */
	bool verb, csdgather;
	int is, it, distx, distz;
	int nz, nx, nt, ns, ng;
	int sxbeg, szbeg, gxbeg, gzbeg, jsx, jsz, jgx, jgz;/*  parameters of acquisition geometery */
	float dx, dz, fm, dt, dtx, dtz, tmp, amp;
	float *dobs, *wlt, *bndr, *trans;
	int *sxz, *gxz;		
	float **vv, **sp0, **sp1, **sp2, **gp0, **gp1, **gp2, **g0, **g1, **ptr=NULL;	
	sf_file vel, shots;//, wfd;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/* set up I/O files */
    	vel=sf_input ("in");   /* initial velocity model, unit=m/s */
    	shots=sf_output("out"); /* updated velocity in iterations */ 
 //   	wfd=sf_output("wfd"); /* updated velocity in iterations */ 

    	
    	/* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;/* vebosity */	
    	if (!sf_histint(vel,"n1",&nz)) sf_error("no n1");/* nz */
    	if (!sf_histint(vel,"n2",&nx)) sf_error("no n2");/* nx */
    	if (!sf_histfloat(vel,"d1",&dz)) sf_error("no d1");/* dz */
   	if (!sf_histfloat(vel,"d2",&dx)) sf_error("no d2");/* dx */

   	if (!sf_getint("nt",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_getint("ng",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_getint("ns",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_getfloat("dt",&dt)) sf_error("no dt");
	/* time sampling interval */
 	if (!sf_getfloat("amp",&amp)) amp=1000;
	/* maximum amplitude of ricker */
    	if (!sf_getfloat("fm",&fm)) fm=10;	
	/* dominant freq of ricker */
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
	
//	int nwfd;
//	nwfd=(nt/100);
//	sf_putint(wfd,"n3",nwfd); 
//	sf_putint(wfd,"o3",0);
	 	
	dtx=dt/dx; 
	dtz=dt/dz; 


	vv=sf_floatalloc2(nz, nx);/* updated velocity */
	sp0=sf_floatalloc2(nz, nx);/* source wavefield p0 */
	sp1=sf_floatalloc2(nz, nx);/* source wavefield p1 */
	sp2=sf_floatalloc2(nz, nx);/* source wavefield p2 */
	gp0=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p0 */
	gp1=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p1 */
	gp2=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p2 */
	g0=sf_floatalloc2(nz, nx);/* gradient at previous step */
	g1=sf_floatalloc2(nz, nx);/* gradient at curret step */




	wlt=(float*)malloc(nt*sizeof(float));/* ricker wavelet */
	sxz=(int*)malloc(ns*sizeof(int)); /* source positions */
	gxz=(int*)malloc(ng*sizeof(int)); /* geophone positions */
	bndr=(float*)malloc(nt*(2*nz+nx)*sizeof(float));/* boundaries for wavefield reconstruction */
	trans=(float*)malloc(ng*nt*sizeof(float));/* transposed one shot */
	dobs=(float*)malloc(ng*nt*sizeof(float));/* observed seismic data */



	/* initialize varibles */
	sf_floatread(vv[0], nz*nx, vel);
	memset(sp0[0], 0, nz*nx*sizeof(float));
	memset(sp1[0], 0, nz*nx*sizeof(float));
	memset(sp2[0], 0, nz*nx*sizeof(float));
	memset(gp0[0], 0, nz*nx*sizeof(float));
	memset(gp1[0], 0, nz*nx*sizeof(float));
	memset(gp2[0], 0, nz*nx*sizeof(float));
	memset(g0[0], 0, nz*nx*sizeof(float));
	memset(g1[0], 0, nz*nx*sizeof(float));




	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
			//sf_warning("wlt[%d]=%g",it,wlt[it]);
	}
	

	
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_warning("sources exceeds the computing zone!\n"); exit(1);}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1); }
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1); }
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
	memset(bndr, 0, nt*(2*nz+nx)*sizeof(float));
	memset(dobs, 0, ng*nt*sizeof(float));	

		memcpy(g0[0], g1[0], nz*nx*sizeof(float));
		memset(g1[0], 0, nz*nx*sizeof(float));
		for(is=0;is<ns;is++)
		{
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
			}
			memset(sp0[0], 0, nz*nx*sizeof(float));
			memset(sp1[0], 0, nz*nx*sizeof(float));
			for(it=0; it<nt; it++)
			{
				add_source(sp1, &wlt[it], &sxz[is], 1, nz, true);			
				step_forward(sp0, sp1, sp2, vv, dtz, dtx, nz, nx);
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
				record_seis(&dobs[it*ng], gxz, sp0, ng, nz);
//				if(it%100==0) sf_floatwrite(sp0[0],nz*nx,wfd);
			}
			matrix_transpose(dobs,trans,ng,nt);
			sf_floatwrite(trans,ng*nt,shots);	
			if(verb)sf_warning("Shot %d/%d is finishes!",is+1,ns);
		
		}

	exit(0);
}
