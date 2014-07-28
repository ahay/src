/* 2-D forward modeling to generate shot records 
Note: 	Clayton-Enquist absorbing boundary condition (A2) is applied!
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

  Reference: Clayton, Robert, and Björn Engquist. "Absorbing boundary 
	conditions for acoustic and elastic wave equations." Bulletin 
	of the Seismological Society of America 67.6 (1977): 1529-1540.
*/

#include <rsf.h>
#include <time.h>


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
    }

    for (ix=1; ix < nx-1; ix++) { 
	/* top boundary */
/*
	iz=0;
	v1=vv[ix][iz]*dtz; 
	v2=vv[ix][iz]*dtx;
	diff1=	(p1[ix][iz+1]-p1[ix][iz])-
		(p0[ix][iz+1]-p0[ix][iz]);
	diff2=	p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
	diff1*=v1;
 	diff2*=0.5*v2*v2;
	p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
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



int main(int argc, char* argv[])
{
	bool csdgather, chk;
	int nz, nx, nt, ns, ng, is, it, kt, distx, distz, sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
	int *sxz, *gxz;
  	float dx, dz, fm, dt, dtx, dtz, amp, tmp, totaltime=0	;
	float *trans, *wlt, *dobs, *bndr, **vv, **p0, **p1, **p2, **ptr=NULL;
	clock_t start, end;
	sf_file vinit, shots, check, time;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
    	shots=sf_output("out");  /* output image with correlation imaging condition */ 
	time=sf_output("time"); /* output total time */ 

    	/* get parameters for forward modeling */
    	if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");
    	if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

    	if(!sf_getbool("chk",&chk)) chk=false;
    	/*check whether GPU-CPU implementation coincide with each other or not */
	if(chk){
    		if (!sf_getint("kt",&kt))  kt=100;/* check it at it=100 */
		check=sf_output("check");/* output shotsnap for correctness checking*/
	}
	if (!sf_getfloat("amp",&amp)) amp=1000;
	/* maximum amplitude of ricker */
    	if (!sf_getfloat("fm",&fm)) fm=10;	
	/* dominant freq of ricker */
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
	sf_putint(time,"n1",1);
	sf_putint(time,"n2",1);

	dtx=dt/dx; 
	dtz=dt/dz; 

	wlt=(float*)malloc(nt*sizeof(float));
	bndr=(float*)malloc(nt*(2*nz+nx)*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	trans=(float*)malloc(ng*nt*sizeof(float));
	vv=sf_floatalloc2(nz, nx);
	p0=sf_floatalloc2(nz, nx);
	p1=sf_floatalloc2(nz, nx);
	p2=sf_floatalloc2(nz, nx);
	sxz=(int*)malloc(ns*sizeof(int));
	gxz=(int*)malloc(ng*sizeof(int));

	for(it=0; it<nt; it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
		wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
	}
	memset(bndr,0,nt*(2*nz+nx)*sizeof(float));
	memset(dobs,0,ng*nt*sizeof(float));
	memset(trans,0,ng*nt*sizeof(float));
	sf_floatread(vv[0],nz*nx,vinit);
	memset(p0[0],0,nz*nx*sizeof(float));
	memset(p1[0],0,nz*nx*sizeof(float));
	memset(p2[0],0,nz*nx*sizeof(float));
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
 	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);

	for(is=0; is<ns; is++)
	{
		start = clock();

		if (csdgather)	{
			gxbeg=sxbeg+is*jsx-distx;
			sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
		}
		memset(p0[0],0,nz*nx*sizeof(float));
		memset(p1[0],0,nz*nx*sizeof(float));
		memset(p2[0],0,nz*nx*sizeof(float));
		for(it=0; it<nt; it++)
		{
			add_source(p1, &wlt[it], &sxz[is], 1, nz, true);			
			step_forward(p0, p1, p2, vv, dtz, dtx, nz, nx);
			ptr=p0; p0=p1; p1=p2; p2=ptr;
			record_seis(&dobs[it*ng], gxz, p0, ng, nz);

			if(it==kt){
				sf_floatwrite(p0[0],nz*nx, check);
			}
		}
		matrix_transpose(dobs, trans, ng, nt);
		sf_floatwrite(trans,ng*nt,shots);
		
 		end = clock();
 		sf_warning("shot %d finished: %f (s)", is+1,((float)(end-start))/CLOCKS_PER_SEC); 
		totaltime+=((float)(end-start))/CLOCKS_PER_SEC;
	}
	totaltime/=ns;
	sf_floatwrite(&totaltime,1,time);

	free(sxz);
	free(gxz);
	free(bndr);
	free(dobs);
	free(trans);
	free(wlt);
	free(*vv); free(vv);
	free(*p0); free(p0);
	free(*p1); free(p1);
	free(*p2); free(p2);


    	exit(0);
}

