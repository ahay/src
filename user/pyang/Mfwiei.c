/* Envelope inversion of FWI
Note: 	Enquist absorbing boundary condition (A2) is applied!
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

  Important references:
    [1] Clayton, Robert, and Björn Engquist. "Absorbing boundary 
	conditions for acoustic and elastic wave equations." Bulletin 
	of the Seismological Society of America 67.6 (1977): 1529-1540.
    [2] Tarantola, Albert. "Inversion of seismic reflection data in the 
	acoustic approximation." Geophysics 49.8 (1984): 1259-1266.
    [3] Pica, A., J. P. Diet, and A. Tarantola. "Nonlinear inversion 
	of seismic reflection data in a laterally invariant medium." 
	Geophysics 55.3 (1990): 284-292.
    [4] Dussaud, E., Symes, W. W., Williamson, P., Lemaistre, L., 
	Singer, P., Denel, B., & Cherrett, A. (2008). Computational 
	strategies for reverse-time migration. In SEG Technical Program 
	Expanded Abstracts 2008 (pp. 2267-2271).
    [5] Hager, William W., and Hongchao Zhang. "A survey of nonlinear
	conjugate gradient methods." Pacific journal of Optimization 
	2.1 (2006): 35-58.
*/

#include <rsf.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hilbert.h"

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

void step_backward(float **illum, float **lap, float **p0, float **p1, float **p2, float **vv, float dtz, float dtx, int nz, int nx)
/*< step backward >*/
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
	    lap[ix][iz]=diff1+diff2;
	    diff1*=v1;
	    diff2*=v2;
	    p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
	    illum[ix][iz]+=p1[ix][iz]*p1[ix][iz];
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



void cal_gradient(float **grad, float **lap, float **gp, int nz, int nx)
/*< calculate gradient >*/
{
  int ix, iz;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      grad[ix][iz]+=lap[ix][iz]*gp[ix][iz];
    }
  }
}

void scale_gradient(float **grad, float **vv, float **illum, int nz, int nx, bool precon)
/*< scale gradient >*/
{
  int ix, iz;
  float a;
  for(ix=1; ix<nx-1; ix++){
    for(iz=1; iz<nz-1; iz++){
	a=vv[ix][iz];
	if (precon) a*=sqrtf(illum[ix][iz]+SF_EPS);/*precondition with residual wavefield illumination*/
	grad[ix][iz]*=2.0/a;
    }
  }

  for(ix=0; ix<nx; ix++){
	grad[ix][0]=grad[ix][1];
	grad[ix][nz-1]=grad[ix][nz-2];
  }

  for(iz=0; iz<nz; iz++){
	grad[0][iz]=grad[1][iz];
	grad[nx-1][iz]=grad[nx-2][iz];
  }
}

float cal_objective(sf_complex *adcal, sf_complex *adobs, int ng, int nt)
/*< calculate the value of envelope objective function >*/
{
  int i;
  float a, obj=0;

  for(i=0; i<ng*nt; i++){
    a=logf(fabsf(adcal[i])/(fabsf(adobs[i])+SF_EPS));
    obj+=a*a;
  }
  return obj;
}

float cal_beta(float **g0, float **g1, float **cg, int nz, int nx)
/*< calculate beta >*/
{
  int ix, iz;
  float a,b,c;

  a=b=c=0;
  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
	a += g1[ix][iz]*(g1[ix][iz]-g0[ix][iz]);// numerator of HS
	b += cg[ix][iz]*(g1[ix][iz]-g0[ix][iz]);// denominator of HS,DY
	c += g1[ix][iz]*g1[ix][iz];		// numerator of DY
    }
  }

  float	beta_HS=(fabsf(b)>0)?(a/b):0.0; 
  float beta_DY=(fabsf(b)>0)?(c/b):0.0;
  return SF_MAX(0.0, SF_MIN(beta_HS, beta_DY));
}


void cal_conjgrad(float **g1, float **cg, float beta, int nz, int nx)
/*< calculate conjugate gradient >*/
{
  int ix, iz;

  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      cg[ix][iz]=-g1[ix][iz]+beta*cg[ix][iz];
    }
  }
}

float cal_epsilon(float **vv, float **cg, int nz, int nx)
/*< calculate epsilcon >*/
{
  int ix, iz;
  float vvmax, cgmax;
  vvmax=cgmax=0.0;

  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      vvmax=SF_MAX(vvmax, fabsf(vv[ix][iz]));
      cgmax=SF_MAX(cgmax, fabsf(cg[ix][iz]));
    }
  }

  return 0.01*vvmax/(cgmax+SF_EPS);
}

void cal_vtmp(float **vtmp, float **vv, float **cg, float epsil, int nz, int nx)
/*< calculate temporary velcity >*/
{
  int ix, iz;

  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      vtmp[ix][iz]=vv[ix][iz]+epsil*cg[ix][iz];
    }
  }
}

void sum_alpha12(float *alpha1, float *alpha2, float *dcaltmp, float *dobs, float *derr, int ng)
/*< calculate numerator and denominator of alpha >*/
{
  int ig;
  float a, b, c;
  for(ig=0; ig<ng; ig++){
	c=derr[ig];
	a=dobs[ig]+c;/* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */
	b=dcaltmp[ig]-a;/* f(mk+epsil*cg)-f(mk) */
	alpha1[ig]-=b*c; 
	alpha2[ig]+=b*b; 
  }
}


float cal_alpha(float *alpha1, float *alpha2, float epsil, int ng)
/*< calculate alpha >*/
{
  int ig;
  float a,b;

  a=b=0;
  for(ig=0; ig<ng; ig++){
    a+=alpha1[ig];
    b+=alpha2[ig];
  }

  return (a*epsil/(b+SF_EPS));
}

void update_vel(float **vv, float **cg, float alpha, int nz, int nx)
/*< update velcity >*/
{
  int ix, iz;

  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){
      vv[ix][iz]+=alpha*cg[ix][iz];
    }
  }
}

void bell_smoothz(float **g, float **smg, int rbell, int nz, int nx)
/*< gaussian bell smoothing for z-axis >*/
{
	int ix, iz, i;
	float s;

	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		s=0.0;
		for(i=-rbell; i<=rbell; i++) if(iz+i>=0 && iz+i<nz) s+=expf(-(2.0*i*i)/rbell)*g[ix][iz+i];
		smg[ix][iz]=s;		
	}
}

void bell_smoothx(float **g, float **smg, int rbell, int nz, int nx)
/*< gaussian bell smoothing for x-axis >*/
{
	int ix, iz, i;
	float s;

	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		s=0.0;
		for(i=-rbell; i<=rbell; i++) if(ix+i>=0 && ix+i<nx) s+=expf(-(2.0*i*i)/rbell)*g[ix+i][iz];
		smg[ix][iz]=s;		
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


void cal_adjsource(sf_complex *adcal, sf_complex *adobs, float *derr, int ng, int nt)
/*< calculate adjoint source and store it in derr[] >*/
{
  int ig, it;
  float a,b,c,d;
  float *h1;
  sf_complex *h2;
  h1=(float *)malloc(ng*nt*sizeof(sf_complex));
  h2=(sf_complex*)malloc(ng*nt*sizeof(sf_complex));

  for(it=0; it<nt; it++){
	  for(ig=0; ig<ng; ig++){
	    a=fabsf(adobs[ig+it*ng])+SF_EPS;
	    b=fabsf(adcal[ig+it*ng])+SF_EPS;
	    c=logf(a/b)/(b*b);
	    //h1[ig+it*ng]=c*crealf(adcal[ig+it*ng]);
	    h1[ig+it*ng]=c*cimagf(adcal[ig+it*ng]);
	  }
  }
  hilbert_trans(h1, h2);

  for(it=0; it<nt; it++){
	  for(ig=0; ig<ng; ig++){
	    a=fabsf(adobs[ig+it*ng])+SF_EPS;
	    b=fabsf(adcal[ig+it*ng])+SF_EPS;
	    c=logf(a/b)/(b*b);
	    d=c*crealf(adcal[ig+it*ng]);
	    derr[ig+it*ng]=-d+cimagf(h2[ig+it*ng]);
	  }
  }

  free(h1);
  free(h2);
}

int main(int argc, char *argv[])
{
	/* variables on host */
	bool verb, precon, csdgather;
	int is, it, iter, niter, distx, distz, csd, rbell;
	int nz, nx, ns, ng, nt;
	int sxbeg, szbeg, gxbeg, gzbeg, jsx, jsz, jgx, jgz;/*  parameters of acquisition geometery */
	float dx, dz, fm, dt, dtx, dtz, tmp, amp, obj1, obj, beta, epsil, alpha;
	float *dobs, *dcal, *derr, *wlt, *bndr, *trans, *objval;
	int *sxz, *gxz;		
	float **vv, **illum, **lap, **vtmp, **sp0, **sp1, **sp2, **gp0, **gp1, **gp2, **g0, **g1, **cg, *alpha1, *alpha2, **ptr=NULL;
	sf_complex *adobs, *adcal;
	clock_t start, stop;/* timer */
	sf_file vinit, shots, vupdates, grads, objs, illums;/* I/O files */

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/* set up I/O files */
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
	shots=sf_input("shots"); /* recorded shots from exact velocity model */
    	vupdates=sf_output("out"); /* updated velocity in iterations */ 
    	grads=sf_output("grads");  /* gradient in iterations */ 
	illums=sf_output("illums");/* source illumination in iterations */
	objs=sf_output("objs");/* values of objective function in iterations */

    	/* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;/* vebosity */
    	if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");/* nz */
    	if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");/* nx */
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");/* dz */
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");/* dx */
	if (!sf_getbool("precon",&precon)) precon=false;/* precondition or not */
    	if (!sf_getint("niter",&niter))   niter=100;	/* number of iterations */
	if (!sf_getint("rbell",&rbell))	  rbell=2;	/* radius of bell smooth */

   	if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
	/* time sampling interval */
   	if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
	/* maximum amplitude of ricker */
   	if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
	/* dominant freq of ricker */
   	if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
	/* source x-axis  jump interval  */
   	if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
	/* source z-axis jump interval  */
   	if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
	/* receiver x-axis jump interval  */
   	if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
	/* receiver z-axis jump interval  */
   	if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
	/* default, common shot-gather; if n, record at every point*/

	sf_putint(vupdates,"n1",nz);	
	sf_putint(vupdates,"n2",nx);
	sf_putfloat(vupdates,"d1",dz);
	sf_putfloat(vupdates,"d2",dx);
	sf_putstring(vupdates,"label1","Depth");
	sf_putstring(vupdates,"label2","Distance");
	sf_putstring(vupdates,"label3","Iteration");
	sf_putint(vupdates,"n3",niter);
	sf_putint(vupdates,"d3",1);
	sf_putint(vupdates,"o3",1);
	sf_putint(grads,"n1",nz);	
	sf_putint(grads,"n2",nx);
	sf_putint(grads,"n3",niter);
	sf_putfloat(grads,"d1",dz);
	sf_putfloat(grads,"d2",dx);
	sf_putint(grads,"d3",1);
	sf_putint(grads,"o3",1);
	sf_putstring(grads,"label1","Depth");
	sf_putstring(grads,"label2","Distance");
	sf_putstring(grads,"label3","Iteration");
	sf_putint(illums,"n1",nz);	
	sf_putint(illums,"n2",nx);
	sf_putfloat(illums,"d1",dz);
	sf_putfloat(illums,"d2",dx);
	sf_putint(illums,"n3",niter);
	sf_putint(illums,"d3",1);
	sf_putint(illums,"o3",1);
	sf_putint(objs,"n1",niter);
	sf_putint(objs,"n2",1);
	sf_putfloat(objs,"d1",1);
	sf_putfloat(objs,"o1",1);

	dtx=dt/dx; 
	dtz=dt/dz; 
	csdgather=(csd>0)?true:false;

	vv=sf_floatalloc2(nz, nx);/* updated velocity */
	vtmp=sf_floatalloc2(nz, nx);/* temporary velocity computed with epsil */
	sp0=sf_floatalloc2(nz, nx);/* source wavefield p0 */
	sp1=sf_floatalloc2(nz, nx);/* source wavefield p1 */
	sp2=sf_floatalloc2(nz, nx);/* source wavefield p2 */
	gp0=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p0 */
	gp1=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p1 */
	gp2=sf_floatalloc2(nz, nx);/* geophone/receiver wavefield p2 */
	g0=sf_floatalloc2(nz, nx);/* gradient at previous step */
	g1=sf_floatalloc2(nz, nx);/* gradient at curret step */
	cg=sf_floatalloc2(nz, nx);/* conjugate gradient */
	lap=sf_floatalloc2(nz, nx);/* laplace of the source wavefield */
	illum=sf_floatalloc2(nz, nx);/* illumination of the source wavefield */
	objval=(float*)malloc(niter*sizeof(float));/* objective/misfit function */
	wlt=(float*)malloc(nt*sizeof(float));/* ricker wavelet */
	sxz=(int*)malloc(ns*sizeof(int)); /* source positions */
	gxz=(int*)malloc(ng*sizeof(int)); /* geophone positions */
	bndr=(float*)malloc(nt*(2*nz+nx)*sizeof(float));/* boundaries for wavefield reconstruction */
	trans=(float*)malloc(ng*nt*sizeof(float));/* transposed one shot */
	dobs=(float*)malloc(ng*nt*sizeof(float));/* observed seismic data */
	dcal=(float*)malloc(ng*nt*sizeof(float));/* calculated/synthetic seismic data */
	adobs=(sf_complex*)malloc(ng*nt*sizeof(sf_complex));/* analytic observed seismic data */
	adcal=(sf_complex*)malloc(ng*nt*sizeof(sf_complex));/* analytic calculated/synthetic seismic data */
	derr=(float*)malloc(ns*ng*nt*sizeof(float));/* residual/error between synthetic and observation */
	alpha1=(float*)malloc(ng*sizeof(float));/* numerator of alpha, length=ng */
	alpha2=(float*)malloc(ng*sizeof(float));/* denominator of alpha, length=ng */

	/* initialize varibles */
	sf_floatread(vv[0], nz*nx, vinit);
	memset(sp0[0], 0, nz*nx*sizeof(float));
	memset(sp1[0], 0, nz*nx*sizeof(float));
	memset(sp2[0], 0, nz*nx*sizeof(float));
	memset(gp0[0], 0, nz*nx*sizeof(float));
	memset(gp1[0], 0, nz*nx*sizeof(float));
	memset(gp2[0], 0, nz*nx*sizeof(float));
	memset(g0[0], 0, nz*nx*sizeof(float));
	memset(g1[0], 0, nz*nx*sizeof(float));
	memset(cg[0], 0, nz*nx*sizeof(float));
	memset(lap[0], 0, nz*nx*sizeof(float));
	memset(vtmp[0], 0, nz*nx*sizeof(float));
	memset(illum[0], 0, nz*nx*sizeof(float));
	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
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
	memset(dcal, 0, ng*nt*sizeof(float));
	memset(adobs, 0, ng*nt*sizeof(sf_complex));
	memset(adcal, 0, ng*nt*sizeof(sf_complex));
	memset(derr, 0, ns*ng*nt*sizeof(float));
	memset(alpha1, 0, ng*sizeof(float));
	memset(alpha2, 0, ng*sizeof(float));
	memset(dobs, 0, ng*nt*sizeof(float));	
	memset(objval, 0, niter*sizeof(float));

/* 
	// test the correctness of hilbert transform
	float a[]={1,2,3,4};
	sf_complex b[4];
	hilbert_init(1, 4);
	hilbert_trans(a,b);
	for(int ii=0; ii<4; ii++) sf_warning("real b[%d]=%g",ii,crealf(b[ii]));
	hilbert_close();
*/

	hilbert_init(ng, nt);
	for(iter=0; iter<niter; iter++)
	{
		obj=0;
		if(verb) start=clock();// record starting time
		sf_seek(shots, 0L, SEEK_SET);
		memcpy(g0[0], g1[0], nz*nx*sizeof(float));
		memset(g1[0], 0, nz*nx*sizeof(float));
		memset(illum[0], 0, nz*nx*sizeof(float));
		for(is=0;is<ns;is++)
		{
			sf_floatread(trans, ng*nt, shots);
			matrix_transpose(trans, dobs, nt, ng);
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
				rw_bndr(&bndr[it*(2*nz+nx)], sp0, nz, nx, true);

				record_seis(&dcal[it*ng], gxz, sp0, ng, nz);				
			}
			hilbert_trans(dcal, adcal);
			hilbert_trans(dobs, adobs);
			obj+=cal_objective(adcal, adobs, ng, nt);
			cal_adjsource(adcal, adobs, &derr[is*ng*nt], ng, nt);

			ptr=sp0; sp0=sp1; sp1=ptr;
			memset(gp0[0], 0, nz*nx*sizeof(float));
			memset(gp1[0], 0, nz*nx*sizeof(float));
			for(it=nt-1; it>-1; it--)
			{
				rw_bndr(&bndr[it*(2*nz+nx)], sp1, nz, nx, false);
				step_backward(illum, lap, sp0, sp1, sp2, vv, dtz, dtx, nz, nx);
				add_source(sp1, &wlt[it], &sxz[is], 1, nz, false);

				add_source(gp1, &derr[is*ng*nt+it*ng], gxz, ng, nz, true);
				step_forward(gp0, gp1, gp2, vv, dtz, dtx, nz, nx);

				cal_gradient(g1, lap, gp1, nz, nx);
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;
				ptr=gp0; gp0=gp1; gp1=gp2; gp2=ptr;
			}
		}
		

		scale_gradient(g1, vv, illum, nz, nx, precon);		
		sf_floatwrite(illum[0], nz*nx, illums);
		bell_smoothz(g1, illum, rbell, nz, nx);
		bell_smoothx(illum, g1, rbell, nz, nx);
		sf_floatwrite(g1[0], nz*nx, grads);


		if (iter>0) beta=cal_beta(g0, g1, cg, nz, nx); else beta=0.0;
		cal_conjgrad(g1, cg, beta, nz, nx);
		epsil=cal_epsilon(vv, cg, nz, nx);

		sf_seek(shots, 0L, SEEK_SET);
		memset(alpha1, 0, ng*sizeof(float));
		memset(alpha2, 0, ng*sizeof(float));
		cal_vtmp(vtmp, vv, cg, epsil, nz, nx);
		for(is=0;is<ns;is++)
		{
			sf_floatread(trans, ng*nt, shots);
			matrix_transpose(trans, dobs, nt, ng);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
			}
			memset(sp0[0], 0, nz*nx*sizeof(float));
			memset(sp1[0], 0, nz*nx*sizeof(float));
			for(it=0; it<nt; it++)
			{
				add_source(sp1, &wlt[it], &sxz[is], 1, nz, true);			
				step_forward(sp0, sp1, sp2, vtmp, dtz, dtx, nz, nx);
				ptr=sp0; sp0=sp1; sp1=sp2; sp2=ptr;

				record_seis(&dcal[it*ng], gxz, sp0, ng, nz);
				sum_alpha12(alpha1, alpha2, &dcal[it*ng], &dobs[it*ng], &derr[is*ng*nt+it*ng], ng);
			}
		}
		alpha=cal_alpha(alpha1, alpha2, epsil, ng);
		update_vel(vv, cg, alpha, nz, nx);
		sf_floatwrite(vv[0], nz*nx, vupdates);

		if(iter==0) {obj1=obj; objval[iter]=1.0;}
		else	objval[iter]=obj/obj1;

		if(verb) {// output important information at each FWI iteration 
			sf_warning("obj=%f  beta=%f  epsil=%f  alpha=%f", obj, beta, epsil, alpha);
			stop=clock();// record ending time 
			sf_warning("iteration %d finished: %f (s)",iter+1, ((float)(stop-start))/CLOCKS_PER_SEC);
		}
	}
	sf_floatwrite(objval, niter, objs);
	hilbert_close();

	free(*vv); free(vv);
	free(*vtmp); free(vtmp);
	free(*sp0); free(sp0);
	free(*sp1); free(sp1);
	free(*sp2); free(sp2);
	free(*gp0); free(gp0);
	free(*gp1); free(gp1);
	free(*gp2); free(gp2);
	free(*g0); free(g0);
	free(*g1); free(g1);
	free(*cg); free(cg);
	free(*lap); free(lap);
	free(*illum); free(illum);
	free(objval);
	free(wlt);
	free(sxz);
	free(gxz);
	free(bndr);
	free(trans);
	free(dobs);
	free(dcal);
	free(adobs);
	free(adcal);
	free(derr);
	free(alpha1);
	free(alpha2);

	exit(0);
}
