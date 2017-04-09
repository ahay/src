/* RTM and angle gather (ADCIG) extraction using poynting vector
SPML boundary condition combined with 4-th order finite difference,
effective boundary saving strategy used!
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University (Pengliang Yang)

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

  Reference: 
    [1] Yoon, Kwangjin, and Kurt J. Marfurt. "Reverse-time migration 
	using the Poynting vector." Exploration Geophysics 37.1 (2006): 
	102-107.
    [2] Costa, J. C., et al. "Obliquity-correction imaging condition 
	for reverse time migration." Geophysics 74.3 (2009): S57-S66.
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static bool 	csdgather;/* common shot gather (CSD) or not */
static int 	nb, nz, nx, nzpad, nxpad, nt, ns, ng, na;
static float 	fm, dt, dz, dx, _dz, _dx, da, var;

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

void  pmlcoeff_init(float *d1z, float *d2x, float vmax)
/*< initialize PML abosorbing coefficients >*/
{
	int ix, iz;
 	float Rc=1.e-5;
    	float x, z, L=nb* SF_MAX(dx,dz);
    	float d0=-3.*vmax*logf(Rc)/(2.*L*L*L);

	for(ix=0; ix<nxpad; ix++)
	{
		x=0;
	    	if (ix>=0 && ix<nb)   x=(ix-nb)*dx;
	    	else if(ix>=nxpad-nb && ix<nxpad) x=(ix-(nxpad-nb-1))*dx;
	    	d2x[ix] = d0*x*x;
	}
	for(iz=0; iz<nzpad; iz++)
	{
		z=0;
	    	if (iz>=0 && iz<nb)   z=(iz-nb)*dz;
	    	else if(iz>=nzpad-nb && iz<nzpad) z=(iz-(nzpad-nb-1))*dz; 
	    	d1z[iz] = d0*z*z;  
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

void wavefield_init(float **p, float **pz, float **px, float **vx, float **vz)
/*< wavefield initialization >*/
{
	memset(p[0],0,nzpad*nxpad*sizeof(float));
	memset(pz[0],0,nzpad*nxpad*sizeof(float));
	memset(px[0],0,nzpad*nxpad*sizeof(float));
	memset(vz[0],0,nzpad*nxpad*sizeof(float));
	memset(vx[0],0,nzpad*nxpad*sizeof(float));
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
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]+=source[is];
		}
	}else{ /* subtract sources */
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(is,sx,sz)		\
	shared(p,source,sxz,nb,ns,nz)
#endif
		for(is=0;is<ns; is++){
			sx=sxz[is]/nz+nb;
			sz=sxz[is]%nz+nb;
			p[sx][sz]-=source[is];
		}
	}
}


void step_forward(float **p, float **pz, float **px, float **vz, float **vx, float **vv, float *d1z, float *d2x)
/*< forward modeling step >*/
{
	int i1, i2;
	float tmp, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1,diff2)			\
	shared(nzpad, nxpad, p, vz, vx, dt, _dz, _dx, d1z, d2x)
#endif
	for(i2=1; i2<nxpad-2; i2++)
	for(i1=1; i1<nzpad-2; i1++)
	{
		diff1=1.125*(p[i2][i1+1]-p[i2][i1])-0.041666666666667*(p[i2][i1+2]-p[i2][i1-1]);
		diff2=1.125*(p[i2+1][i1]-p[i2][i1])-0.041666666666667*(p[i2+2][i1]-p[i2-1][i1]);
		vz[i2][i1]=((1.-0.5*dt*d1z[i1])*vz[i2][i1]+dt*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
		vx[i2][i1]=((1.-0.5*dt*d2x[i2])*vx[i2][i1]+dt*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1, diff2, tmp)		\
	shared(nzpad, nxpad, vv, p, pz, px, vz, vx, dt, _dz, _dx, d1z, d2x)
#endif
	for(i2=2; i2<nxpad-1; i2++)
	for(i1=2; i1<nzpad-1; i1++)
	{
		tmp=vv[i2][i1]; tmp=tmp*tmp;
		diff1=1.125*(vz[i2][i1]-vz[i2][i1-1])-0.041666666666667*(vz[i2][i1+1]-vz[i2][i1-2]);
		diff2=1.125*(vx[i2][i1]-vx[i2-1][i1])-0.041666666666667*(vx[i2+1][i1]-vx[i2-2][i1]);
		pz[i2][i1]=((1.-0.5*dt*d1z[i1])*pz[i2][i1]+dt*tmp*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
		px[i2][i1]=((1.-0.5*dt*d2x[i2])*px[i2][i1]+dt*tmp*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
		p[i2][i1]=px[i2][i1]+pz[i2][i1];
	}
}


void step_backward(float **p, float **vz, float **vx, float **vv)
/*< backward reconstruction step >*/
{
	int i1, i2;
	float tmp, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1, diff2, tmp)		\
	shared(nz, nx, nb, vv, p, vz, vx, dt, _dz, _dx)
#endif
	for(i2=nb; i2<nx+nb; i2++)
	for(i1=nb; i1<nz+nb; i1++)
	{
		tmp=vv[i2][i1]; tmp=tmp*tmp;
		diff1=1.125*(vz[i2][i1]-vz[i2][i1-1])-0.041666666666667*(vz[i2][i1+1]-vz[i2][i1-2]);
		diff2=1.125*(vx[i2][i1]-vx[i2-1][i1])-0.041666666666667*(vx[i2+1][i1]-vx[i2-2][i1]);
		tmp=tmp*(_dz*diff1+_dx*diff2);
		p[i2][i1]-=dt*tmp;
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1,diff2)			\
	shared(nz, nx, nb, p, vz, vx, dt, _dz, _dx)
#endif
	for(i2=nb; i2<nx+nb; i2++)
	for(i1=nb; i1<nz+nb; i1++)
	{
		diff1=1.125*(p[i2][i1+1]-p[i2][i1])-0.041666666666667*(p[i2][i1+2]-p[i2][i1-1]);
		diff2=1.125*(p[i2+1][i1]-p[i2][i1])-0.041666666666667*(p[i2+2][i1]-p[i2-1][i1]);
		diff1*=_dz;
		diff2*=_dx;
		vz[i2][i1]-=dt*diff1;
		vx[i2][i1]-=dt*diff2;
	}
}

void bndr_rw(bool read, float **vz, float **vx, float *bndr)
/*< read boundaries into vx and vz or write vx and vz for saving >*/
{
	int i1, i2;

	if(read){
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1,i2)			\
	shared(nz,nx,nb,vz,bndr)
#endif
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<4; i1++)
		{	
			vz[i2+nb][i1+nb-2]=bndr[i1+8*i2];
			vz[i2+nb][i1+nz+nb-2]=bndr[i1+4+8*i2];
		}
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1,i2)			\
	shared(nz,nx,nb,vx,bndr)
#endif
		for(i2=0; i2<4; i2++)
		for(i1=0; i1<nz; i1++)
		{
			vx[i2+nb-2][i1+nb]=bndr[8*nx+i1+nz*i2];
			vx[i2+nx+nb-2][i1+nb]=bndr[8*nx+i1+nz*(i2+4)];
		}
	}else{	
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1,i2)			\
	shared(nz,nx,nb,vz,bndr)
#endif
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<4; i1++)
		{	
			bndr[i1+8*i2]=vz[i2+nb][i1+nb-2];
			bndr[i1+4+8*i2]=vz[i2+nb][i1+nz+nb-2];
		}
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1,i2)			\
	shared(nz,nx,nb,vx,bndr)
#endif
		for(i2=0; i2<4; i2++)
		for(i1=0; i1<nz; i1++)
		{
			bndr[8*nx+i1+nz*i2]=vx[i2+nb-2][i1+nb];
			bndr[8*nx+i1+nz*(i2+4)]=vx[i2+nx+nb-2][i1+nb];
		}
	}
}


void record_seis(float *seis_it, int *gxz, float **p, int ng)
/*< record seismogram at time it into a vector length of ng >*/
{
	int ig, gx, gz;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(ig,gx,gz)		\
	shared(seis_it,p,gxz,nb,ng,nz)
#endif
	for(ig=0;ig<ng; ig++)
	{
		gx=gxz[ig]/nz+nb;
		gz=gxz[ig]%nz+nb;
		seis_it[ig]=p[gx][gz];
	}
}


void matrix_transpose(float **mat, int n1, int n2)
/*< transpose a matrix n1(row) x n2(column) into n2xn1 >*/
{
	int i1, i2;
	float **tmp;
	tmp=sf_floatalloc2(n2, n1);

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		tmp[i1][i2]=mat[i2][i1];
	}
	memcpy(mat[0], tmp[0], n1*n2*sizeof(float));

	free(*tmp);free(tmp);
}

void cross_correlation(float ***num, float **den, float **sp, float **gp, float **svz, float **svx, float **gvz, float **gvx)
/*< cross correlation and poynting vector computing >*/
{
	int i1, i2, ia;
	float Ssz,Ssx,Sgz,Sgx,b1, b2, a,tmp;

	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
		Ssz=sp[i2+nb][i1+nb]*svz[i2+nb][i1+nb];
		Ssx=sp[i2+nb][i1+nb]*svx[i2+nb][i1+nb];
		Sgz=gp[i2+nb][i1+nb]*gvz[i2+nb][i1+nb];
		Sgx=gp[i2+nb][i1+nb]*gvx[i2+nb][i1+nb];
		b1=Ssz*Ssz+Ssx*Ssx;//|Ss|^2
		b2=Sgz*Sgz+Sgx*Sgx;//|Sg|^2
		a=Ssx*Sgx+Ssz*Sgz; //<Ss,Sg>
		a=a/sqrtf(b1*b2);//+SF_EPS);	
	
		a=0.5*acosf(a);
		ia=(int)(a/da);
		if(ia<0) ia=0;
		if(ia==na) ia=ia-1;
		tmp=sp[i2+nb][i1+nb]*gp[i2+nb][i1+nb];//numerator 
		//tmp=tmp*expf(-(a-ia*da)*(a-ia*da)/var); //Gaussian smoothing in Fresnel zone
		num[ia][i2][i1]+=tmp;
		den[i2][i1]+=gp[i2+nb][i1+nb]*gp[i2+nb][i1+nb];//denominator
	}
}



int main(int argc, char* argv[])
{
  int it,kt,ia,is,ig,i1,i2,jsx,jsz,jgx,jgz,sxbeg,szbeg,gxbeg,gzbeg, distx, distz;
	int *sxz, *gxz;
	float tmp, amp, vmax;
	float *wlt, *d2x, *d1z, *bndr, *adjsource;
	float **v0, **vv, **vvs, **dcal, **dobs, **den;
	float **sp, **spz, **spx, **svz, **svx, **gp, **gpz, **gpx, **gvz, **gvx;
	float ***num, ***adcig;
    	sf_file vmodl,vmods, rtmadcig, vecx, vecz; /* I/O files */

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

    	/*< set up I/O files >*/
    	vmodl = sf_input ("in");   /* velocity model, unit=m/s */
	vmods = sf_input("velsmooth");/* smooth background velocity model for muting direct arrival*/
    	rtmadcig = sf_output("out");  /* ADCIG obtained by Poynting vector */
	vecx=sf_output("vecx");
	vecz=sf_output("vecz");

    	/* get parameters for RTM */
    	if (!sf_histint(vmodl,"n1",&nz)) sf_error("no n1");
    	if (!sf_histint(vmodl,"n2",&nx)) sf_error("no n2");
    	if (!sf_histfloat(vmodl,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vmodl,"d2",&dx)) sf_error("no d2");

    	if (!sf_getfloat("amp",&amp)) amp=1.e3;	
	/* maximum amplitude of ricker wavelet*/
    	if (!sf_getfloat("fm",&fm)) sf_error("no fm");	
	/* dominant freq of ricker */
    	if (!sf_getfloat("dt",&dt)) sf_error("no dt");	
	/* time interval */
    	if (!sf_getint("nt",&nt))   sf_error("no nt");	
	/* total modeling time steps */
    	if (!sf_getint("ns",&ns))   sf_error("no ns");	
	/* total shots */
    	if (!sf_getint("ng",&ng))   sf_error("no ng");	
	/* total receivers in each shot */
    	if (!sf_getint("nb",&nb))   nb=20; 
	/* thickness of split PML */
    	if (!sf_getint("na",&na)) na=30;
	/* number of angles*/
    	if (!sf_getint("kt",&kt))   kt=200;
	/* record poynting vector at kt */
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
	if (!sf_getbool("csdgather",&csdgather)) csdgather=true;
	/* default, common shot-gather; if n, record at every point*/

	_dx=1./dx;
	_dz=1./dz;
	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	da=SF_PI/(2.*na);/* angle unit, rad */
	var=da/3.;
	var=2.0*var*var;

    	sf_putint(rtmadcig,"n1",nz);
    	sf_putint(rtmadcig,"n2",nx);
	sf_putfloat(rtmadcig,"n3",na);
    	sf_putfloat(rtmadcig,"d1",dz);
    	sf_putfloat(rtmadcig,"d2",dx);
	sf_putfloat(rtmadcig,"d3",90./(float)na);

	/* allocate variables */
	wlt=sf_floatalloc(nt);
	v0=sf_floatalloc2(nz,nx); 	
	vv=sf_floatalloc2(nzpad, nxpad);
	vvs=sf_floatalloc2(nzpad, nxpad);
	sp =sf_floatalloc2(nzpad, nxpad);
	spz=sf_floatalloc2(nzpad, nxpad);
	spx=sf_floatalloc2(nzpad, nxpad);
	svz=sf_floatalloc2(nzpad, nxpad);
	svx=sf_floatalloc2(nzpad, nxpad);
	gp =sf_floatalloc2(nzpad, nxpad);
	gpz=sf_floatalloc2(nzpad, nxpad);
	gpx=sf_floatalloc2(nzpad, nxpad);
	gvz=sf_floatalloc2(nzpad, nxpad);
	gvx=sf_floatalloc2(nzpad, nxpad);
	d1z=sf_floatalloc(nzpad);
	d2x=sf_floatalloc(nxpad);
	sxz=sf_intalloc(ns);
	gxz=sf_intalloc(ng);
	dobs=sf_floatalloc2(ng,nt);
	dcal=sf_floatalloc2(ng,nt);
	adjsource=sf_floatalloc(ng);
	bndr=(float*)malloc(nt*8*(nx+nz)*sizeof(float));
	den=sf_floatalloc2(nz,nx);
	num=sf_floatalloc3(nz,nx,na);
	adcig=sf_floatalloc3(nz,nx,na);

	/* initialize variables */
	for(it=0;it<nt;it++){
		tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
		wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
	}
	/* time integration for true amplitude RTM */
	for(it=1;it<nt;it++) wlt[it]=wlt[it]+wlt[it-1];
	sf_floatread(v0[0],nz*nx,vmodl);
	expand2d(vv, v0);
	sf_floatread(v0[0],nz*nx,vmods);
	expand2d(vvs, v0);
	memset(sp [0],0,nzpad*nxpad*sizeof(float));
	memset(spx[0],0,nzpad*nxpad*sizeof(float));
	memset(spz[0],0,nzpad*nxpad*sizeof(float));
	memset(svx[0],0,nzpad*nxpad*sizeof(float));
	memset(svz[0],0,nzpad*nxpad*sizeof(float));
	memset(gp [0],0,nzpad*nxpad*sizeof(float));
	memset(gpx[0],0,nzpad*nxpad*sizeof(float));
	memset(gpz[0],0,nzpad*nxpad*sizeof(float));
	memset(gvx[0],0,nzpad*nxpad*sizeof(float));
	memset(gvz[0],0,nzpad*nxpad*sizeof(float));
	vmax=v0[0][0];
	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
		vmax=SF_MAX(v0[i2][i1],vmax);
	pmlcoeff_init(d1z, d2x, vmax);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_error("sources exceeds the computing zone!"); exit(1);}
	sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_error("geophones exceeds the computing zone!"); exit(1);}
	}else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ sf_error("geophones exceeds the computing zone!"); exit(1);}
	}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
	memset(adcig[0][0], 0, na*nz*nx*sizeof(float));

	for(is=0; is<ns; is++)
	{
	        sf_warning("source:%d",is+1);

		memset(dobs[0], 0, ng*nt*sizeof(float));
		memset(dcal[0], 0, ng*nt*sizeof(float));
		memset(adjsource, 0, ng*sizeof(float));
		wavefield_init(sp, spz, spx, svz, svx);
		wavefield_init(gp, gpz, gpx, gvz, gvx);
		if (csdgather)	{
			gxbeg=sxbeg+is*jsx-distx;
			sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
		}
		for(it=0; it<nt; it++)
		{
			add_source(&sxz[is], sp, 1, &wlt[it], true);
			step_forward(sp, spz, spx, svz, svx, vv, d1z, d2x);
			bndr_rw(false, svz, svx, &bndr[it*8*(nx+nz)]);		
			record_seis(dobs[it], gxz, sp, ng);

			//remodeling with a smooth background velocity model
			add_source(&sxz[is], gp, 1, &wlt[it], true);
			step_forward(gp, gpz, gpx, gvz, gvx, vvs, d1z, d2x);
			record_seis(dcal[it], gxz, gp, ng);
		}
		sf_warning("forward done!");

		wavefield_init(gp, gpz, gpx, gvz, gvx);
		memset(num[0][0], 0, na*nz*nx*sizeof(float));
		memset(den[0], 0, nz*nx*sizeof(float));
		for(it=nt-1; it>-1; it--)
		{	
                        //muting direct arrival
                        for(ig=0; ig<ng; ig++) adjsource[ig]=dobs[it][ig]-dcal[it][ig];

			add_source(gxz, gp, ng, adjsource, true);
			step_forward(gp, gpz, gpx, gvz, gvx, vv, d1z, d2x);	

			if(it==kt)
			{
				window2d(v0,svx);
				sf_floatwrite(v0[0],nz*nx,vecx);
				window2d(v0,svz);
				sf_floatwrite(v0[0],nz*nx,vecz);
			}

			bndr_rw(true, svz, svx, &bndr[it*8*(nx+nz)]);	
			cross_correlation(num, den, sp, gp, svz, svx, gvz, gvx);

			step_backward(sp, svz, svx, vv);
			add_source(&sxz[is], sp, 1, &wlt[it], false);
		}	
		sf_warning("backward done!");

		for(ia=0; ia<na; ia++)
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<nz; i1++)
			adcig[ia][i2][i1]+=num[ia][i2][i1]/(den[i2][i1]+SF_EPS);
	}
	sf_floatwrite(adcig[0][0], na*nz*nx,rtmadcig);

	free(wlt);
	free(*dobs); free(dobs);
	free(*dcal); free(dcal);
	free(*v0); free(v0);
	free(*vv); free(vv);
	free(*vvs); free(vvs);
	free(*sp); free(sp);
	free(*spx); free(spx);
	free(*spz); free(spz);
	free(*svx); free(svx);
	free(*svz); free(svz);
	free(*gp); free(gp);
	free(*gpx); free(gpx);
	free(*gpz); free(gpz);
	free(*gvx); free(gvx);
	free(*gvz); free(gvz);
	free(d1z);
	free(d2x);
	free(sxz);
	free(gxz);
	free(bndr);
	free(adjsource);
	free(*den); free(den);
	free(**num); free(*num); free(num);
	free(**adcig); free(*adcig); free(adcig);

    	exit(0);
}

