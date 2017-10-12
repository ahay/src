/* RTM output ODCIG with extended images
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
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static bool 	csdgather;/* common shot gather (CSD) or not */
static int 	nb, nz, nx, nzpad, nxpad, nt, ns, ng, nh;
static float 	fm, dt, dz, dx, _dz, _dx, vmute;

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
#pragma omp parallel for default(none)		\
  private(is,sx,sz)				\
  shared(p,source,sxz,nb,ns,nz)
#endif
    for(is=0;is<ns; is++){
      sx=sxz[is]/nz+nb;
      sz=sxz[is]%nz+nb;
      p[sx][sz]+=source[is];
    }
  }else{ /* subtract sources */
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  private(is,sx,sz)				\
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
#pragma omp parallel for default(none)				\
  private(i1,i2,diff1,diff2)					\
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
#pragma omp parallel for default(none)					\
  private(i1,i2,diff1, diff2, tmp)					\
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
  private(i1,i2,diff1, diff2, tmp)			\
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
#pragma omp parallel for default(none)		\
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
#pragma omp parallel for default(none)		\
  private(i1,i2)				\
  shared(nz,nx,nb,vz,bndr)
#endif
    for(i2=0; i2<nx; i2++)
      for(i1=0; i1<4; i1++)
	{	
	  vz[i2+nb][i1+nb-2]=bndr[i1+8*i2];
	  vz[i2+nb][i1+nz+nb-2]=bndr[i1+4+8*i2];
	}
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  private(i1,i2)				\
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
#pragma omp parallel for default(none)		\
  private(i1,i2)				\
  shared(nz,nx,nb,vz,bndr)
#endif
    for(i2=0; i2<nx; i2++)
      for(i1=0; i1<4; i1++)
	{	
	  bndr[i1+8*i2]=vz[i2+nb][i1+nb-2];
	  bndr[i1+4+8*i2]=vz[i2+nb][i1+nz+nb-2];
	}
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  private(i1,i2)				\
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
#pragma omp parallel for default(none)		\
  private(ig,gx,gz)				\
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

void muting(float *seis_kt, int gzbeg, int szbeg, int gxbeg, int sxc, int jgx, int it, int tdmute)
/*< muting the direct arrivals >*/
{
  int id, kt;
  float a,b,t0;

#ifdef _OPENMP
#pragma omp parallel for default(none)					\
  private(id,a,b,t0,kt)							\
  shared(seis_kt,vmute,ng,it,dt,dx,dz,tdmute,gzbeg,szbeg,gxbeg,sxc,jgx)
#endif
  for(id=0;id<ng;id++)
    {
      a=dx*abs(gxbeg+id*jgx-sxc);
      b=dz*(gzbeg-szbeg);
      t0=sqrtf(a*a+b*b)/vmute;
      kt=t0/dt+tdmute;// tdmute manually added to obtain the best muting effect.
      if (it<kt) seis_kt[id]=0.;
    }
}


void cross_correlation(float ***odcig, float **sp, float **gp, float **svz, float **svx, float **gvz, float **gvx)
/*< cross correlation and poynting vector computing >*/
{
  int i1, i2,ih;

  for(i2=0; i2<nx; i2++){
    for(i1=0; i1<nz; i1++){
      for(ih=-nh; ih<=nh;ih++){
	if(i2-ih>=0 && i2-ih<nx && i2+ih>=0 && i2+ih<nx) odcig[ih+nh][i2][i1]+=sp[i2-ih+nb][i1+nb]*gp[i2+ih+nb][i1+nb];
      }
    }
  }
}


int main(int argc, char* argv[])
{
  int it,kt,is,ig,i1,i2,tdmute,jsx,jsz,jgx,jgz,sxbeg,szbeg,gxbeg,gzbeg, distx, distz;
  int *sxz, *gxz;
  float tmp, amp, vmax;
  float *wlt, *d2x, *d1z, *bndr, *adjsource;
  float **v0, **vv, **vv1st,**dcal,**dobs,***odcig;
  float **sp, **spz, **spx, **svz, **svx;//source wavefield
  float **gp, **gpz, **gpx, **gvz, **gvx;//receiver wavefield

  sf_file vmodl, vmods, Fodcig;

  sf_init(argc,argv);
#ifdef _OPENMP
  omp_init();
#endif

  /*< set up I/O files >*/
  vmodl = sf_input("in");   /* velocity model, unit=m/s */
  vmods = sf_input("vel1stlayer"); /* velocity model at 1st layer*/
  Fodcig = sf_output("out");  /* ODCIG by subsurface offset*/

  /* get parameters for RTM */
  if (!sf_histint(vmodl,"n1",&nz)) sf_error("no n1");
  if (!sf_histint(vmodl,"n2",&nx)) sf_error("no n2");
  if (!sf_histfloat(vmodl,"d1",&dz)) sf_error("no d1");
  if (!sf_histfloat(vmodl,"d2",&dx)) sf_error("no d2");

  if (!sf_getfloat("amp",&amp)) amp=1.;	
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
  if (!sf_getint("nh",&nh)) nh=30;
  /* number of points in offset coordinate*/
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
  if (!sf_getfloat("vmute",&vmute))   vmute=1500;
  /* muting velocity to remove the low-freq noise, unit=m/s*/
  if (!sf_getint("tdmute",&tdmute))   tdmute=2./(fm*dt);
  /* number of deleyed time samples to mute */

  sf_putint(Fodcig,"n1",nz);
  sf_putint(Fodcig,"n2",nx);
  sf_putfloat(Fodcig,"n3",2*nh+1);
  sf_putfloat(Fodcig,"d1",dz);
  sf_putfloat(Fodcig,"d2",dx);
  sf_putfloat(Fodcig,"d3",dx);
  sf_putfloat(Fodcig,"o3",-nh*dx);

  _dx=1./dx;
  _dz=1./dz;
  nzpad=nz+2*nb;
  nxpad=nx+2*nb;

  /* allocate variables */
  wlt=sf_floatalloc(nt);
  v0=sf_floatalloc2(nz,nx); 	
  vv=sf_floatalloc2(nzpad, nxpad);
  vv1st=sf_floatalloc2(nzpad, nxpad);
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
  dcal=sf_floatalloc2(ng,nt);
  dobs=sf_floatalloc2(ng,nt);
  adjsource=sf_floatalloc(ng);
  bndr=(float*)malloc(nt*8*(nx+nz)*sizeof(float));
  odcig=sf_floatalloc3(nz,nx,2*nh+1);

  /* initialize variables */
  for(it=0;it<nt;it++){
    tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
    wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
  }
  sf_floatread(v0[0],nz*nx,vmodl);
  expand2d(vv, v0);//true input velocity
  sf_floatread(v0[0],nz*nx,vmods);
  expand2d(vv1st, v0);//velocity at source location/ 1st layer
  wavefield_init(sp, spz, spx, svz, svx);
  wavefield_init(gp, gpz, gpx, gvz, gvx);
  vmax=v0[0][0];
  for(i2=0; i2<nx; i2++)
    for(i1=0; i1<nz; i1++)
      vmax=SF_MAX(v0[i2][i1],vmax);
  pmlcoeff_init(d1z, d2x, vmax);
  if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))    {
    sf_error("sources exceeds the computing zone!"); exit(1);
  }
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
  memset(dobs[0],0,ng*nt*sizeof(float));
  memset(dcal[0],0,ng*nt*sizeof(float));
  memset(adjsource,0,ng*sizeof(float));
  memset(odcig[0][0], 0, (2*nh+1)*nz*nx*sizeof(float));

  for(is=0; is<ns; is++)    {
      sf_warning("shot %d",is);

      wavefield_init(sp, spz, spx, svz, svx);
      wavefield_init(gp, gpz, gpx, gvz, gvx);
      if (csdgather)	{
	gxbeg=sxbeg+is*jsx-distx;
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
      }
      for(it=0; it<nt; it++)	{
	//source wavefield simulation with true velocity: vv
	add_source(&sxz[is], sp, 1, &wlt[it], true);
	step_forward(sp, spz, spx, svz, svx, vv, d1z, d2x);
	bndr_rw(false, svz, svx, &bndr[it*8*(nx+nz)]);
	record_seis(dobs[it], gxz, sp, ng);

	//source wavefield simulation with velocity at 1st layer: vv1st
	add_source(&sxz[is], gp, 1, &wlt[it], true);
	step_forward(gp, gpz, gpx, gvz, gvx, vv1st, d1z, d2x);		
	record_seis(dcal[it], gxz, gp, ng);
      }

      //re-initialize the wavefield on receiver/geophone side
      wavefield_init(gp, gpz, gpx, gvz, gvx);
      for(it=nt-1; it>-1; it--)	{	
	//mute the direct arrival by subtracting the above two recorded seismograms
	for(ig=0; ig<ng; ig++) adjsource[ig]=dobs[it][ig]-dcal[it][ig];

	add_source(gxz, gp, ng, adjsource, true);
	step_forward(gp, gpz, gpx, gvz, gvx, vv, d1z, d2x);	

	bndr_rw(true, svz, svx, &bndr[it*8*(nx+nz)]);	
	cross_correlation(odcig, sp, gp, svz, svx, gvz, gvx);

	step_backward(sp, svz, svx, vv);
	add_source(&sxz[is], sp, 1, &wlt[it], false);
      }	
    }
  sf_floatwrite(odcig[0][0], (2*nh+1)*nz*nx,Fodcig);

  free(wlt);
  free(*dobs); free(dobs);
  free(*dcal); free(dcal);
  free(*v0); free(v0);
  free(*vv); free(vv);
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
  free(**odcig); free(*odcig); free(odcig);

  exit(0);
}

