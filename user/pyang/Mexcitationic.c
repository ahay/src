/* Demo for excitation imaging condition
Note that excitation imaging condition has some multipathing artifacts.
The significance of this imaging condition is the cheap computation and
low memory requirement. (1) Cheap computation: only 2 times of simulation,
one for source wavefield the other for receiver wavefield, are needed for 
single shot imaging before stacking. (2) Low memory request: this imaging 
condition only asks for the excitation time and the amplitude. Therefore,
it differs from cross-correlation imaging condition which needs storing 
or reconstructing the source wavefield to cross-correlate with receiver
wavefield at each time step.
 */
/*
  Copyright (C) 2016  Univ. Grenoble Alpes, Pengliang Yang

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
  [1]  Wen‐Fong Chang and George A. McMechan (1986). ”Reverse‐time migration
       of offset vertical seismic profiling data using the excitation‐time 
       imaging condition.” GEOPHYSICS, 51(1), 67-84. doi: 10.1190/1.1442041
  [2] Bao D. Nguyen, George A. McMechan. (2013) Excitation amplitude imaging
      condition for prestack reverse-time migration. GEOPHYSICS 78:1, S37-S46.
*/

#include <rsf.h>
#include <time.h>

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


void step_forward(float **p0, float **p1)
/*< one step of forward modeling >*/
{
  int ix,iz;
  float tmp;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
  private(ix,iz,tmp)					\
  shared(vv,p1,p0,nxpad,nzpad,c0,c11,c12,c21,c22)  
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

void apply_sponge(float **p0)
/* apply absorbing boundary condition */
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
	  p[ix+nb][iz-2+nb]=spo[iz+4*ix];
	  p[ix+nb][iz+nz+nb]=spo[iz+2+4*ix];
	}
#ifdef _OPENMP
#pragma omp parallel for			\
  private(ix,iz)				\
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
#pragma omp parallel for			\
  private(ix,iz)				\
  shared(p,spo,nx,nz,nb)  
#endif	
    for(ix=0; ix<nx; ix++)
      for(iz=0; iz<2; iz++)
	{
	  spo[iz+4*ix]=p[ix+nb][iz-2+nb];
	  spo[iz+2+4*ix]=p[ix+nb][iz+nz+nb];
	}
#ifdef _OPENMP
#pragma omp parallel for			\
  private(ix,iz)				\
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
/*< add source or subtract source >*/
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

//update excitation time and amplitude
void update_excitation_time(int it, int **ext, float **umax, float **p, int nz, int nx, int nb)
{
  int i1, i2;

  for(i2=0; i2<nx; i2++){
    for(i1=0; i1<nz; i1++){
      if(umax[i2][i1]< fabsf(p[i2+nb][i1+nb]) ) {
	umax[i2][i1]=fabsf(p[i2+nb][i1+nb]);
	ext[i2][i1]=it;
      }
    }
  }
}

void imaging_condition(int it, int **ext, float **umax, float **p0, float **img, int nz, int nx, int nb)
{
  int i1,i2;

  for(i2=0; i2<nx; i2++){
    for(i1=0; i1<nz; i1++){
      if(ext[i2][i1]==it) img[i2][i1]+=p0[i2+nb][i1+nb]/umax[i2][i1];	
    }
  }
}


int main(int argc, char* argv[])
{
  bool csdgather;
  int it, kt, ns, ng, iz, ix, ib, is;
  int jsx,jsz,sxbeg,szbeg,jgx,jgz,gxbeg,gzbeg,distx,distz;
  int *sxz, *gxz, **ext; //excitation time index
  float tmp;
  float *wlt, **v0, **dcal, **umax, **img;
  sf_file Fv, Fout;

  sf_init(argc,argv);
#ifdef _OPENMP
  omp_init();
#endif

  Fv = sf_input("in");/* veloctiy model */
  Fout = sf_output("out");/* excitation amplitude image */

  if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
  /* veloctiy model: nz */
  if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
  /* veloctiy model: nx */
  if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
  /* veloctiy model: dz */
  if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
  /* veloctiy model: dx */
  if (!sf_getint("nb",&nb)) nb=20;
  /* thickness of sponge ABC */
  if (!sf_getint("nt",&nt)) sf_error("nt required");
  /* number of time steps */
  if (!sf_getfloat("dt",&dt)) sf_error("dt required");
  /* time sampling interval */
  if (!sf_getfloat("fm",&fm)) fm=20.0;
  /*dominant freq of Ricker wavelet */
  if (!sf_getint("kt",&kt)) kt=300;	
  /* output wavefield at time kt */
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
  if (!sf_getint("ns",&ns)) ns=1;	/* number of shots */
  if (!sf_getint("ng",&ng)) ng=nx;/* number of receivers */

  /*< initialize 4-th order fd coefficients >*/
  tmp = 1.0/(dz*dz);
  c11 = 4.0*tmp/3.0;
  c12= -tmp/12.0;
  tmp = 1.0/(dx*dx);
  c21 = 4.0*tmp/3.0;
  c22= -tmp/12.0;
  c0=-2.0*(c11+c12+c21+c22);

  nzpad=nz+2*nb;
  nxpad=nx+2*nb;

  v0=sf_floatalloc2(nz,nx); 
  vv=sf_floatalloc2(nzpad, nxpad);
  bndr=sf_floatalloc(nb);
  wlt=(float*)malloc(nt*sizeof(float));
  sxz=(int*)malloc(ns*sizeof(int));
  gxz=(int*)malloc(ng*sizeof(int));	
  p0=sf_floatalloc2(nzpad, nxpad);
  p1=sf_floatalloc2(nzpad, nxpad);
  dcal=sf_floatalloc2(ng,nt);
  ext=sf_intalloc2(nz,nx);//excitation time index
  umax=sf_floatalloc2(nz,nx);
  img=sf_floatalloc2(nz,nx);

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
  for(it=0;it<nt;it++){
    tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
    wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
  }	
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
  memset(img[0],0,nz*nx*sizeof(float));

  for(is=0; is<ns; is++){
    // determine receiver locations 
    if (csdgather)	{
      gxbeg=sxbeg+is*jsx-distx;
      sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
    }
    //source wavefield simulation
    memset(p0[0],0,nzpad*nxpad*sizeof(float));
    memset(p1[0],0,nzpad*nxpad*sizeof(float));
    memset(dcal[0],0,nt*ng*sizeof(float));
    memset(umax[0],0,nz*nx*sizeof(float));
    memset(ext[0],0,nz*nx*sizeof(int));

    for(it=0; it<nt; it++)	{
      add_source(sxz, p1, 1, &wlt[it], true);
      step_forward(p0, p1);
      apply_sponge(p0);
      apply_sponge(p1);
      ptr=p0; p0=p1; p1=ptr;
      record_seis(dcal[it], gxz, p1, ng);

      //update excitation time and amplitude
      update_excitation_time(it, ext, umax, p1, nz, nx,nb);
    }

    //receiver wavefield simulation
    memset(p0[0],0,nzpad*nxpad*sizeof(float));
    memset(p1[0],0,nzpad*nxpad*sizeof(float));
    for(it=nt-1; it>-1; it--){
      add_source(gxz, p1, ng, dcal[it], true);
      step_forward(p0, p1);
      apply_sponge(p0);
      apply_sponge(p1);
      ptr=p0; p0=p1; p1=ptr;

      //apply excitation amplitude imaging condition
      imaging_condition(it, ext, umax, p1, img, nz, nx, nb);
    }
  }
  sf_floatwrite(img[0],nz*nx,Fout);


  free(sxz);
  free(gxz);
  free(wlt);
  free(*v0); free(v0);
  free(*vv); free(vv);
  free(*p0); free(p0);
  free(*p1); free(p1);
  free(bndr);
  free(*dcal); free(dcal);
  free(*ext); free(ext);
  free(*img); free(img);
  free(*umax); free(umax);

  exit(0);
}

