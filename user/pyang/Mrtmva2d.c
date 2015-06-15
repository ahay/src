/* RTM with checkpointing in 2D acoustic media
  The real value of checkpointing technology resides in the backpropagation with
  viscoacoustic and viscoelastic wave equation, where the wavefield 
  reconstruction method using saved boundaries fails. Here, we only
  demonstrate how to implement it in acoustic media without dissipation.
*/
/*
  Copyright (C) 2015 Xi'an Jiaotong University (Pengliang Yang)

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
  
  Reference: William Symes, Reverse time migration with optimal checkpointing,
  Geophysics, v. 72 no. 5 p. SM213-SM221 doi: 10.1190/1.2742686 
*/
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int nb, nz, nx, nzpad, nxpad, nt, ns, ng;
static float dz, dx, _dz, _dx, dt, fm;


typedef struct {
  int no;// number/index of the checkpoint in time axis
  float **p;
  float **r;
  float **vx;
  float **vz;
} checkpoint;


static int nb, nz, nx, nzpad, nxpad, nt, ns;
static float dz, dx, _dz, _dx, dt, fm;

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

void apply_sponge(float **u, float *bndr)
/*< apply absorbing boundary condition >*/
{
	int ix,iz;

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,u)
#endif
	for(ix=0; ix<nxpad; ix++)
	{
		for(iz=0;iz<nb;iz++){	// top ABC			
			u[ix][iz]=bndr[iz]*u[ix][iz];
		}
		for(iz=nz+nb;iz<nzpad;iz++){// bottom ABC			
			u[ix][iz]=bndr[nzpad-iz-1]*u[ix][iz];
		}
	}

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(ix,iz)		    \
    shared(bndr,u)
#endif
	for(iz=0; iz<nzpad; iz++)
	{
		for(ix=0;ix<nb;ix++){	// left ABC			
			u[ix][iz]=bndr[ix]*u[ix][iz];
		}	
		for(ix=nx+nb;ix<nxpad;ix++){// right ABC			
			u[ix][iz]=bndr[nxpad-ix-1]*u[ix][iz];
		}	
	}
}


void variable_inverse(float **rho, float **tauo)
/*< inverse of variables >*/
{
	int i1, i2;
	
	for(i2=0; i2<nxpad; i2++)
	for(i1=0; i1<nzpad; i1++)
	{
		rho[i2][i1]=1./rho[i2][i1];
		tauo[i2][i1]=1./tauo[i2][i1];
	}
}

void average_variable(float **rho, float **tau)
/*< average the parameters >*/
{

}

void step_forward(float **p, float **r, float **vz, float **vx, float **vv, float **rho, float **tau, float **tauo)
/*< forward modeling step >*/
{
	int i1, i2;
	float tmp, tmp2, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1,diff2)			\
	shared(nzpad, nxpad, rho, p, vz, vx, dt, _dz, _dx)
#endif
	for(i2=3; i2<nxpad-4; i2++)
	for(i1=3; i1<nzpad-4; i1++)
	{
		diff1=	 1.196289062500000*(p[i2][i1+1]-p[i2][i1])
			-0.079752604166667*(p[i2][i1+2]-p[i2][i1-1])
			+0.009570312500000*(p[i2][i1+3]-p[i2][i1-2])
			-0.000697544642857*(p[i2][i1+4]-p[i2][i1-3]);
		diff2=	 1.196289062500000*(p[i2+1][i1]-p[i2][i1])
			-0.079752604166667*(p[i2+2][i1]-p[i2-1][i1])
			+0.009570312500000*(p[i2+3][i1]-p[i2-2][i1])
			-0.000697544642857*(p[i2+4][i1]-p[i2-3][i1]);
		vz[i2][i1]-=dt*_dz*diff1/rho[i2][i1];
		vx[i2][i1]-=dt*_dx*diff2/rho[i2][i1];
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)			\
	private(i1,i2,diff1, diff2, tmp, tmp2)		\
	shared(nzpad, nxpad, rho, tau, tauo, vv, p, r, vz, vx, dt, _dz, _dx)
#endif
	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		tmp=vv[i2][i1]; tmp=tmp*tmp;
		diff1=	 1.196289062500000*(vz[i2][i1]-vz[i2][i1-1])
			-0.079752604166667*(vz[i2][i1+1]-vz[i2][i1-2])
			+0.009570312500000*(vz[i2][i1+2]-vz[i2][i1-3])
			-0.000697544642857*(vz[i2][i1+3]-vz[i2][i1-4]);
		diff2=	 1.196289062500000*(vx[i2][i1]-vx[i2-1][i1])
			-0.079752604166667*(vx[i2+1][i1]-vx[i2-2][i1])
			+0.009570312500000*(vx[i2+2][i1]-vx[i2-3][i1])
			-0.000697544642857*(vx[i2+3][i1]-vx[i2-4][i1]);
		tmp=tmp*rho[i2][i1]*(_dz*diff1+_dx*diff2);
		tmp2=dt/tauo[i2][i1];
		r[i2][i1]=((1.-0.5*tmp2)*r[i2][i1]-tmp2*tau[i2][i1]*tmp)/(1.+0.5*tmp2);
		p[i2][i1]-=dt*((1.+tau[i2][i1])*tmp+r[i2][i1]);
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


void muting(float *seis_kt, int gzbeg, int szbeg, int gxbeg, int sxc, int jgx, int it, int tdmute, float vmute)
/*< muting the direct arrivals >*/
{
  int id, kt;
  float a,b,t0;

  for(id=0;id<ng;id++)
    {
      a=dx*abs(gxbeg+id*jgx-sxc);
      b=dz*(gzbeg-szbeg);
      t0=sqrtf(a*a+b*b)/vmute;
      kt=t0/dt+tdmute;// tdmute manually added to obtain the best muting effect.
      if (it<kt) seis_kt[id]=0.;
    }
}


void buffer_init(int *checkpt,/*time index/location of checkpoints*/
		 checkpoint *buffer,/* buffer to store snapshots*/
		 int noc,/* number of checkpoints */
		 int nob,/* number of buffer points*/
		 int nt/* number of time steps*/)      
/*< initialize checkpoints and buffers
  Note that points in buffer are a subset of checkpoints.>*/
{
  int ib,ic;

  //set checkpoints in time axis
  for(ic=0; ic<noc; ic++)
    checkpt[ic]=nt*ic/noc;
  
  //set buffer points: allocate and initialize variables at each buffer point
  for(ib=0; ib<nob; ib++) {
    buffer[ib].p=sf_floatalloc2(nzpad, nxpad);
    buffer[ib].r=sf_floatalloc2(nzpad, nxpad);
    buffer[ib].vz=sf_floatalloc2(nzpad, nxpad);
    buffer[ib].vx=sf_floatalloc2(nzpad, nxpad);

    ic=(int)(noc-1)*ib/(nob-1);// index in checkpt[]   
    //ib==nob-1, ic==noc-1, the last buffer is the last checkpoint!
    buffer[ib].no=checkpt[ic];//set buffer to be one of checkpoints
    memset(buffer[ib].p[0], 0, nzpad*nxpad*sizeof(float));
    memset(buffer[ib].r[0], 0, nzpad*nxpad*sizeof(float));
    memset(buffer[ib].vz[0], 0, nzpad*nxpad*sizeof(float));
    memset(buffer[ib].vx[0], 0, nzpad*nxpad*sizeof(float));
  }
}

void buffer_remodeling(float **p,
		       float **r,
		       float **vz,
		       float **vx,
		       float **vv,
		       float **rho,
		       float **tau,
		       float **tauo,
		       float *bndr,
		       float *wlt,
		       checkpoint *buffer,/*structure checkpoint, size=nob*/
		       int *checkpt,/*array of checkpoint, size=noc*/
		       int *ib_,
		       int *ic_,
		       int *sxz,
		       int is,
		       int it)
/*<remodeling from current buffer to step it>*/
{
  int ib, ic, itt;

  ib=*ib_;
  ic=*ic_;
  if(it>buffer[ib].no){
    //remodelling from buffer[ib].no to step it

    memcpy(p[0],buffer[ib].p[0], nzpad*nxpad*sizeof(float));
    memcpy(r[0],buffer[ib].r[0], nzpad*nxpad*sizeof(float));
    memcpy(vz[0],buffer[ib].vz[0], nzpad*nxpad*sizeof(float));
    memcpy(vx[0],buffer[ib].vx[0], nzpad*nxpad*sizeof(float));
    for(itt=buffer[ib].no+1; itt<=it; itt++){
      add_source(&sxz[is], p, 1, &wlt[itt], true);
      step_forward(p, r, vz, vx, vv, rho, tau, tauo);
      apply_sponge(p, bndr);
      apply_sponge(r, bndr);
      apply_sponge(vx, bndr);
      apply_sponge(vz, bndr);
    }
  }else if(it==buffer[ib].no){
    //read directly from buffer[ib].no
    memcpy(p[0],buffer[ib].p[0], nzpad*nxpad*sizeof(float));
    memcpy(r[0],buffer[ib].r[0], nzpad*nxpad*sizeof(float));
    memcpy(vz[0],buffer[ib].vz[0], nzpad*nxpad*sizeof(float));
    memcpy(vx[0],buffer[ib].vx[0], nzpad*nxpad*sizeof(float));

    //while checking checkpt[ic-1] and buffer[ib-1].no
    if(ic>=1){//buffer[ib]/checkpt[ic] is not the first buffer/checkpoint
      if(checkpt[ic-1]>buffer[ib-1].no){
	//modeling from buffer[ib-1].no to checkpt[ic-1], store it to buffer[ib]
	memcpy(buffer[ib].p[0],buffer[ib-1].p[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].r[0],buffer[ib-1].r[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].vz[0],buffer[ib-1].vz[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].vx[0],buffer[ib-1].vx[0], nzpad*nxpad*sizeof(float));
	for(itt=buffer[ib-1].no+1; itt<=checkpt[ic-1]; itt++){
	  add_source(&sxz[is], buffer[ib].p, 1, &wlt[itt], true);
          step_forward(buffer[ib].p, buffer[ib].r, buffer[ib].vz, buffer[ib].vx, vv, rho, tau, tauo);
          apply_sponge(buffer[ib].p, bndr);
          apply_sponge(buffer[ib].r, bndr);
          apply_sponge(buffer[ib].vx, bndr);
          apply_sponge(buffer[ib].vz, bndr);
	}	
	//update buffer[ib] with checkpt[ic-1], deprecate checkpt[ic]
	buffer[ib].no=checkpt[ic-1];
	ic=ic-1;
      }else if(checkpt[ic-1]==buffer[ib-1].no){
	//deprecate last checkpt[ic] and buffer[ib] simultaneously
	ib=ib-1;
	ic=ic-1;
      }
    }
  }
  *ib_=ib;
  *ic_=ic;
}

void buffer_free(checkpoint *buffer, int nob)
/*< free buffer points  >*/
{
  int ib;
  
  for(ib=0; ib<nob; ib++)   {
    free(buffer[ib].p[0]); free(buffer[ib].p);
    free(buffer[ib].r[0]); free(buffer[ib].r);
    free(buffer[ib].vz[0]); free(buffer[ib].vz);
    free(buffer[ib].vx[0]); free(buffer[ib].vx);
  }
}


int main(int argc, char* argv[])
{
  bool verb, csdgather;
  int it, kt, is, i1, i2, *sxz, *gxz;
  int distz, distx, jsx, jsz, jgx, jgz, sxbeg, szbeg, gxbeg, gzbeg, tdmute;
  int noc, nob, ib, ic, *checkpt;
  checkpoint *buffer;
  float tmp, vmute, *wlt, *bndr;
  float **rho, **tau, **tauo, **v0, **vv, **sp, **sr, **svz, **svx, **gp, **gr, **gvz, **gvx, **dcal, **imag1, **imag2, **imag;
  sf_file Fv, Frho, Ftau, Ftauo, Fw, Fp1, Fp2;

  sf_init(argc,argv);
#ifdef _OPENMP
  omp_init();
#endif

  Fv=sf_input("in");/* veloctiy model */
  Frho=sf_input("rho");/* density */
  Ftau=sf_input("tau");/* tau, computed according to quality factor Q */
  Ftauo=sf_input("tauo");/* tauo, computed according to quality factor Q */
  Fw=sf_output("out");/* image */
  Fp1=sf_output("p1");/* forward wavefield at kt */
  Fp2=sf_output("p2");/* backward reconstructed wavefied by checkpoints at kt */

  if(!sf_getbool("verb",&verb)) verb=false;   
  /* verbosity */
  if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
  /* veloctiy model: nz */
  if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
  /* veloctiy model: nx */
  if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
  /* veloctiy model: dz */
  if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
  /* veloctiy model: dx */
  if (!sf_getint("nb",&nb)) nb=20; 
  /* thickness of PML ABC */
  if (!sf_getint("nt",&nt)) sf_error("nt required");
  /* number of time steps */
  if (!sf_getfloat("dt",&dt)) sf_error("dt required");
  /* time sampling interval */
  if (!sf_getfloat("fm",&fm)) fm=20.0; 
  /*dominant freq of Ricker wavelet */
  if (!sf_getint("ns",&ns))   sf_error("no ns");	
  /* number of shots */
  if (!sf_getint("ng",&ng))   sf_error("no ng");	
  /* number of geophones/receivers per shot */
  if (!sf_getint("kt",&kt)) sf_error("kt required");
  /* output px and pz component at kt */
  if (!sf_getint("nob",&nob)) nob=(int)log2f(nt);
  /*number of buffers, default=optimal value */
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

  noc=10*nob;//number of checkpoints 
  if (nob>=nt ||noc>=nt) sf_error("make sure: nob << noc << nt");
  _dx=1./dx;
  _dz=1./dz;
  nzpad=nz+2*nb;
  nxpad=nx+2*nb;

  buffer=(checkpoint*)sf_alloc(nob,sizeof(checkpoint));
  checkpt=sf_intalloc(noc);
  wlt=sf_floatalloc(nt);
  v0=sf_floatalloc2(nz,nx); 	
  rho=sf_floatalloc2(nzpad, nxpad);
  tau=sf_floatalloc2(nzpad, nxpad);
  tauo=sf_floatalloc2(nzpad, nxpad);
  vv=sf_floatalloc2(nzpad, nxpad);
  sp =sf_floatalloc2(nzpad, nxpad);
  sr =sf_floatalloc2(nzpad, nxpad);
  svz=sf_floatalloc2(nzpad, nxpad);
  svx=sf_floatalloc2(nzpad, nxpad);
  gp =sf_floatalloc2(nzpad, nxpad);
  gr =sf_floatalloc2(nzpad, nxpad);
  gvz=sf_floatalloc2(nzpad, nxpad);
  gvx=sf_floatalloc2(nzpad, nxpad);
  imag1=sf_floatalloc2(nz, nx);
  imag2=sf_floatalloc2(nz, nx);
  imag=sf_floatalloc2(nz, nx);
  bndr=sf_floatalloc(nb);
  sxz=sf_intalloc(ns);
  gxz=sf_intalloc(ng);
  dcal=sf_floatalloc2(ng,nt);

  buffer_init(checkpt, buffer, noc, nob, nt);
  for(it=0;it<nt;it++){
    tmp=SF_PI*fm*(it*dt-1.0/fm);
    tmp*=tmp;
    wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
  }
  sf_floatread(v0[0],nz*nx,Fv);
  expand2d(vv, v0);
  sf_floatread(v0[0],nz*nx,Frho);
  expand2d(rho, v0);
  sf_floatread(v0[0],nz*nx,Ftau);
  expand2d(tau, v0);
  sf_floatread(v0[0],nz*nx,Ftauo);
  expand2d(tauo, v0);
  memset(sp[0],0,nzpad*nxpad*sizeof(float));
  memset(sr[0],0,nzpad*nxpad*sizeof(float));
  memset(svx[0],0,nzpad*nxpad*sizeof(float));
  memset(svz[0],0,nzpad*nxpad*sizeof(float));
  memset(gp[0],0,nzpad*nxpad*sizeof(float));
  memset(gr[0],0,nzpad*nxpad*sizeof(float));
  memset(gvx[0],0,nzpad*nxpad*sizeof(float));
  memset(gvz[0],0,nzpad*nxpad*sizeof(float));
  memset(imag1[0],0,nz*nx*sizeof(float));
  memset(imag2[0],0,nz*nx*sizeof(float));
  memset(imag[0],0,nz*nx*sizeof(float));
  for(ib=0;ib<nb;ib++){
    tmp=0.015*(nb-ib);
    bndr[ib]=expf(-tmp*tmp);
  }
  if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))  { 
    sf_error("sources exceeds the computing zone!"); 
  }
  sg_init(sxz, szbeg, sxbeg, jsz, jsx, ns);
  distx=sxbeg-gxbeg;
  distz=szbeg-gzbeg;
  if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz)){ 
    sf_error("geophones exceeds the computing area!"); 
  }
  if (csdgather && !((sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  
    && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz)) {
    sf_error("geophones exceeds the computing area!"); 
  }
  sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);

  for(is=0; is<ns; is++){
    memset(imag1[0],0,nz*nx*sizeof(float));
    memset(imag2[0],0,nz*nx*sizeof(float));

    //forward modeling while recording snapshots in buffer
    memset(sp[0],0,nzpad*nxpad*sizeof(float));
    memset(sr[0],0,nzpad*nxpad*sizeof(float));
    memset(svx[0],0,nzpad*nxpad*sizeof(float));
    memset(svz[0],0,nzpad*nxpad*sizeof(float));
    if (csdgather){
      gxbeg=sxbeg+is*jsx-distx;
      sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng);
    }
    for(ib=0,it=0; it<nt; it++){
      add_source(&sxz[is], sp, 1, &wlt[it], true);
      step_forward(sp, sr, svz, svx, vv, rho, tau, tauo);
      apply_sponge(sp, bndr);
      apply_sponge(sr, bndr);
      apply_sponge(svx, bndr);
      apply_sponge(svz, bndr);

      record_seis(dcal[it], gxz, sp, ng);
      muting(dcal[it], gzbeg,szbeg, gxbeg,sxbeg+is*jsx,jgx, it, tdmute, vmute);

      if(it==buffer[ib].no && ib<nob){//record snapshots in buffer
	memcpy(buffer[ib].p[0], sp[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].r[0], sr[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].vz[0], svz[0], nzpad*nxpad*sizeof(float));
	memcpy(buffer[ib].vx[0], svx[0], nzpad*nxpad*sizeof(float));
	ib++;
      }

      if(it==kt){//record in forward procedure
	window2d(v0,sp);    
	sf_floatwrite(v0[0],nz*nx,Fp1);
      }
    }

    //backward/reverse mode, apply imaging condition
    memset(sp[0],0,nzpad*nxpad*sizeof(float));
    memset(sr[0],0,nzpad*nxpad*sizeof(float));
    memset(svx[0],0,nzpad*nxpad*sizeof(float));
    memset(svz[0],0,nzpad*nxpad*sizeof(float));
    memset(gp[0],0,nzpad*nxpad*sizeof(float));
    memset(gr[0],0,nzpad*nxpad*sizeof(float));
    memset(gvx[0],0,nzpad*nxpad*sizeof(float));
    memset(gvz[0],0,nzpad*nxpad*sizeof(float));
    for(ib=nob-1,ic=noc-1,it=nt-1; it>-1; it--){
      buffer_remodeling(sp, sr, svz, svx, vv, rho, tau, tauo, bndr, wlt, buffer, checkpt, &ib, &ic, sxz, is, it);

	//the backpropagation operator should be the adjoint of forward modeling!
	//here we just use forward modeling operator for the time being
      add_source(gxz, gp, ng, dcal[it], true);
      step_forward(gp, gr, gvz, gvx, vv, rho, tau, tauo);
      apply_sponge(gp, bndr);
      apply_sponge(gr, bndr);
      apply_sponge(gvx, bndr);
      apply_sponge(gvz, bndr);

      if(it==kt){//record in backward procedure
	window2d(v0,sp);
	sf_floatwrite(v0[0],nz*nx,Fp2);
      }

      //apply correlation imaging condition
      for(i2=0; i2<nx; i2++){
	for(i1=0; i1<nz; i1++){
	  imag1[i2][i1]+=sp[i2+nb][i1+nb]*gp[i2+nb][i1+nb];
	  imag2[i2][i1]+=sp[i2+nb][i1+nb]*sp[i2+nb][i1+nb];
	}
      }
    }

    //normalized correlation imaging condition
    for(i2=0; i2<nx; i2++){
      for(i1=0; i1<nz; i1++){
	imag[i2][i1]+=imag1[i2][i1]/(imag2[i2][i1]+SF_EPS);
      }
    }
  }
  //output final image
  sf_floatwrite(imag[0],nz*nx,Fw);
  

  free(checkpt);
  buffer_free(buffer,nob);
  free(wlt);
  free(*rho); free(rho);
  free(*tau); free(tau);
  free(*tauo); free(tauo);
  free(*v0); free(v0);
  free(*vv); free(vv);
  free(*sp); free(sp);
  free(*sr); free(sr);
  free(*svx); free(svx);
  free(*svz); free(svz);
  free(*gp); free(gp);
  free(*gr); free(gr);
  free(*gvx); free(gvx);
  free(*gvz); free(gvz);
  free(bndr);
  free(sxz);
  free(gxz);
  free(*dcal); free(dcal);
  free(*imag1); free(imag1);
  free(*imag2); free(imag2);
  free(*imag); free(imag);

  exit(0);
}

