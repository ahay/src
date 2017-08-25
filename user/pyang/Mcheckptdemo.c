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

static int nb, nz, nx, nzpad, nxpad, nt;
static float dz, dx, _dz, _dx, dt, fm;


typedef struct {
    int no;/* number/index of the checkpoint in time coordinate */
    float **p;/* p is redundant since we have px and pz: p=px+pz */
    float **px;
    float **pz;
    float **vx;
    float **vz;
} checkpoint;


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

    for(ix=0; ix<nxpad; ix++){
	x=0;
	if (ix>=0 && ix<nb)   x=(ix-nb)*dx;
	else if(ix>=nxpad-nb && ix<nxpad) x=(ix-(nxpad-nb-1))*dx;
	d2x[ix] = d0*x*x;
    }
    for(iz=0; iz<nzpad; iz++){
	z=0;
	if (iz>=0 && iz<nb)   z=(iz-nb)*dz;
	else if(iz>=nzpad-nb && iz<nzpad) z=(iz-(nzpad-nb-1))*dz; 
	d1z[iz] = d0*z*z;  
    }
}



void step_forward(float **p, float **pz, float **px, float **vz, float **vx, float **vv, float *d1z, float *d2x)
{
    int i1, i2;
    float tmp, diff1, diff2;

#ifdef _OPENMP
#pragma omp parallel for default(none)				\
    private(i1,i2,diff1,diff2)					\
    shared(nzpad, nxpad, p, vz, vx, dt, _dz, _dx, d1z, d2x)
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
	    vz[i2][i1]=((1.-0.5*dt*d1z[i1])*vz[i2][i1]+dt*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
	    vx[i2][i1]=((1.-0.5*dt*d2x[i2])*vx[i2][i1]+dt*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
	}

#ifdef _OPENMP
#pragma omp parallel for default(none)					\
    private(i1,i2,diff1, diff2, tmp)					\
    shared(nzpad, nxpad, vv, p, pz, px, vz, vx, dt, _dz, _dx, d1z, d2x)
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
	    pz[i2][i1]=((1.-0.5*dt*d1z[i1])*pz[i2][i1]+dt*tmp*_dz*diff1)/(1.+0.5*dt*d1z[i1]);
	    px[i2][i1]=((1.-0.5*dt*d2x[i2])*px[i2][i1]+dt*tmp*_dx*diff2)/(1.+0.5*dt*d2x[i2]);
	    p[i2][i1]=px[i2][i1]+pz[i2][i1];
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

    /*set checkpoints in time axis */
    for(ic=0; ic<noc; ic++)
	checkpt[ic]=nt*ic/noc;
  
    /* set buffer points: allocate and initialize variables at each buffer point */
    for(ib=0; ib<nob; ib++) {
	buffer[ib].p=sf_floatalloc2(nzpad, nxpad);
	buffer[ib].pz=sf_floatalloc2(nzpad, nxpad);
	buffer[ib].px=sf_floatalloc2(nzpad, nxpad);
	buffer[ib].vz=sf_floatalloc2(nzpad, nxpad);
	buffer[ib].vx=sf_floatalloc2(nzpad, nxpad);

	ic=(int)(noc-1)*ib/(nob-1);/* index in checkpt[]   */
	/* ib==nob-1, ic==noc-1, the last buffer is the last checkpoint! */
	buffer[ib].no=checkpt[ic];/*set buffer to be one of checkpoints */
	memset(buffer[ib].p[0], 0, nzpad*nxpad*sizeof(float));
	memset(buffer[ib].pz[0], 0, nzpad*nxpad*sizeof(float));
	memset(buffer[ib].px[0], 0, nzpad*nxpad*sizeof(float));
	memset(buffer[ib].vz[0], 0, nzpad*nxpad*sizeof(float));
	memset(buffer[ib].vx[0], 0, nzpad*nxpad*sizeof(float));
    }
}

void buffer_remodeling(float **p,
		       float **pz,
		       float **px,
		       float **vz,
		       float **vx,
		       float **vv,
		       float *d1z,
		       float *d2x,
		       float *wlt,
		       checkpoint *buffer,/*structure checkpoint, size=nob*/
		       int *checkpt,/*array of checkpoint, size=noc*/
		       int *ib_,
		       int *ic_,
		       int sz,
		       int sx,
		       int it)
/*<remodeling from current buffer to step it>*/
{
    int ib, ic, itt;

    ib=*ib_;
    ic=*ic_;
    if(it>buffer[ib].no){
	/*remodelling from buffer[ib].no to step it */

	memcpy(p[0],buffer[ib].p[0], nzpad*nxpad*sizeof(float));
	memcpy(pz[0],buffer[ib].pz[0], nzpad*nxpad*sizeof(float));
	memcpy(px[0],buffer[ib].px[0], nzpad*nxpad*sizeof(float));
	memcpy(vz[0],buffer[ib].vz[0], nzpad*nxpad*sizeof(float));
	memcpy(vx[0],buffer[ib].vx[0], nzpad*nxpad*sizeof(float));
	for(itt=buffer[ib].no+1; itt<=it; itt++){
	    p[sx][sz]+=wlt[itt];
	    step_forward(p, pz, px, vz, vx, vv, d1z, d2x);
	}
    }else if(it==buffer[ib].no){
	/*read directly from buffer[ib].no */
	memcpy(p[0],buffer[ib].p[0], nzpad*nxpad*sizeof(float));
	memcpy(pz[0],buffer[ib].pz[0], nzpad*nxpad*sizeof(float));
	memcpy(px[0],buffer[ib].px[0], nzpad*nxpad*sizeof(float));
	memcpy(vz[0],buffer[ib].vz[0], nzpad*nxpad*sizeof(float));
	memcpy(vx[0],buffer[ib].vx[0], nzpad*nxpad*sizeof(float));

	/*while checking checkpt[ic-1] and buffer[ib-1].no */
	if(ic>=1){/*buffer[ib]/checkpt[ic] is not the first buffer/checkpoint */
	    if(checkpt[ic-1]>buffer[ib-1].no){
		/*modeling from buffer[ib-1].no to checkpt[ic-1], store it to buffer[ib] */
		memcpy(buffer[ib].p[0],buffer[ib-1].p[0], nzpad*nxpad*sizeof(float));
		memcpy(buffer[ib].pz[0],buffer[ib-1].pz[0], nzpad*nxpad*sizeof(float));
		memcpy(buffer[ib].px[0],buffer[ib-1].px[0], nzpad*nxpad*sizeof(float));
		memcpy(buffer[ib].vz[0],buffer[ib-1].vz[0], nzpad*nxpad*sizeof(float));
		memcpy(buffer[ib].vx[0],buffer[ib-1].vx[0], nzpad*nxpad*sizeof(float));
		for(itt=buffer[ib-1].no+1; itt<=checkpt[ic-1]; itt++){
		    buffer[ib].p[sx][sz]+=wlt[itt];
		    step_forward(buffer[ib].p, buffer[ib].pz, buffer[ib].px, buffer[ib].vz, buffer[ib].vx, vv, d1z, d2x);
		}	
		/* update buffer[ib] with checkpt[ic-1], deprecate checkpt[ic] */
		buffer[ib].no=checkpt[ic-1];
		ic=ic-1;
	    }else if(checkpt[ic-1]==buffer[ib-1].no){
		/* deprecate last checkpt[ic] and buffer[ib] simultaneously */
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
	free(buffer[ib].pz[0]); free(buffer[ib].pz);
	free(buffer[ib].px[0]); free(buffer[ib].px);
	free(buffer[ib].vz[0]); free(buffer[ib].vz);
	free(buffer[ib].vx[0]); free(buffer[ib].vx);
    }
}


int main(int argc, char* argv[])
{
    bool verb;
    int jt, ft, it, kt, i2, i1, sx, sz;
    int noc, nob, ib, ic, *checkpt;
    checkpoint *buffer;
    float tmp, vmax;
    float *wlt, *d2x, *d1z;
    float **v0, **vv, **p, **pz, **px, **vz, **vx;
    sf_file Fv, Fw, Fp1, Fp2;

    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    Fv = sf_input("in");/* veloctiy model */
    Fw = sf_output("out");/* wavefield snaps */
    Fp1 = sf_output("p1");/* forward wavefield at kt */
    Fp2 = sf_output("p2");/* backward reconstructed wavefied by checkpoints at kt */

    if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
    if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
    if (!sf_getint("nb",&nb)) nb=20; /* thickness of PML ABC */
    if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    if (!sf_getfloat("fm",&fm)) fm=20.0; /*dominant freq of Ricker wavelet */
    if (!sf_getint("ft",&ft)) ft=0; /* first recorded time */
    if (!sf_getint("jt",&jt)) jt=1;	/* time interval */
    if(!sf_getbool("verb",&verb)) verb=false;    /* verbosity, if y, output px and pz */
    if (!sf_getint("kt",&kt)) sf_error("kt required"); /* output px and pz component at kt */
    if (!sf_getint("nob",&nob)) nob=(int)log2f(nt);/*number of buffers, default=optimal value */

    sf_putint(Fw,"n1",nz);
    sf_putint(Fw,"n2",nx);
    sf_putint(Fw,"n3",(nt-ft)/jt);
    sf_putfloat(Fw,"d3",jt*dt);
    sf_putfloat(Fw,"o3",ft*dt);

    noc=10*nob;/* number of checkpoints */
    if (noc>=nt) sf_error("make sure: nob << noc << nt");
    _dx=1./dx;
    _dz=1./dz;
    nzpad=nz+2*nb;
    nxpad=nx+2*nb;
    sx=nxpad/2;
    sz=nzpad/2;

    buffer=(checkpoint*)sf_alloc(nob,sizeof(checkpoint));
    checkpt=sf_intalloc(noc);
    wlt=sf_floatalloc(nt);
    v0=sf_floatalloc2(nz,nx); 	
    vv=sf_floatalloc2(nzpad, nxpad);
    p =sf_floatalloc2(nzpad, nxpad);
    pz=sf_floatalloc2(nzpad, nxpad);
    px=sf_floatalloc2(nzpad, nxpad);
    vz=sf_floatalloc2(nzpad, nxpad);
    vx=sf_floatalloc2(nzpad, nxpad);
    d1z=sf_floatalloc(nzpad);
    d2x=sf_floatalloc(nxpad);

    buffer_init(checkpt, buffer, noc, nob, nt);
    for(it=0;it<nt;it++){
	tmp=SF_PI*fm*(it*dt-1.0/fm);
	tmp*=tmp;
	wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
    }
    sf_floatread(v0[0],nz*nx,Fv);
    expand2d(vv, v0);
    memset(p [0],0,nzpad*nxpad*sizeof(float));
    memset(px[0],0,nzpad*nxpad*sizeof(float));
    memset(pz[0],0,nzpad*nxpad*sizeof(float));
    memset(vx[0],0,nzpad*nxpad*sizeof(float));
    memset(vz[0],0,nzpad*nxpad*sizeof(float));
    vmax=v0[0][0];
    for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	    vmax=SF_MAX(v0[i2][i1],vmax);
    pmlcoeff_init(d1z, d2x, vmax);

    /* forward modeling while recording snapshots in buffer */
    for(ib=0,it=0; it<nt; it++){
	p[sx][sz]+=wlt[it];
	step_forward(p, pz, px, vz, vx, vv, d1z, d2x);

	if(it==buffer[ib].no && ib<nob){/* record snapshots in buffer */
	    memcpy(buffer[ib].p[0], p[0], nzpad*nxpad*sizeof(float));
	    memcpy(buffer[ib].pz[0], pz[0], nzpad*nxpad*sizeof(float));
	    memcpy(buffer[ib].px[0], px[0], nzpad*nxpad*sizeof(float));
	    memcpy(buffer[ib].vz[0], vz[0], nzpad*nxpad*sizeof(float));
	    memcpy(buffer[ib].vx[0], vx[0], nzpad*nxpad*sizeof(float));
	    ib++;
	}

	if(it==kt){
	    window2d(v0,p);    
	    sf_floatwrite(v0[0],nz*nx,Fp1);
	}

	if(it>=ft) {
	    window2d(v0,p);
	    sf_floatwrite(v0[0],nz*nx,Fw);
	}
    }

    memset(p [0],0,nzpad*nxpad*sizeof(float));
    memset(px[0],0,nzpad*nxpad*sizeof(float));
    memset(pz[0],0,nzpad*nxpad*sizeof(float));
    memset(vx[0],0,nzpad*nxpad*sizeof(float));
    memset(vz[0],0,nzpad*nxpad*sizeof(float));
    /* backward */
    for(ib=nob-1,ic=noc-1,it=nt-1; it>-1; it--){
	/* source wavefield re-modeling in the interval between two checkpoints */
	buffer_remodeling(p, pz, px, vz, vx, vv, d1z, d2x, wlt, buffer, checkpt, &ib, &ic, sz, sx, it);

	if(it==kt){
	    window2d(v0,p);
	    sf_floatwrite(v0[0],nz*nx,Fp2);
	}
    }

    free(checkpt);
    buffer_free(buffer,nob);
    free(wlt);
    free(*v0); free(v0);
    free(*vv); free(vv);
    free(*p); free(p);
    free(*px); free(px);
    free(*pz); free(pz);
    free(*vx); free(vx);
    free(*vz); free(vz);
    free(d1z);
    free(d2x);

    exit(0);
}

