/*A 2D demo of elliptically-anisotropic wave propagation (4th order)
  Note: It is adapted according to Seregy Fomel's lecture on Seismic imaging.
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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

static int nb, nz, nx, nt, nzpad, nxpad;
static float dz, dx, dt, fm, c10, c20, c11, c12, c21, c22;
static float *bndr, *wlt;
static float **vz, **vx, **p0, **p1, **ptr=NULL;

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
/*< forward modeling step >*/
{
    int ix,iz;
    float tmp1, tmp2;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz,tmp1,tmp2)			\
    shared(p1,p0,vz,vx,nzpad,nxpad)
#endif	
    for (ix=2; ix < nxpad-2; ix++) 
	for (iz=2; iz < nzpad-2; iz++) 
	{
	    tmp1 =	c11*(p1[ix][iz-1]+p1[ix][iz+1])+
		c12*(p1[ix][iz-2]+p1[ix][iz+2])+
		c10*p1[ix][iz];
	    tmp2 =	c21*(p1[ix-1][iz]+p1[ix+1][iz])+
		c22*(p1[ix-2][iz]+p1[ix+2][iz])+
		c20*p1[ix][iz];
	    p0[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vz[ix][iz]*tmp1+vx[ix][iz]*tmp2;
	}
}

void apply_sponge(float**p0, float **p1)
/*< apply absorbing boundary condition >*/
{
    int ix,iz;

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(bndr,p0,p1)
#endif
    for(ix=0; ix<nxpad; ix++)
    {
	for(iz=0;iz<nb;iz++){	/* top ABC */
	    p0[ix][iz]=bndr[iz]*p0[ix][iz];
	    p1[ix][iz]=bndr[iz]*p1[ix][iz];
	}
	for(iz=nz+nb;iz<nzpad;iz++){/* bottom ABC*/
	    p0[ix][iz]=bndr[nzpad-iz-1]*p0[ix][iz];
	    p1[ix][iz]=bndr[nzpad-iz-1]*p1[ix][iz];
	}
    }

#ifdef _OPENMP
#pragma omp parallel for			\
    private(ix,iz)				\
    shared(bndr,p0,p1)
#endif
    for(iz=0; iz<nzpad; iz++)
    {
	for(ix=0;ix<nb;ix++){	/* left ABC */
	    p0[ix][iz]=bndr[ix]*p0[ix][iz];
	    p1[ix][iz]=bndr[ix]*p1[ix][iz];
	}	
	for(ix=nx+nb;ix<nxpad;ix++){/* right ABC */
	    p0[ix][iz]=bndr[nxpad-ix-1]*p0[ix][iz];
	    p1[ix][iz]=bndr[nxpad-ix-1]*p1[ix][iz];
	}	
    }
}


int main(int argc, char* argv[])
{
    bool verb;
    int jt, ft, it, ib, ix, iz, sx, sz;
    float tmp;
    float **vz0, **vx0;
    sf_file Fvz, Fvx, Fw;

    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    Fvz = sf_input("in");/* veloctiy-z */
    Fvx = sf_input("vx");/* veloctiy-x */
    Fw = sf_output("out");/* wavefield snaps */

    if(!sf_getbool("verb",&verb))  verb=false;   
    /* verbosity */
    if (!sf_histint(Fvz,"n1",&nz)) sf_error("No n1= in input");
    /* veloctiy model: nz */
    if (!sf_histint(Fvz,"n2",&nx)) sf_error("No n2= in input");
    /* veloctiy model: nx */
    if (!sf_histfloat(Fvz,"d1",&dz)) sf_error("No d1= in input");
    /* veloctiy model: dz */
    if (!sf_histfloat(Fvz,"d2",&dx)) sf_error("No d2= in input");
    /* veloctiy model: dx */
    if (!sf_getint("nb",&nb)) nb=30; 
    /* thickness of sponge ABC */
    if (!sf_getint("nt",&nt)) sf_error("nt required");
    /* number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("dt required");
    /* time sampling interval */
    if (!sf_getfloat("fm",&fm)) fm=20.0;
    /*dominant freq of Ricker wavelet */
    if (!sf_getint("ft",&ft)) ft=0;
    /* first recorded time */
    if (!sf_getint("jt",&jt)) jt=1;
    /* time interval */

    sf_putint(Fw,"n1",nz);
    sf_putint(Fw,"n2",nx);
    sf_putint(Fw,"n3",(nt-ft)/jt);
    sf_putfloat(Fw,"d3",jt*dt);
    sf_putfloat(Fw,"o3",ft*dt);

    nzpad=nz+2*nb;
    nxpad=nx+2*nb;
    sx=nxpad/2;
    sz=nzpad/2;

    /*< initialize 4-th order fd coefficients >*/
    tmp = 1.0/(dz*dz);
    c11 = 4.0*tmp/3.0;
    c12= -tmp/12.0;
    tmp = 1.0/(dx*dx);
    c21 = 4.0*tmp/3.0;
    c22= -tmp/12.0;
    c10=-2.0*(c11+c12);
    c20=-2.0*(c21+c22);

    wlt=sf_floatalloc(nt);
    bndr=sf_floatalloc(nb);
    vz0=sf_floatalloc2(nz,nx); 	
    vx0=sf_floatalloc2(nz,nx); 
    vz=sf_floatalloc2(nzpad, nxpad);
    vx=sf_floatalloc2(nzpad, nxpad);
    p0=sf_floatalloc2(nzpad, nxpad);
    p1=sf_floatalloc2(nzpad, nxpad);

    for(it=0;it<nt;it++){
	tmp=SF_PI*fm*(it*dt-1.0/fm);tmp*=tmp;
	wlt[it]=(1.0-2.0*tmp)*expf(-tmp);
    }
    for(ib=0;ib<nb;ib++){
	tmp=0.015*(nb-ib);
	bndr[ib]=expf(-tmp*tmp);
    }
    sf_floatread(vz0[0],nz*nx,Fvz);
    sf_floatread(vx0[0],nz*nx,Fvx);
    expand2d(vz, vz0);
    expand2d(vx, vx0);
    for(ix=0;ix<nxpad;ix++){
	for(iz=0;iz<nzpad;iz++){
	    tmp=vz[ix][iz]*dt; vz[ix][iz]=tmp*tmp;
	    tmp=vx[ix][iz]*dt; vx[ix][iz]=tmp*tmp;
	}
    }
    memset(p0[0],0,nzpad*nxpad*sizeof(float));
    memset(p1[0],0,nzpad*nxpad*sizeof(float));

    for(it=0; it<nt; it++)
    {
	if(it>=ft)
	{
	    window2d(vx0,p0);
	    sf_floatwrite(vx0[0],nz*nx,Fw);
	}

	p1[sx][sz]+=wlt[it];
	step_forward(p0, p1);
	apply_sponge(p0,p1);
	ptr=p0; p0=p1; p1=ptr;

	if (verb) sf_warning("it=%d;",it);
    }

    free(wlt);
    free(*vz0); free(vz0);
    free(*vx0); free(vx0);
    free(*vz); free(vz);
    free(*vx); free(vx);
    free(*p0); free(p0);
    free(*p1); free(p1);
    free(bndr);

    exit(0);
}

