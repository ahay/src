/* 2D coustic time-domain FD modeling with different boundary conditions using OpenMP
*/
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

void oneway_abc(float **uo, float **um, float **vv, int nx, int nz, int nb, float dx, float dz, float dt, bool free);
void sponge_abc(float**   uu, int nx, int nz, int nb);

int main(int argc, char* argv[])
{
    /* Laplacian coefficients */
    float c0=-30./12.,c1=+16./12.,c2=- 1./12.;
 
    bool verb,free, ifoneway, ifsponge;           /* verbose flag */
    sf_file Fw=NULL,Fv=NULL,Fr=NULL,Fo=NULL; /* I/O files */
    sf_axis at,az,ax;    /* cube axes */
    int it,iz,ix,nb;        /* index variables */
    int nt,nz,nx;
    float dt,dz,dx,idx,idz,dt2;
 
    float  *ww,**vv,**rr;     /* I/O arrays*/
    float **um,**uo,**up,**ud;/* tmp arrays */

	/* Testing for OpenMP */
	double start_time, end_time;
 
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=0;
    if(! sf_getbool("free",&free)) free=false;
    if(! sf_getbool("ifoneway",&ifoneway)) ifoneway=true;
    if(! sf_getbool("ifsponge",&ifsponge)) ifsponge=true;
    if(! sf_getint("nb",&nb)) nb=5;
 
    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fv = sf_input ("vel");
    Fr = sf_input ("ref");
 
    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);
 
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);
 
    dt2 =    dt*dt;
    idz = 1/(dz*dz);
    idx = 1/(dx*dx);
 
    /* read wavelet, velocity & reflectivity */
    ww=sf_floatalloc(nt);     sf_floatread(ww   ,nt   ,Fw);
    vv=sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fv);
    rr=sf_floatalloc2(nz,nx); sf_floatread(rr[0],nz*nx,Fr);
 
    /* allocate temporary arrays */
    um=sf_floatalloc2(nz,nx);
    uo=sf_floatalloc2(nz,nx);
    up=sf_floatalloc2(nz,nx);
    ud=sf_floatalloc2(nz,nx);
 
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	}
    }

    /* MAIN LOOP */
    if(verb) fprintf(stderr,"\n");

	/* Starting timer */
    start_time = omp_get_wtime();

    for (it=0; it<nt; it++) 
{
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
 
	/* 4th order laplacian */
	if(ifoneway)
	{
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,iz)						\
    shared(ud,nb,uo,c0,c1,c2,nx,nz,idx,idz)
#endif
	for (iz=nb; iz<nz-nb; iz++) {
	    for (ix=nb; ix<nx-nb; ix++) {
		ud[ix][iz] = 
		    c0* uo[ix  ][iz  ] * (idx+idz) + 
		    c1*(uo[ix-1][iz  ] + uo[ix+1][iz  ])*idx +
		    c2*(uo[ix-2][iz  ] + uo[ix+2][iz  ])*idx +
		    c1*(uo[ix  ][iz-1] + uo[ix  ][iz+1])*idz +
		    c2*(uo[ix  ][iz-2] + uo[ix  ][iz+2])*idz;	  
	    }
	}
 	}else
	{
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,iz)						\
    shared(ud,uo,c0,c1,c2,nx,nz,idx,idz)
#endif
	for (iz=2; iz<nz-2; iz++) {
	    for (ix=2; ix<nx-2; ix++) {
		ud[ix][iz] = 
		    c0* uo[ix  ][iz  ] * (idx+idz) + 
		    c1*(uo[ix-1][iz  ] + uo[ix+1][iz  ])*idx +
		    c2*(uo[ix-2][iz  ] + uo[ix+2][iz  ])*idx +
		    c1*(uo[ix  ][iz-1] + uo[ix  ][iz+1])*idz +
		    c2*(uo[ix  ][iz-2] + uo[ix  ][iz+2])*idz;	  
	    }
	}	
	}

	/* inject wavelet */
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,iz)						\
    shared(nx,nz,ww,rr,ud,it)
#endif
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] -= ww[it] * rr[ix][iz];
	    }
	}
 
	/* scale by velocity */
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,iz)						\
    shared(ud,vv,nx,nz)
#endif
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] *= vv[ix][iz]*vv[ix][iz];
	    }
	}
 
	/* time step */
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,iz)						\
    shared(up,uo,um,ud,nx,nz,dt2)
#endif
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		up[ix][iz] = 
		    2*uo[ix][iz] 
		    - um[ix][iz] 
		    + ud[ix][iz] * dt2; 
 
		um[ix][iz] = uo[ix][iz];
		uo[ix][iz] = up[ix][iz];
	    }
	}
 
	/* one-way abc apply */
	if(ifoneway)
	{oneway_abc(uo,um,vv,nx,nz,nb,dx,dz,dt,free);}
	if(ifsponge)
	{sponge_abc(um,nx,nz,nb);
	 sponge_abc(uo,nx,nz,nb);}
 
	/* write wavefield to output */
	sf_floatwrite(uo[0],nz*nx,Fo);
    }

	/* Ending timer */
    end_time = omp_get_wtime();
	sf_warning("Elapsed time is %f.",end_time-start_time);

    if(verb) fprintf(stderr,"\n");    
    sf_close();
    exit (0);
}


void oneway_abc(float **uo, float **um, float **vv, int nx, int nz, int nb, float dx, float dz, float dt, bool free)
/* oneway - absorbing condition */
{
    int iz,ix,ib;
    float q;

#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,ib,iz,q)						\
    shared(nx,nz,nb,free,uo,um,dt,dz,vv)
#endif
    for(ix=0;ix<nx;ix++) {
	for(ib=0;ib<nb;ib++) {

	    /* top BC */
	    if(!free) { /* not free surface, apply ABC */
		iz = nb-ib;
		q = vv[ix][           nb  ] *dt/dz; 
		uo      [ix][iz  ] 
		    = um[ix][iz+1] 
		    +(um[ix][iz  ]
		    - uo[ix][iz+1]) *(1-q)/(1+q);
	    }

	    /* bottom BC */
	    iz = nz-nb+ib-1;
	    q = vv[ix][nz-nb-1] *dt/dz; 
	    uo      [ix][iz  ] 
		= um[ix][iz-1]
		+(um[ix][iz  ]
		- uo[ix][iz-1]) *(1-q)/(1+q);
	}
    }

#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix,ib,iz,q)						\
    shared(nx,nz,nb,uo,um,dt,dx,vv)
#endif
    for(iz=0;iz<nz;iz++) {
	for(ib=0;ib<nb;ib++) {

	    /* left BC */
	    ix = nb-ib;
	    q = vv[           nb  ][iz] *dt/dx;
	    uo      [ix  ][iz] 
		= um[ix+1][iz] 
		+(um[ix  ][iz]
		- uo[ix+1][iz]) *(1-q)/(1+q);

	    /* right BC */
	    ix = nx-nb+ib-1;
	    q = vv[nx-nb-1][iz] *dt/dx; 
	    uo      [ix  ][iz] 
		= um[ix-1][iz]
		+(um[ix  ][iz]
		- uo[ix-1][iz]) *(1-q)/(1+q);
	}
    }

}

void sponge_abc(float**   uu, int nx, int nz, int nb)
/* sponge - absorbing condition 
 Sponge boundary conditions multiply incoming wavefields
by smaller coefficients to attenuate the wavefield over time and space.

The sponge coefficients need to deviate from 1 very gradually to ensure
that there are no induced reflections caused by large impedance 
contrasts */
{
    int iz,ix,ib,ibz,ibx,sb;
    float *w,fb;
    w = sf_floatalloc(nb);
    sb = 4.0*nb;

#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ib,fb)						\
    shared(nb,sb,w)
#endif          
    for(ib=0; ib<nb; ib++) {
	fb = ib/(sqrt(2.0)*sb);
	w[ib] = exp(-fb*fb);
    }

    for(ib=0; ib<nb; ib++) {
	ibz = nz-ib-1;
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(ix)						\
    shared(nb,nx,ibz,ib,w,uu)
#endif   
	for(ix=0; ix<nx; ix++) {
	    uu[ix][ib ] *= w[nb-ib-1]; /*    top sponge */
	    uu[ix][ibz] *= w[nb-ib-1]; /* bottom sponge */
	}

	ibx = nx-ib-1;
#ifdef _OPENMP
#pragma omp parallel for default(none)						\
    private(iz)						\
    shared(nb,nz,ibx,ib,w,uu)
#endif   
	for(iz=0; iz<nz; iz++) {
	    uu[ib ][iz] *= w[nb-ib-1]; /*   left sponge */
	    uu[ibx][iz] *= w[nb-ib-1]; /*  right sponge */
	}

    }
}

