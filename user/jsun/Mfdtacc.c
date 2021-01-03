/* 2-D Fourth-order Finite-difference wave extrapolation with timing option (no ABC)*/
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
#include <math.h>
#include "timer.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    bool verb,timer;                         /* verbose and timer flag */
    sf_axis at,az,ax;                        /* cube axes */
    int nt, nz, nx, it, iz, ix;              /* dimension and index variables */
    float dt, dz, dx, idz, idx, dt2;         /* intervals */
    float **next, **curr, **prev, **lapl;    /* tmp arrays */
    float *ww,**vv,**rr;                     /* I/O arrays*/
    double time=0.,t0=0.,t1=0.;
    sf_file Fw=NULL,Fv=NULL,Fr=NULL,Fo=NULL; /* I/O files */
    /* Laplacian coefficients */
    float c0=-30./12.,c1=+16./12.,c2=- 1./12.;

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=0;
    if(! sf_getbool("timer",&timer)) timer=0;

    /* setup I/O files */
    Fw = sf_input("in");   /* source wavlet*/
    Fv = sf_input("vel");  /* velocity */
    Fr = sf_input("ref");  /* source location */
    Fo = sf_output("out"); /* output wavefield */

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    //sf_oaxa(Fo,at,3);

    dt2 =    dt*dt;
    idz = 1/(dz*dz);
    idx = 1/(dx*dx);
 
    /* allocate memory for wavelet, velocity & source location */
    ww = sf_floatalloc(nt);     sf_floatread(ww,nt,Fw);
    vv = sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fv);
    rr = sf_floatalloc2(nz,nx); sf_floatread(rr[0],nz*nx,Fr);
    
    /* allocate temporary arrays */
    next=sf_floatalloc2(nz,nx);
    curr=sf_floatalloc2(nz,nx);
    prev=sf_floatalloc2(nz,nx);
    lapl=sf_floatalloc2(nz,nx);
 
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    next[ix][iz]=0.;
	    curr[ix][iz]=0.;
	    prev[ix][iz]=0.;
	    lapl[ix][iz]=0.;
	}
    }

    /* Main loop: propagation in time */
    if(verb) fprintf(stderr,"\n");
    if(timer) t0 = gtod_timer();

    for (it=0; it < nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

#ifdef _OPENMP
#pragma omp parallel private(iz,ix) shared(next,curr,prev,lapl)
#endif
	{
	    /* 4th order laplacian */
#ifdef _OPENMP
#pragma omp for schedule(static) 
#endif
	    for (iz=2; iz<nz-2; iz++) {
		for (ix=2; ix<nx-2; ix++) {
		    lapl[ix][iz] = 
			c0* curr[ix  ][iz  ] * (idx+idz) + 
			c1*(curr[ix-1][iz  ] + curr[ix+1][iz  ])*idx +
			c2*(curr[ix-2][iz  ] + curr[ix+2][iz  ])*idx +
			c1*(curr[ix  ][iz-1] + curr[ix  ][iz+1])*idz +
			c2*(curr[ix  ][iz-2] + curr[ix  ][iz+2])*idz;	  
		}
	    }

	    
	    /* time stepping */
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
	    for (iz=0; iz<nz; iz++) {
		for (ix=0; ix<nx; ix++) {
		    next[ix][iz] = 
			2*curr[ix][iz] 
			- prev[ix][iz] 
			+ (lapl[ix][iz]+ww[it]*rr[ix][iz])*vv[ix][iz]*vv[ix][iz]*dt2; 
		    prev[ix][iz] = curr[ix][iz]; 
		    curr[ix][iz] = next[ix][iz];
		}
	    }
	}
        //sf_floatwrite(curr[0],nz*nx,Fo);
    }
    if(timer) 
    {
	t1 = gtod_timer();
	time = t1-t0;
	sf_warning("Time = %lf\n",time);
    }
    if(verb) fprintf(stderr,"\n"); 
    /* write final wavefield to output */
    sf_floatwrite(curr[0],nz*nx,Fo);
    //sf_close();
    exit(0); 
}           
