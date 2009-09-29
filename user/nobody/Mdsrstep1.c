/* 2-D prestack modeling/migration with split-step DSR. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "dsr1.h"

int main (int argc, char *argv[])
{
    int nt;		/* number of time samples */
    int nz;		/* number of migrated time samples */
    int nx;		/* number of midpoints 	*/
    int nh;             /* number of offsets */
    int ix,it,iz;         /* loop counters 	*/
    int ntfft;	        /* fft size		*/
    int nw;		/* number of frequencies */	

    float dt;		/* time sampling interval 	*/
    float dz;		/* migrated time sampling interval */
    float dw;	        /* frequency sampling interval */
    float x0, dx;	/* spatial origin, sampling interval	*/
    float h0, dh;       /* offset origin, sampling interval */
    float **vt, *v, v0;	/* velocities		*/
    float *p,**q;	/* input, output data		*/

    float complex ***cp;	   /* complex input		*/

    bool inv;             /* modeling or migration        */
    bool depth;           /* depth or time                */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant          */  
    kiss_fftr_cfg fft;
    sf_file in, out, vel;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool("depth",&depth)) depth = false;
    /* depth or time migration */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stability parameter */

    if (NULL == sf_getstring("velocity")) {
	/* velocity file */
	if (!sf_getfloat("vel",&v0)) sf_error ("Need vel=");
	/* constant velocity (if no velocity file) */
	vel = NULL;	
    } else {
	vel = sf_input("velocity");
    }

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error ("No d1= in input");
	if (!sf_getint("nt",&nt)) sf_error ("Need nt=");
	/* Length of time axis (for modeling) */ 
	if (!sf_getfloat("dt",&dt)) sf_error ("Need dt=");
	/* Time sampling (for modeling) */
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);

	if (!sf_histint(in,"n2",&nx)) nx=1;
	if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&x0)) x0=0.;

	if (!sf_histint(in,"nh",&nh)) sf_error("Need nh=");
	/* Number of offsets (for modeling) */
	if (!sf_histfloat(in,"dh",&dh)) sf_error ("Need dh=");
	/* Offset sampling (for modeling) */
	if (!sf_histfloat(in,"h0",&h0)) h0=0.;
        /* Offset origin (for modeling) */
	
	sf_putint(out,"n2",nh);
	sf_putfloat(out,"d2",dh);
	sf_putfloat(out,"o2",h0);

	sf_putint(out,"n3",nx);
	sf_putfloat(out,"d3",dx);
	sf_putfloat(out,"o3",x0);
    } else { /* migration */
	if (!sf_histint(in,"n1",&nt)) sf_error ("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error ("No d1= in input");
	if (NULL == vel) {
	    if (!sf_getint("nz",&nz)) {
		/* number of steps in depth 
		   (for constant-velocity depth migration) */
		if (depth) sf_error ("Need nz=");
		nz = nt;
	    }
	    if (!sf_getfloat("dz",&dz)) {
		/* sampling in depth 
		   (for constant-velocity depth migration) */
		if (depth) sf_error ("Need dz=");
		dz = dt*v0;
	    }
	} else {
	    if (!sf_histint(vel,"n1",&nz)) 
		sf_error ("No n1= in velocity");
	    if (!sf_histfloat(vel,"d1",&dz)) 
		sf_error ("No d1= in velocity");
	}
	sf_putint(out,"n1",nz);
	sf_putfloat(out,"d1",dz);
	
	if (!sf_histint(in,"n2",&nh)) nh=1;
	if (!sf_histfloat(in,"d2",&dh)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&h0)) h0=0.;
	
	if (!sf_histint(in,"n3",&nx)) nx=1;
	if (!sf_histfloat(in,"d3",&dx)) sf_error ("No d3= in input");
	if (!sf_histfloat(in,"o3",&x0)) x0=0.;

	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);

	sf_putint(out,"n3",1);
    }
    
    vt = sf_floatalloc2(nz,nx);
    v = sf_floatalloc(nz);

    if (NULL == vel) {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		vt[ix][iz] = v0;
	    }
	}
    } else {
	sf_floatread(vt[0],nx*nz,vel);
	sf_fileclose(vel);
    }
    for (iz=0; iz < nz; iz++) {
	v[iz] = 0.;
	for (ix=0; ix < nx; ix++) {
	    vt[ix][iz] = 1./(vt[ix][iz]*vt[ix][iz]);
	    v[iz] += vt[ix][iz];
	}
	v[iz] /= nx;
    }

    /* determine wavenumber sampling (for real to complex FFT) */
    ntfft = nt*2;
    nw = ntfft/2+1;
    dw = 2.0*SF_PI/(ntfft*dt);
	
    /* allocate space */
    p = sf_floatalloc(ntfft);
    q = sf_floatalloc2(nz,nx);
    cp = sf_complexalloc3(nw,nh,nx);

    fft = kiss_fftr_alloc(ntfft,inv? 1: 0,NULL,NULL);

    for (ix=0; ix<nx; ix++) {
	if (inv) {
	    sf_floatread(q[ix],nz,in);
	} else { 
	    for (ih=0; ih < nh; ih++) {
		sf_floatread(p,nt,in);

		/* pad with zeros and Fourier transform t to w */
		for (it=nt; it<ntfft; it++) {
		    p[it] = 0.0;
		}

		kiss_fftr(fft, p, (kiss_fft_cpx *) cp[ix][ih]);

		for (it=0; it<nw; it++)
		    cp[ix][ih][it] /= ntfft;
	    }
	}
    }

    dsr1 (verb, inv, eps,  
	  nw, dw, 
	  nz, dz, 
	  nx, dx,
	  nh, dh, h0,
	  vt, v,
	  cp, q);
    
    for (ix=0; ix<nx; ix++) {
	if (inv) {
	    for (ih=0; ih < nh; ih++) {		
		/* Fourier transform w to t (including FFT scaling) */
		kiss_fftri(fft,(const kiss_fft_cpx *) cp[ix][ih], p);
		for (it=0; it<nt; it++)
		    p[it] /= ntfft;
		sf_floatwrite (p,nt,out);
	    }
	} else {
	    sf_floatwrite (q[ix],nz,out);
	}
    }

    exit (0);
}

/* 	$Id: Msstep1.c 703 2004-07-12 18:16:36Z fomels $	 */
