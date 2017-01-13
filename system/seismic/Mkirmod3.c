/* Kirchhoff 3-D modeling with analytical Green's functions. */
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
#include <float.h>
#include <math.h>

#include <rsf.h>

#include "kirmod.h"
#include "kirmod3.h"
#include "ricker.h"

int main(int argc, char* argv[]) 
{
    int nx,ny, nt, nsx,nsy, nhx,nhy, nc, nxyc, isx,isy, ihx,ihy, ix,iy, ic;
    int ns,nh, two, is, ih;
    bool absoff, verb;
    float ***rfl, ***rgd, ***crv, ***dipx, ***dipy, *trace;
    float slow, dx, x0, dy, y0, dt, t0, aper, x, y, dx1, dy1, dx2, dy2;
    float dsx, dsy, s0x, s0y, dhx, h0x, dhy, h0y, r0;
    float theta, ava, amp, obl, ***time, ***ampl, ***delt, **geom, freq;
    const char *type;
    velocity3 vel;
    ktable ts, tg;
    sf_file refl, curv, modl, head;
    
    sf_init(argc,argv);
    curv = sf_input("in");
    modl = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&ny)) sf_error("No n2= in input");
    if (!sf_histfloat(curv,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(curv,"o2",&y0)) sf_error("No o2= in input");
    if (!sf_histint  (curv,"n3",&nc)) nc=1; /* number of reflectors */
    nxyc = nx*ny*nc;

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */
    
    /*** Initialize trace ***/
    
    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time sampling */
    if (!sf_getfloat("t0",&t0)) t0=0.;
    /* time origin */

    if (!sf_getbool("absoff",&absoff)) absoff=false;
    /* y - h0x, h0y - are not in shot coordinate system */

    trace = sf_floatalloc(nt+1);

    sf_putint  (modl,"n1",nt);
    sf_putfloat(modl,"d1",dt);
    sf_putfloat(modl,"o1",t0);
    
    /*** Initialize shots ***/
    
    if (NULL != sf_getstring("head")) {
	/* source-receiver geometry (optional) */

	/* possibly irregular geometry in a file */
	head = sf_input("head");

	if (SF_FLOAT != sf_gettype(head)) sf_error("Need float type in head");
	if (!sf_histint(head,"n1",&two) || two != 2) 
	    sf_error("Need n1=2 in head");
	if (!sf_histint(head,"n2",&nh) || nh < 2) 
	    sf_error("Need n2 >=2 in head"); 
	nh--; /* nh includes shot */
	if (!sf_histint(head,"n3",&ns)) ns=1;

	sf_putint(modl,"n2",nh);
	sf_putint(modl,"n3",ns);	    
    } else {
	/* assume regular geometry */
	head = NULL;

	if (!sf_getint("nsx",&nsx)) nsx=nx;
	/* number of inline shots */
	if (!sf_getfloat("s0x",&s0x)) s0x=x0;
	/* first inline shot */
	if (!sf_getfloat("dsx",&dsx)) dsx=dx;
	/* inline shot increment */
    
	sf_putfloat(modl,"o4",s0x);
	sf_putfloat(modl,"d4",dsx);
	sf_putint(modl,"n4",nsx);
    
	if (!sf_getint("nsy",&nsy)) nsy=ny;
	/* number of crossline shots */
	if (!sf_getfloat("s0y",&s0y)) s0y=y0;
	/* first crossline shot */
	if (!sf_getfloat("dsy",&dsy)) dsy=dy;
	/* crossline shot increment */
	
	sf_putfloat(modl,"o5",s0y);
	sf_putfloat(modl,"d5",dsy);
	sf_putint(modl,"n5",nsy);
    
	ns = nsx*nsy;

	/*** Initialize offsets ***/
	
	if (!sf_getint  ("nhx",&nhx)) nhx=nx;
	/* number of inline offsets */
	if (!sf_getfloat("h0x",&h0x)) h0x=0.;
	/* first inline offset */
	if (!sf_getfloat("dhx",&dhx)) dhx=dx;
	/* inline offset increment */

	sf_putint  (modl,"n2",nhx);
	sf_putfloat(modl,"o2",h0x);
	sf_putfloat(modl,"d2",dhx);
	
	if (!sf_getint  ("nhy",&nhy)) nhy=ny;
	/* number of crossline offsets */
	if (!sf_getfloat("h0y",&h0y)) h0y=0.;
	/* first crossline offset */
	if (!sf_getfloat("dhy",&dhy)) dhy=dy;
	/* crossline offset increment */

	sf_putint  (modl,"n3",nhy);
	sf_putfloat(modl,"o3",h0y);
	sf_putfloat(modl,"d3",dhy);

	nh = nhx*nhy;
    }

    /*** Initialize reflector ***/

    crv = sf_floatalloc3(nx,ny,nc);
    rfl = sf_floatalloc3(nx,ny,nc);
    rgd = sf_floatalloc3(nx,ny,nc);
    dipx = sf_floatalloc3(nx,ny,nc);
    dipy = sf_floatalloc3(nx,ny,nc);

    sf_floatread(crv[0][0],nxyc,curv);
    
    /* reflectivity (A) */
    if (NULL != sf_getstring("refl")) {
	refl = sf_input("refl");
	sf_floatread(rfl[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	if (!sf_getfloat("r0",&r0)) r0=1.;
	/* constant reflectivity */
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][iy][ix] = r0;
		}
	    }
	}
    }

    /* AVO gradient (B/A) */
    if (NULL != sf_getstring("rgrad")) {
	refl = sf_input("rgrad");
	sf_floatread(rgd[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    rgd[ic][iy][ix] = 0.;
		}
	    }
	}
    }

    /* reflector inline dip */
    if (NULL != sf_getstring("dipx")) {
	refl = sf_input("dipx");
	sf_floatread(dipx[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    dipx[ic][iy][ix] = 0.;
		}
	    }
	}
    }

    /* reflector crossline dip */
    if (NULL != sf_getstring("dipy")) {
	refl = sf_input("dipy");
	sf_floatread(dipy[0][0],nxyc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    dipy[ic][iy][ix] = 0.;
		}
	    }
	}
    }

    /*** Initialize velocity ***/
    vel  = (velocity3) sf_alloc(1,sizeof(*vel));

    if (!sf_getfloat("vel",&(vel->v0))) sf_error("Need vel=");
    /* velocity */
    
    if (!sf_getfloat("gradx",&(vel->gx))) (vel->gx)=0.;
    if (!sf_getfloat("grady",&(vel->gy))) (vel->gy)=0.;
    if (!sf_getfloat("gradz",&(vel->gz))) (vel->gz)=0.;
    /* velocity gradient */

    type = sf_getstring("type");
    /* type of velocity ('c': constant, 's': linear sloth, 'v': linear velocity) */
    if (NULL==type) {
	type= ((vel->gx)==0. && (vel->gy)==0. && (vel->gz)==0.)?"const":"veloc";
    } else if ((vel->gx)==0. && (vel->gy)==0. && (vel->gz)==0.) {
	type = "const"; 
    } else if ('s'==type[0]) {
	/* linear slowness squared */
	
	slow = 1./((vel->v0)*(vel->v0));
	/* slowness squared */
	(vel->gx) *= -2.*slow/(vel->v0);
	(vel->gy) *= -2.*slow/(vel->v0);
	(vel->gz) *= -2.*slow/(vel->v0);
	(vel->v0) = slow;     
    } else if ('v' != type[0]) {
	sf_error("Unknown type=%s",type);
    }
	
    if (!sf_getfloat("refx",&(vel->x0))) (vel->x0)=x0;
    if (!sf_getfloat("refy",&(vel->y0))) (vel->y0)=y0;
    if (!sf_getfloat("refz",&(vel->z0))) (vel->z0)=0.;
    /* reference coordinates for velocity */

    if (!sf_getfloat("aper",&aper)) aper=hypotf(nx*dx,ny*dy);
    /* aperture */
    aper *= aper;
    
    /*** Allocate space ***/    
    kirmod3_init(x0, dx, y0, dy, vel, type[0], crv, dipx, dipy);
    ts = (ktable) sf_alloc(1,sizeof(*ts));
    tg = (ktable) sf_alloc(1,sizeof(*tg));


    /*** Initialize stretch ***/
    sf_aastretch_init (false, nt, t0, dt, nxyc);

    time = sf_floatalloc3(nx,ny,nc);
    ampl = sf_floatalloc3(nx,ny,nc);
    delt = sf_floatalloc3(nx,ny,nc);
    geom = sf_floatalloc2(2,nh+1);

    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    /*** Main loop ***/
    /* loop over sources */
    for (is=0; is < ns; is++) { 
	if (NULL == head) { /* regular */
	    isy = is/nsx;
	    isx = is - isy*nsx;

	    geom[nh][0] = s0x + isx*dsx;
	    geom[nh][1] = s0y + isy*dsy;
	} else { /* irregular */
	    sf_floatread(geom[0],2*(nh+1),head);
	}


	/* loop over offsets */
	for (ih=0; ih < nh; ih++) { 
            if (verb && 0==ih%(nh/10+1)) sf_warning("source %d of %d, receiver %d of %d;",is+1,ns,ih+1,nh);
	    if (NULL == head) { /* regular */		
		ihy = ih/nhx;
		ihx = ih - ihy*nhx;

		if (absoff) {
                    geom[ih][0] = h0x + ihx*dhx;
                    geom[ih][1] = h0y + ihy*dhy;
		} else {
		    geom[ih][0] = geom[nh][0] + h0x + ihx*dhx;
  		    geom[ih][1] = geom[nh][1] + h0y + ihy*dhy;
		}
	    }	

	    /* loop over surface */
	    for (iy=0; iy < ny; iy++) { for (ix=0; ix < nx; ix++) {
		x = x0+ix*dx; dx1 = x-geom[nh][0]; dx2 = x-geom[ih][0];
		y = y0+iy*dy; dy1 = y-geom[nh][1]; dy2 = y-geom[ih][1];

		/* skip if outside the aperture */
		if (aper < dx1*dx1 + dy1*dy1 ||
		    aper < dx2*dx2 + dy2*dy2) {
		    for (ic=0; ic < nc; ic++) {
			time[ic][iy][ix] = t0+2.*nt*dt;
			ampl[ic][iy][ix] = 0.;
			delt[ic][iy][ix] = dt;
		    }
		    continue;
		}

		/* loop over surface number */
		for (ic=0; ic < nc; ic++) {
		    kirmod3_map(ts,geom[nh],ix,iy,ic);
		    kirmod3_map(tg,geom[ih],ix,iy,ic);
		    
		    time[ic][iy][ix] = ts->t + tg->t;
		    
		    tg->an /= sqrtf(1.-(tg->tn)*(tg->tn));
		    ts->an /= sqrtf(1.-(ts->tn)*(ts->tn));
                    if (rgd[ic][iy][ix] != 0.) { 
			theta = hypotf((tg->an)*(tg->tx)-(ts->an)*(ts->tx),
					(tg->an)*(tg->ty)-(ts->an)*(ts->ty));

			/* AVA */
			theta = sinf(0.5*theta);
			ava = rfl[ic][iy][ix]+rgd[ic][iy][ix]*theta*theta;
		    } else
			ava = rfl[ic][iy][ix];

		    /* obliguity */
		    obl = 0.5*(ts->tn + tg->tn);
		    
		    /* Geometrical spreading */
		    amp = ts->a * tg->a + FLT_EPSILON;
		    
		    ampl[ic][iy][ix] = ava*obl*dx/amp;
		    delt[ic][iy][ix] = SF_MAX(fabsf(ts->tx+tg->tx)*dx,
					      fabsf(ts->ty+tg->ty)*dy); 
		}
	    }}
	    sf_aastretch_define (time[0][0],delt[0][0],NULL);
	    sf_aastretch_lop (false,false,nxyc,nt,ampl[0][0],trace);
	    
	    /* convolve with Ricker wavelet */
	    sf_freqfilt(nt,trace);
	    
	    sf_floatwrite(trace,nt,modl);
	} /* ih */

    } /* is */
    if (verb) sf_warning(".");


    exit(0);
}
