/* Normal moveout.

Takes: < gather.rsf velocity=velocity.rsf [offset=offset.rsf] > nmod.rsf
*/

#include <math.h>

#include <rsf.h>

#include "fint1.h"

int main (int argc, char* argv[])
{
    fint1 nmo;
    bool half;
    int it,ix,iz,ih, nt,nx,nw, nh, CDPtype;
    float dt, t0, h, h0, f, dh, dy;
    float *trace, *vel, *off;
    sf_file cmp, nmod, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");

    off = sf_floatalloc(nh);

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	sf_read (off,sizeof(float),nh,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");

	if (!sf_getbool("half",&half)) half=true;
	/* if y, the second axis is half-offset instead of full offset */
	
	if (half) {
	    dh *= 2.;
	    h0 *= 2.;
	}

	if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=0.5+0.5*dh/dy;
	    if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
	} 	    
	sf_warning("CDPtype=%d",CDPtype);

	for (ih = 0; ih < nh; ih++) {
	    off[ih] = h0 + ih*dh; 
	}	
    }

    nx = sf_leftsize(cmp,2);

    if (!sf_getfloat ("h0",&h0)) h0=0.;
    /* reference offset */
    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);

    nmo = fint1_init (nw, nt);
    
    for (ix = 0; ix < nx; ix++) {
	sf_read (vel,sizeof(float),nt,velocity);	

	for (ih = 0; ih < nh; ih++) {
	    sf_read (trace,sizeof(float),nt,cmp);
	    fint1_set(nmo,trace);
	    
	    h = off[ih] + (dh/CDPtype)*(ix%CDPtype); 
	    h = h*h - h0*h0;
	    
	    for (it=0; it < nt; it++) {
		f = t0 + it*dt;
		f = f*f + h/(vel[it]*vel[it]);
		if (f < 0.) {
		    trace[it]=0.;
		} else {
		    f = (sqrtf(f) - t0)/dt;
		    iz = f;
		    if (iz >= 0 && iz < nt) {
			trace[it] = fint1_apply(nmo,iz,f-iz,false);
		    } else {
			trace[it] = 0.;
		    }
		}
	    }
	    
	    sf_write (trace,sizeof(float),nt,nmod);
	}
    }
	
    sf_close();
    exit (0);
}

/* 	$Id: Mnmo.c,v 1.7 2004/04/02 02:23:02 fomels Exp $	 */

