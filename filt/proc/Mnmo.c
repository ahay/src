/* Normal moveout.

Takes: < gather.rsf velocity=velocity.rsf offset=offset.rsf > nmod.rsf
*/

#include <math.h>

#include <rsf.h>

#include "fint1.h"

int main (int argc, char* argv[])
{
    fint1 nmo;
    int it, ix, iz, nt,nx,nw;
    float dt, t0, x, h, h0, f;
    float *trace, *vel, *off;
    sf_file cmp, nmod, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    offset = sf_input("offset");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    nx = sf_leftsize(cmp,1);

    if (!sf_getfloat ("h0",&h0)) h0=0.;
    /* reference offset */
    if (!sf_getint("extend",&nw)) nw=4;
    /* interpolation accuracy */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);
    off = sf_floatalloc(nx);

    nmo = fint1_init (nw, nt);

    sf_read (off,sizeof(float),nx,offset);
    for (ix = 0; ix < nx; ix++) {
	sf_read (trace,sizeof(float),nt,cmp);
	fint1_set(nmo,trace);

	sf_read (vel,sizeof(float),nt,velocity);	
	x = off[ix]; 
	h = x*x - h0*h0;

	for (it=0; it < nt; it++) {
	    f = t0 + it*dt;
	    f = f*f + h*vel[it]*vel[it];
	    if (f < 0.) {
		trace[it]=0.;
	    } else {
		f = sqrtf(f);
		iz = f;
		trace[it] = fint1_apply(nmo,iz,f-iz,false);
	    }
	}
	
	sf_write (trace,sizeof(float),nt,nmod);
    }

    exit (0);
}

/* 	$Id: Mnmo.c,v 1.3 2003/10/01 14:38:31 fomels Exp $	 */

