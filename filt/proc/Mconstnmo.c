/* NMO with a constant velocity.

Takes: < gather.rsf > nmoed.rsf
*/
#include <math.h>

#include <rsf.h>
 
#include "fint1.h"

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool inv;
    int it,ih,ix, nt,nh,nx, nw, iz;
    float dt, dh, t0, h0, v0, h, t, *trace;
    sf_file in, nmod;

    sf_init(argc,argv);
    in = sf_input("in");
    nmod = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input"); 
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input"); 
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input"); 
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input"); 

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse NMO */

    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* NMO velocity */
    sf_putfloat(nmod,"v0",v0);
    v0 *= 0.5;
    h0 /= v0; 
    dh /= v0;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt);

    for (ix=0; ix < nx; ix++) {
	for (ih=0; ih < nh; ih++) {
	    h = h0 + ih*dh;
	    h *= h;

	    sf_read (trace,sizeof(float),nt,in);
	    fint1_set(nmo,trace);
	
	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		t *= t;
		if (inv) {
		    t -= h;
		} else {
		    t += h;
		}
		if (t < 0.) {
		    trace[it] = 0.;
		} else {
		    t = (sqrtf(t) - t0)/dt;
		    iz = t;
		    if (iz >= 0 && iz < nt) { 
			trace[it] = fint1_apply(nmo,iz,t-iz,false);
		    } else {
			trace[it] = 0.;
		    }
		}
	    }
	    
	    sf_write (trace,sizeof(float),nt,nmod);
	} /* h */
    } /* x */

    sf_close();
    exit(0);
}

/* 	$Id: Mconstnmo.c,v 1.2 2004/03/22 05:43:24 fomels Exp $	 */


