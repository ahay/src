/* Normal moveout in tau-p domain.

Takes: < input.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

#include "fint1.h"

int main (int argc, char* argv[])
{
    fint1 nmo;
    int it, iz, ip, ix, nt, np, nx, nw;
    float dt, t0, dp, p0, p, f, ft;
    float *trace, *vel;
    sf_file taup, nmod, velocity;

    sf_init (argc,argv);
    taup = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(taup)) sf_error("Need float input");
    if (!sf_histint(taup,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(taup,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(taup,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histint(taup,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(taup,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(taup,"o2",&p0)) sf_error("No o2= in input");   

    /* half velocity */
    p0 *= 0.5;
    dp *= 0.5;

    nx = sf_leftsize(taup,2);

    if (!sf_getint("extend",&nw)) nw=4;
    /* interpolation accuracy */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);

    nmo = fint1_init (nw, nt);

    for (ix=0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	

	for (ip=0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p *= p;

	    sf_floatread (trace,nt,taup);
	    fint1_set(nmo,trace);

	    f = 0.;
	    for (it=0; it < nt; it++) {
		ft = vel[it];
		ft = 1.-p*ft*ft;
		if (ft < 0.) {
		    for (; it < nt; it++) {
			trace[it]=0.;			
		    }
		    break;
		} 

		f += sqrtf(ft);

		ft = (f-t0)/dt;
		iz = f;

		if (iz >= 0 && iz < nt) {
		    trace[it] = fint1_apply(nmo,iz,f-iz,false);
		} else {
		    trace[it] = 0.;
		}
	    }
	    sf_floatwrite (trace,nt,nmod);
	}
    }
    
    sf_close();
    exit (0);
}

/* 	$Id: Mtaupmo.c,v 1.5 2004/04/19 21:51:46 fomels Exp $	 */
