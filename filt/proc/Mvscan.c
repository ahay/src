/* Velocity analysis.

Takes: < cmp.rsf > scan.rsf
*/

#include <math.h>

#include <rsf.h>

#include "fint1.h"

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool sembl;
    int it,ih,ix,iv, nt,nh,nx,nv, ib,ie,nb,i, nw, iz, CDPtype;
    float dt, dh, t0, h0, v0, dv, h, v, num, den, t, dy;
    float *trace, **stack, **stack2;
    sf_file cmp, scan;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");

    CDPtype=1;
    if (sf_histfloat(cmp,"d3",&dy)) {
	CDPtype=0.5+dh/dy;
	if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    } 	    
    sf_warning("CDPtype=%d",CDPtype);

    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    if (!sf_getint("nv",&nv)) sf_error("Need nv=");

    sf_putfloat(scan,"o2",v0+dv);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);
    sf_putstring(scan,"label2","velocity");

    trace = sf_floatalloc(nt);
    stack =  sf_floatalloc2(nt,nv);
    stack2 = sembl? sf_floatalloc2(nt,nv) : NULL;

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt);

    for (ix=0; ix < nx; ix++) {
	sf_warning("cmp %d of %d",ix+1,nx);

	for (it=0; it < nt*nv; it++) {
	    stack[0][it] = 0.;
	    if (sembl) stack2[0][it] = 0.;
	}

	for (ih=0; ih < nh; ih++) {
	    h = h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    h = 4. * h * h; 
	    sf_read(trace,sizeof(float),nt,cmp); 

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*nh;
	    }
	    fint1_set(nmo,trace);

	    for (iv=0; iv < nv; iv++) {
		v = v0 + iv * dv;
		v = 1./(v*v);

		for (it=0; it < nt; it++) {
		    t = t0+it*dt;
		    t = sqrtf(t*t + h*v);
		    t = (t-t0)/dt;
		    iz = t;

		    if (iz >=0 && iz < nt) {
			t = fint1_apply(nmo,iz,t-iz,false);
			stack[iv][it] += t;
			if (sembl) stack2[iv][it] += t*t;
		    }
		}
	    } /* v */
	} /* h */
	
	if (sembl) {
	    for (iv=0; iv < nv; iv++) {
		for (it=0; it < nt; it++) {
		    ib = it-nb;
		    ie = it+nb+1;
		    if (ib < 0) ib=0;
		    if (ie > nt) ie=nt;
		    num = 0.;
		    den = 0.;
		    for (i=ib; i < ie; i++) {
			num += stack[iv][i]*stack[iv][i];
			den += stack2[iv][i];
		    }
		    trace[it] = (den > 0.)? num/den: 0.;
		}
		sf_write(trace,sizeof(float),nt,scan);
	    }
	} else {
	    sf_write (stack[0],sizeof(float),nt*nv,scan);
	}
    } /* x */

    sf_close();
    exit(0);
}

/* 	$Id: Mvscan.c,v 1.2 2004/03/22 05:43:25 fomels Exp $	 */


