/* Time power gain 

Takes: < in.rsf > out.rsf

Does not estimate tpow automatically.
*/
#include <math.h>

#include <rsf.h>

int main(int argc, char *argv[])
{
    int nt, nx, n2, i2, it, ix;
    float tpow, xpow, *trace, *tgain=NULL, *xgain=NULL, dt, t0, dx, x0, x;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) nx=1;
    n2 = sf_leftsize(in,2);

    if (!sf_getfloat("tpow",&tpow)) tpow=2.;
    /* power of time */
    if (!sf_getfloat("xpow",&xpow)) xpow=0.;
    /* power of space */

    trace = sf_floatalloc(nt);

    if (0. != tpow) {
	if (!sf_histfloat(in,"d1",&dt)) dt=1.;
	if (!sf_histfloat(in,"o1",&t0)) t0=0.;

	tgain = sf_floatalloc(nt);
	for (it=0; it < nt; it++) {
	    tgain[it] = powf(t0+(it+1)*dt,tpow);
	}
    }

    if (0. != xpow) {
	if (!sf_histfloat(in,"d2",&dx)) dx=1.;
	if (!sf_histfloat(in,"o2",&x0)) x0=0.;

	xgain = sf_floatalloc(nx);
	for (ix=0; ix < nx; ix++) {
	    xgain[ix] = powf(x0+(ix+1)*dx,xpow);
	}
    }


    for (i2=0; i2 < n2; i2++) {
	for (ix=0; ix < nx; ix++) {
	    sf_floatread(trace,nt,in);
	    
	    if (0. != tpow) {
		for (it=0; it < nt; it++) {
		    trace[it] *= tgain[it];
		}
	    } 

	    if (0. != xpow) {
		x = xgain[ix];
		for (it=0; it < nt; it++) {
		    trace[it] *= x;
		}
	    }

	    sf_floatwrite(trace,nt,out);
	}
    }

    sf_close();
    exit(0);
}
