#include <rsf.h>

int main(int argc, char* argv[])
{
    int it, nt, ix, nx, ia, na;
    float *trace, *ofunc;
    float a, a0, da, t, t0, dt, s;
    sf_file in, out;

    /* initialization */
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* get trace parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("Need n1=");
    if (!sf_histfloat(in,"d1",&dt)) dt=1.;
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    /* get number of traces */
    nx = sf_leftsize(in,1);

    if (!sf_getint("na",&na)) na=1;
    /* number of alpha values */
    if (!sf_getfloat("da",&da)) da=0.;
    /* increment in alpha */
    if (!sf_getfloat("a0",&a0)) a0=0.;
    /* first value of alpha */

    /* change output data dimensions */
    sf_putint(out,"n1",na);
    sf_putint(out,"n2",1);
    sf_putfloat(out,"d1",da);
    sf_putfloat(out,"o1",a0);

    trace = sf_floatalloc(nt);
    ofunc = sf_floatalloc(na);

    /* initialize */
    for (ia=0; ia < na; ia++) {
	ofunc[ia] = 0.;
    }

    /* loop over traces */
    for (ix=0; ix < nx; ix++) {

	/* read data */
	sf_floatread(trace,nt,in);

	/* loop over alpha */
	for (ia=0; ia < na; ia++) {
	    a = a0+ia*da;

	    /* loop over time samples */
	    for (it=0; it < nt; it++) {
		t = t0+it*dt;

		/* apply gain t^alpha */
		s = trace[it]*powf(t,a);
		
                /* !!! MODIFY THE NEXT LINE !!! */
		ofunc[ia] += s*s; 
	    }
	}
    }

    /* write output */
    sf_floatwrite(ofunc,na,out);

    exit(0);
}
