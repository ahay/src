#include <rsf.h>

#include "tridiagonal.h"

int main(int argc, char* argv[])
{
    tris slv;
    int nt, ns, nx, ix, is, it, *ipick;
    float *trace, *pick, *ampl, *diag, *offd, s0, ds, eps, lam, t, asum;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ns)) ns=1;
    nx = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o2",&s0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");
    sf_putint(out,"n2",1);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getfloat("lam",&lam)) lam=0.01;

    trace = sf_floatalloc(nt);
    pick = sf_floatalloc(nt);
    ipick = sf_intalloc(nt);
    offd = sf_floatalloc(nt);
    diag = sf_floatalloc(nt);
    ampl = sf_floatalloc(nt);

    slv = tridiagonal_init (nt);

    for (it=0; it < nt; it++) {
	offd[it] = -eps;
    }

    for (ix = 0; ix < nx; ix++) {
	for (it = 0; it < nt; it++) {
	    ampl[it] = 0.;
	    ipick[it] = 0;
	}
	for (is = 0; is < ns; is++) {
	    sf_read (trace,sizeof(float),nt,in);

	    for (it = 0; it < nt; it++) {
		t = trace[it]*trace[it];
		if (t > ampl[it]) {
		    ampl[it] = t;
		    ipick[it] = is;
		}
	    }
	}
	    
	/* normalize */
	asum = 0.;
	for (it = 0; it < nt; it++) {
	    ampl[it] *= ampl[it];
	    asum += ampl[it];
	}
	for (it = 0; it < nt; it++) {
	    t = ampl[it]/asum;
	    diag[it] = 2.*eps + t;
	    pick[it] = t*(s0+ipick[it]*ds);
	}
	/* boundary conditions */
	diag[0] -= eps;
	diag[nt-1] -= eps;

	if (ix > 0) { /* lateral smoothing */
	    for (it = 0; it < nt; it++) {
		diag[it] += lam;
		pick[it] += lam*pick[it];
	    }
	}
	
	tridiagonal_define (slv,diag,offd);
	tridiagonal_solve (slv, pick);    
	sf_write (pick,sizeof(float),nt,out);
    }

    exit (0);
}


