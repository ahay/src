#include <rsf.h>

#include "tridiagonal.h"

int main(int argc, char* argv[])
{
    tris slv;
    int nt, ns, nx, ix, is, it, iter, niter;
    float **slice, *pick, *pick0, *ampl, *diag, *offd;
    float s0, ds, eps, lam, t, asum;
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
    if (!sf_getint("niter",&niter)) niter=1;

    slice = sf_floatalloc2(nt,ns);
    pick = sf_floatalloc(nt);
    pick0 = sf_floatalloc(nt);
    offd = sf_floatalloc(nt);
    diag = sf_floatalloc(nt);
    ampl = sf_floatalloc(nt);

    slv = tridiagonal_init (nt);

    for (it=0; it < nt; it++) {
	offd[it] = -eps;
    }

    for (ix = 0; ix < nx; ix++) {
	sf_read (slice[0],sizeof(float),nt*ns,in);

	if (ix > 0) {
	    for (it = 0; it < nt; it++) { /* previous trace */
		pick0[it] = pick[it];
	    }
	}

	for (iter=0; iter < niter; iter++) {
	    if (0==iter) { /* pick blind maximum */
		for (it = 0; it < nt; it++) {
		    ampl[it] = 0.;
		    pick[it] = s0+0.5*(ns-1)*ds;
		}

		for (is = 0; is < ns; is++) {
		    for (it = 0; it < nt; it++) {
			t = slice[is][it];
			t *= t;
			if (t > ampl[it]) {
			    ampl[it] = t;
			    pick[it] = s0+is*ds;
			}
		    }
		}
	    } else {
		for (it = 0; it < nt; it++) {
		    /* nearest neighbor interpolation */
		    is = 0.5+(pick[it]-s0)/ds;
		    if (is < 0) {
			is=0;
			pick[it]=s0;
		    } else if (is > ns-1) {
			is=ns-1;
			pick[it] = s0 + (ns-1)*ds;
		    }
		    t = slice[is][it];
		    ampl[it] = t*t;
		}
	    }

	    /* normalize amplitudes */
	    asum = 0.;
	    for (it = 0; it < nt; it++) {
		ampl[it] *= ampl[it];
		asum += ampl[it];
	    }
	    for (it = 0; it < nt; it++) {
		t = ampl[it]/asum;
		diag[it] = 2.*eps + t;
		pick[it] *= t;
	    }
	    /* boundary conditions */
	    diag[0] -= eps;
	    diag[nt-1] -= eps;
	    
	    if (ix > 0) { /* lateral smoothing */
		for (it = 0; it < nt; it++) {
		    diag[it] += lam;
		    pick[it] += lam*pick0[it];
		}
	    }
	    
	    tridiagonal_define (slv,diag,offd);
	    tridiagonal_solve (slv, pick);    
	} /* iterations */

	sf_write (pick,sizeof(float),nt,out);	
    }

    exit (0);
}


