/* Automatic picking from semblance-like panels using shaping regularization.

Takes: < semblance.rsf [ampl=ampl.rsf] > pick.rsf
*/
#include <math.h>
#include <float.h>

#include <rsf.h>

#include "divn.h"

int main(int argc, char* argv[])
{
    int dim, n[SF_MAX_DIM], rect[SF_MAX_DIM], nd;
    int ns, nt, nm, nx, i, ix, is, it, niter, imax;
    float **slice, *pick0, *pick, *ampl;
    float s0, ds, t, asum, ai, amax, ap, am, num, den, dx;
    char key[6];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    dim = sf_filedims (in,n);
    if (dim < 2) sf_error("Need at least two dimensions");

    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	} else {
	    rect[i]=1;
	}
    }
    nt = n[0];
    ns = n[1];
    nm = nd/ns;
    nx = nm/nt;

    if (!sf_histfloat(in,"o2",&s0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");
    sf_putint(out,"n2",1);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    divn_init(dim,nd,n,rect,niter);

    slice = sf_floatalloc2(nt,ns);

    pick0 = sf_floatalloc(nm);
    pick = sf_floatalloc(nm);
    ampl = sf_floatalloc(nm);

    for (ix = 0; ix < nx; ix++) {
	sf_floatread (slice[0],nt*ns,in);

	/* pick blind maximum */
	for (it = 0; it < nt; it++) {
	    i = ix*nt+it;

	    imax = 0;
	    amax = 0.;
	    for (is = 0; is < ns; is++) {
		t = slice[is][it];
		if (t > amax) {
		    amax = t;
		    imax = is;
		}
	    }

	    /* quadratic interpolation for sub-pixel accuracy */
	    if (imax > 0 && imax < ns-1) {
		am = slice[imax-1][it];
		ap = slice[imax+1][it];
		num = 0.5*(am-ap);
		den = am+ap-2.*amax;
		dx = num*den/(den*den + FLT_EPSILON);
		if (fabsf(dx) >= 1.) dx=0.;

		ampl[i] = amax - 0.5*dx*dx*den;
		pick0[i] = s0+(imax+dx)*ds;
	    } else {
		ampl[i] = amax;
		pick0[i] = s0+imax*ds;
	    }
	}
    }
    
    /* normalize amplitudes */
    asum = 0.;
    for (i = 0; i < nm; i++) {
	ai = ampl[i];
	asum += ai*ai;
    }
    asum = sqrtf (asum/nm);

    for(i=0; i < nm; i++) {
	ampl[i] /= asum;
	pick0[i] *= ampl[i];
    }

    divn(pick0,ampl,pick);

    sf_floatwrite (pick,nm,out);	
    
    sf_close();
    exit (0);
}

/* 	$Id: Mblindpick.c,v 1.8 2004/05/07 03:40:04 fomels Exp $	 */

