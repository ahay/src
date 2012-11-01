/* Multi-dimensional missing data interpolation. */

#include <rsf.h>

int main(int argc, char* argv[])
{
    int na,ia, niter, n,i;
    float a0, *mm, *zero;
    bool prec, *known;
    char *lagfile;
    sf_filter aa;
    sf_file in, out, filt, lag, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    /* filter for inverse model covariance */
    out = sf_output("out");

    n = sf_filesize(in);

    if (!sf_getbool("prec",&prec)) prec=true;
    /* If y, use preconditioning */
    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1=");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
 
    /* Get filter lags */
    if (NULL == (lagfile = sf_histstring(filt,"lag"))) {
	for (ia=0; ia < na; ia++) {
	    aa->lag[ia]=ia+1;
	}
    } else {
	lag = sf_input(lagfile);
	sf_intread(aa->lag,na,lag);
	sf_fileclose(lag);
    }

    /* Get filter values */
    sf_floatread(aa->flt,na,filt);
    sf_fileclose(filt);

    /* Normalize */
    for (ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    mm = sf_floatalloc(n);
    zero = sf_floatalloc(n);
    known = sf_boolalloc(n);

    /* specify known locations */
    mask = sf_input("mask");
    sf_floatread(mm,n,mask);
    sf_fileclose(mask);
	
    for (i=0; i < n; i++) {
	known[i] = (bool) (mm[i] != 0.0f);
	zero[i] = 0.0f;
    }

    /* Read input data */
    sf_floatread(mm,n,in);

    if (prec) {                          
	sf_mask_init(known);
	sf_polydiv_init(n, aa);
	sf_solver_prec(sf_mask_lop, sf_cgstep, sf_polydiv_lop, 
		       n, n, n, mm, mm, niter, 0., "end");
    } else {                             
	sf_helicon_init(aa);
	sf_solver (sf_helicon_lop, sf_cgstep, n, n, mm, zero, 
		   niter, "known", known, "x0", mm, "end");
    }

    sf_floatwrite(mm,n,out);

    exit (0);
}
