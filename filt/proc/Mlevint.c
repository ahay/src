/* Leveler inverse interpolation in 1-D. */

#include <rsf.h>

#include "levint.h"

int main(int argc, char* argv[])
{
    int nd, m1, na, nr, niter;
    float *rr, *dd, *coord, o1, d1, eps;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* create model */
    if (!sf_getint ("nx",&m1)) sf_error("Need n1=");
    /* number of bins */
    sf_putint(out,"n1",m1);

    if (!sf_getfloat("x0",&o1)) sf_error("Need o1=");
    /* grid origin */
    sf_putfloat (out,"o1",o1);

    if (!sf_getfloat("dx",&d1)) sf_error("Need d1=");
    /* grid sampling */
    sf_putfloat (out,"d1",d1);

    if (!sf_getint("niter",&niter)) niter=1+m1*3/2; niter *= 2;
    /* number of conjugate-gradient iterations */
    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* regularization parameter */

    /* create filter */
    if (!sf_getint("na",&na)) na=3;
    nr = m1 + na;

    rr = sf_floatalloc(nr);

    if (!sf_histint(in,"n1",&nd)) nd=1;
    coord = sf_floatalloc(nd);
    dd = sf_floatalloc(nd);

    head = sf_input("head");
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");
    sf_floatread (coord,nd,head);
    sf_fileclose (head);

    sf_floatread (dd,nd,in);

    levint1 (niter, m1, nr, nd, coord, dd, o1, d1, rr, eps);
    
    sf_floatwrite (rr,m1,out);

    sf_close();
    exit(0);
}





