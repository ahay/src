#include <math.h>

#include <rsf.h>

#include "shotfill.h"
/*
#include "cburg.h"
*/

int main(int argc, char* argv[])
{
    int ns, nh, nw, iw, ih, is;
    bool sign;
    float ds, h0, dh, w0, dw, w, eps;
    float complex *s, **ss; /* a[3]; */
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ns)) sf_error("No n2= in input");

    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation */

    sf_putint(out,"n2",2*ns-1);
    sf_putfloat(out,"d2",0.5*ds);

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    ss = sf_complexalloc2(nh,ns);
    s = sf_complexalloc(nh);

    shotfill_init(nh,h0,dh,sign? ds: -ds, eps);
/*    cburg_init(nh,ns,3); */

    for (iw=0; iw < nw; iw++) {
	w = w0 + iw*dw;

	sf_read (ss[0],sizeof(float complex),nh*ns,in);

	if(fabsf(w) < dw) { /* dc */
	    for (ih=0; ih < nh; ih++) {
		s[ih] = 0.;
	    }
	    sf_write(s,sizeof(float complex),nh,out);
	    for (is=1; is < ns; is++) {
		sf_write(s,sizeof(float complex),nh,out);
		sf_write(s,sizeof(float complex),nh,out);
	    }
	    continue;
	}

/*
	cburg_apply (ss, a);
	a[0] *= eps;
	a[1] *= eps;
	a[2] *= eps;
*/

	shotfill_define(w);

	sf_write(ss[0],sizeof(float complex),nh,out);

	for (is=1; is < ns; is++) {
	    shotfill_apply(ss[is-1],ss[is],s);

	    sf_write(s,sizeof(float complex),nh,out);
	    sf_write(ss[is],sizeof(float complex),nh,out);
	} /* s */
    } /* w */
  
    sf_close();
    exit(0);
}

/* 	$Id: Minfill.c,v 1.4 2004/03/26 03:30:36 fomels Exp $	 */

