#include <math.h>

#include <rsf.h>

#include "cburg.h"
#include "shotfill.h"

int main(int argc, char* argv[])
{
    int ns, nh, nw, iw, ih;
    float ds, h0, dh, w0, dw, w, eps;
    float complex *s, **ss, a[3];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ns)) sf_error("No n2= in input");
    
    if (ns != 2) sf_error("Need n2=2");
    /* change later */

    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    sf_putint(out,"n3",1);
    sf_putint(out,"n2",nw);
    sf_putfloat(out,"d2",dw);
    sf_putfloat(out,"o2",w0);

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* regularization parameter */

    ss = sf_complexalloc2(nh,2);
    s = sf_complexalloc(nh);

    shotfill_init(nh,h0,dh,ds*0.5);
    cburg_init(nh,2,3);

    for (iw=0; iw < nw; iw++) {
	w = w0 + iw*dw;

	sf_read (ss[0],sizeof(float complex),nh*2,in);

	if(fabsf(w) < dw) {
	    for (ih=0; ih < nh; ih++) {
		s[ih] = 0.;
	    }
	} else {
	    cburg_apply (ss, a); 
	    a[0] *= eps;
	    a[1] *= eps;
	    a[2] *= eps;
	    shotfill_define(w, a);
	    shotfill_apply(ss[0],ss[1],s);
	}

	sf_write(s,sizeof(float complex),nh,out);
    }
  
    exit(0);
}


