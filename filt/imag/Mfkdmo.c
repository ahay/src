/* Offset continuation by log-stretch F-K operator.

Takes < input.rsf > output.rsf

*/

#include <math.h>
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nw, nk, iw, ik, nh, ih;
    float complex *oper;
    float dw,dk, ow, k,w,h, h0, dh, eps,amp,phase;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nk)) sf_error("No n2= in input");

    if (!sf_histfloat (in,"o1",&ow)) sf_error("No o1= in input");
    if (!sf_histfloat (in,"d1",&dw)) sf_error("No d1= in input");
    if (!sf_histfloat (in,"d2",&dk)) sf_error("No d2= in input");

    if (!sf_getfloat("h",&h)) sf_error("Need h=");
    /* final offset */
    if (!sf_getint("nh",&nh)) nh=1;
    /* number of offset steps */
    if (!sf_getfloat("h0",&h0)) h0=0.;
    /* initial offset */
    dh = (h-h0)/nh;

    sf_putint(out,"n3",nh);
    sf_putfloat(out,"d3",dh*2.);
    sf_putfloat(out,"o3",(h0+dh)*2.);

    oper = sf_complexalloc (nw);
    oper[0] = 0.; /* skip odd-ball frequency */

    for (ih = 0; ih < nh; ih++) {
	h = h0 + (ih+1)*dh;

	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk;
  
	    for (iw = 1; iw < nw; iw++) {
		w = ow + iw*dw;

		if (fabsf (w) > FLT_EPSILON) {
		    eps = 2.*k*h/w;
		    eps = sqrtf (1+eps*eps);
		    amp = sqrtf(0.5*(1/eps+eps))*expf(0.5*(1-eps));
		    phase = 1-eps+logf(0.5*(1+eps));

		    oper[iw] = amp*cexpf(SF_PI*I*phase*w);
		} else {
		    oper[iw] = 0.;
		}
	    }

	    sf_complexwrite(oper,nw,out);
	}
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mfkdmo.c,v 1.5 2004/04/19 21:51:16 fomels Exp $	 */
