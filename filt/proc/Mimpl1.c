/* 1-D anisotropic diffusion.

Takes: < in.rsf > out.rsf
*/

#include <rsf.h>

#include "impl1.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2;
    float r, tau, pclip, *dat;
    bool up;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getfloat("rect1",&r)) sf_error("Need rect1=");
    /* smoothing radius */
    if (!sf_getfloat("tau",&tau)) tau=0.1;
    /* smoothing time */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */
    if (!sf_getbool ("up",&up)) up=false;
    /* smoothing style */

    impl1_init (r, n1, tau, pclip, up);

    dat = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (dat,n1,in);
	impl1_apply (dat);
	sf_floatwrite (dat,n1,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mimpl1.c,v 1.3 2004/04/19 21:51:46 fomels Exp $	 */

