/* 2-D anisotropic diffusion.

Takes: < in.rsf > out.rsf
*/

#include <rsf.h>

#include "impl2.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i3;
    float r1, r2, tau, pclip, **dat;
    bool up;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getfloat("rect1",&r1)) sf_error("Need rect1=");
    /* vertical smoothing */
    if (!sf_getfloat("rect2",&r2)) sf_error("Need rect2=");
    /* horizontal smoothing */
    if (!sf_getfloat("tau",&tau)) tau=0.1;
    /* smoothing time */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */
    if (!sf_getbool ("up",&up)) up=false;
    /* smoothing style */

    impl2_init (r1, r2, n1, n2, tau, pclip, up);

    dat = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	sf_read (dat[0],sizeof(float),n1*n2,in);
	impl2_apply (dat,true,false);
	sf_write (dat[0],sizeof(float),n1*n2,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mimpl2.c,v 1.3 2004/04/09 13:17:10 fomels Exp $	 */

