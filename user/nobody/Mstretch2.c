/* Smooth inverse interpolation.

Takes: < data.rsf head=header.rsf > stretched.rsf
*/

#include <rsf.h>

#include "stretch2.h"

int main(int argc, char* argv[]) 
{
    int nx, nd;
    float x0, dx, eps, lam;
    float *dat, *coord, *mod;
    map2 str;
    sf_file in, out, head;

    sf_init(argc,argv);
    in = sf_input("in");
    head = sf_input("head");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) sf_error("No n1= in input");

    dat = sf_floatalloc(nd);
    coord = sf_floatalloc(nd);

    sf_floatread (dat,nd,in);
    sf_floatread (coord,nd,head);
    
    if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
    /* number of output samples */
    if (!sf_getfloat("x0",&x0)) x0=coord[0];
    /* output origin */
    if (!sf_getfloat("dx",&dx)) dx=(coord[nd-1]-x0)/(nd-1);
    /* output sampling */

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putfloat(out,"o1",x0);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=0.5;
    /* horizontal smoothness */

    mod = sf_floatalloc(nx);

    str = stretch2_init (nx,x0,dx,nd,eps,lam);
    stretch2_define (str, coord, true);
    stretch2_apply (str, dat, mod);

    sf_floatwrite (mod,nx,out);


    exit(0);
}

/* 	$Id$	 */

