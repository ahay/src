/* Forward interpolation in 2-D slices.

Takes: < binned.rsf head=header.rsf > output.rsf
*/

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int2.h"
#include "interp.h"
#include "laplac2.h"

int main (int argc, char* argv[])
{
    int id, nk, nd, nm, nt, it, nx, ny, xkey, ykey, interp;
    float *mm, *dd, **xy, *hdr, x0, y0, dx, dy;
    char *xk, *yk;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nx)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&ny)) sf_error("Need n2= in in");
    nt = sf_leftsize(in,2);
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (NULL != (xk = sf_getstring("xk"))) {
	/* x key name */
	xkey = sf_segykey(xk);
    }  else if (!sf_getint("xkey",&xkey)) {
	/* x key number (if no xk), default is sx */
	xkey = sf_segykey("sx");
    }
    if (NULL != (yk = sf_getstring("yk"))) {
	/* y key name */
	ykey = sf_segykey(yk);
    }  else if (!sf_getint("ykey",&ykey)) {
	/* y key number (if no yk), default is sy */
	ykey = sf_segykey("sy");
    }

    /* create coordinates */
    head = sf_input("head");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&nd)) sf_error("Need n2= in head");
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");

    xy = sf_floatalloc2(2,nd);
    hdr = sf_floatalloc(nk);

    for (id=0; id<nd; id++) {	
	sf_read (hdr,sizeof(float),nk,head);
	xy[id][0] = hdr[xkey];
	xy[id][1] = hdr[ykey];
    }

    sf_fileclose (head);

    sf_putint(out,"n1",nd);
    sf_putint(out,"n2",1);

    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&y0)) sf_error("No o2= in input");
    
    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* [1,2] interpolation method, 1: nearest neighbor, 2: bi-linear */

    switch (interp) {
	case 1:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, bin_int, 1, nd);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, lin_int, 2, nd);
	    sf_warning("Using linear interpolation");
	    break;
	case 3:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_read (mm,sizeof(float),nm,in);

	int2_lop(false,false,nm,nd,mm,dd);

	sf_write (dd,sizeof(float),nd,out);
    }

    exit(0);
}

/* 	$Id: Mextract.c,v 1.1 2004/02/14 07:01:42 fomels Exp $	 */
