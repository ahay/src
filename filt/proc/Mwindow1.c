#include <rsf.h>

#include "window1.h"

int main(int argc, char* argv[])
{
    int n2, i2, i1, n, nw, w, iw, i0;
    float h, *dat, *win, *dat2;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint ("nw",&nw)) sf_error("Need nw=");
    if (!sf_getint ("w",&w)) sf_error("Need w=");
    if (!sf_getfloat ("h",&h)) h = (nw*w - n)/(nw-1.);

    sf_putint(out,"n2",nw);
    sf_putint(out,"n3",n2);

    dat = sf_floatalloc (n);
    win = sf_floatalloc (w);
    dat2 = sf_floatalloc (n);
  
    window1_init (w,nw,n,h);

    for (i2=0; i2 < n2; i2++) {
	sf_read (dat, sizeof(float), n, in);
	for (iw=0; iw < nw; iw++) {
	    i0 = window1_apply(iw,dat,(iw > 0),(iw < nw-1),win);
	    for (i1=0; i1 < i0; i1++) {
		dat2[i1] = 0.;
	    }
	    for (i1=i0; i1 < i0+w; i1++) {
		dat2[i1] = win[i1-i0];
	    }
	    for (i1=i0+w; i1 < n; i1++) {
		dat2[i1] = 0.;
	    }
	    sf_write (dat2,sizeof(float),n,out);
	}
    }

    exit (0);
}

