#include <rsf.h>

#include "xcorr.h"
#include "window1.h"
#include "stretch2.h"

int main(int argc, char* argv[]) 
{
    int n, nw, w, iw, i0, maxshift, i, order, i2, n2;
    float dt,h, eps, lam;
    bool verb;
    float **dat, **dat2, **win, **win2, *coord, *shift, *warp;
    map2 str;
    sf_file in, out, other;

    sf_init(argc,argv);
    in = sf_input("in");
    other = sf_input("other");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");

    n2 = sf_leftsize(in,1);

    if (!sf_getint("nw",&nw)) sf_error ("Need nw=");
    if (!sf_getint("w",&w)) sf_error ("Need w=");
    if (!sf_getfloat("h",&h)) h=0.5*(w-1);
    if (!sf_getint("maxshift",&maxshift)) maxshift=w/2;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getfloat("lam",&lam)) lam=0.5;
    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getint("order",&order)) order=nw;

    dat = sf_floatalloc2(n,n2);
    dat2 = sf_floatalloc2(n,n2);
    win = sf_floatalloc2(w,n2);
    win2 = sf_floatalloc2(w,n2);
    coord = sf_floatalloc(nw);
    shift = sf_floatalloc(nw);
    warp = sf_floatalloc(n);
  
    window1_init (w,nw,n,h);
    xcorr_init (w, n2, maxshift);
    str = stretch2_init (n,1.,1.,nw,eps,lam);

    sf_read (dat[0],sizeof(float),n*n2,in);
    sf_read (dat2[0],sizeof(float),n*n2,other);
  
    for (iw=0; iw < nw; iw++) {
	for (i2=0; i2 < n2; i2++) {
	    i0 = window1_apply(iw,dat[i2] ,true,true,win[i2] );
	    i0 = window1_apply(iw,dat2[i2],true,true,win2[i2]);
	}
	coord[iw] = i0 + 0.5*(w+1.);
	shift[iw] = xcorr (win[0], win2[0]);
	if (verb) sf_warning("shift[%d]=%g",iw,shift[iw]);
    }
  
    stretch2_define (str, coord, false);
    stretch2_apply (str, shift, warp);

    for (i=0; i < n; i++) {
	warp[i] *= dt;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_write (warp,sizeof(float),n,out);
    }

    exit(0);
}
