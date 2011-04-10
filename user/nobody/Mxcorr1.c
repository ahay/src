/* Cross-correlation analysis.

Takes: < data.rsf other=other.rsf > out.rsf xcorr=xcorr.rsf
*/

#include <rsf.h>

#include "xcorr.h"
#include "window1.h"
#include "stretch2.h"

int main(int argc, char* argv[]) 
{
    int n, nw, w, iw, i0=0, maxshift, i, i2, n2, nc;
    float dt,h, eps, lam;
    bool verb, taper;
    float **dat, **dat2, **win, **win2, *coord, *shift, *warp, *xc;
    map2 str;
    sf_file in, out, other, xcr;

    sf_init(argc,argv);
    in = sf_input("in");
    other = sf_input("other");
    out = sf_output("out");
    xcr = sf_output("xcorr");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");

    n2 = sf_leftsize(in,1);

    if (!sf_getint("nw",&nw)) sf_error ("Need nw=");
    /* number of windows */
    if (!sf_getint("w",&w)) sf_error ("Need w=");
    /* window size */
    if (!sf_getfloat("h",&h)) h=0.5*(w-1);
    /* window overlap */
    if (!sf_getint("maxshift",&maxshift)) maxshift=w/2; nc=2*maxshift-1;
    /* maximum correlation lag */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* vertical smoothness (for picking) */
    if (!sf_getfloat("lam",&lam)) lam=0.5;
    /* horizontal smoothness (for picking) */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getbool("taper",&taper)) taper=true;
    /* window tapering */

    dat = sf_floatalloc2(n,n2);
    dat2 = sf_floatalloc2(n,n2);
    win = sf_floatalloc2(w,n2);
    win2 = sf_floatalloc2(w,n2);
    coord = sf_floatalloc(nw);
    shift = sf_floatalloc(nw);
    warp = sf_floatalloc(n);
    xc = sf_floatalloc(nc);
  
    window1_init (w,nw,n,h,w-h);
    str = stretch2_init (n,1.,1.,nw,eps,lam);

    sf_floatread (dat[0],n*n2,in);
    sf_floatread (dat2[0],n*n2,other);
  
    sf_putint(xcr,"n1",nc);
    sf_putint(xcr,"n2",nw);
    sf_putfloat(xcr,"o2",(0.5*w+1.)*dt);
    sf_putfloat(xcr,"d2",(n-w)*dt/(nw-1.));
    sf_putfloat(xcr,"o1",-maxshift*dt);

    for (iw=0; iw < nw; iw++) {
	for (i2=0; i2 < n2; i2++) {
	    i0 = window1_apply(iw,dat[i2] ,taper,taper,win[i2] );
	    i0 = window1_apply(iw,dat2[i2],taper,taper,win2[i2]);
	}
	coord[iw] = i0 + 0.5*(w+1.);
	shift[iw] = xcorr (w,n2,win[0], win2[0],nc,xc);
	sf_floatwrite(xc,nc,xcr);
	if (verb) sf_warning("shift[%d]=%g",iw,shift[iw]);
    }
  
    stretch2_define (str, coord, false);
    stretch2_apply (str, shift, warp);

    for (i=0; i < n; i++) {
	warp[i] *= dt;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatwrite (warp,n,out);
    }


    exit(0);
}

/* 	$Id$	 */
