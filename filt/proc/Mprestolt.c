/* Prestack Stolt modeling/migration.

Takes: < input.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"
#include "cosft.h"

int main(int argc, char* argv[])
{
    bool inv, stack, depth;
    int nt,nw,nx,ny,nh,nf, it,iw,ix,iy,ih;
    float dw,dx,dy,dh, x,y,h,xh, vel, w0, wh, w2, sq;
    float *str, *trace2, *trace, *keep=NULL;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* y: modeling, n: migration */
 
    if (!sf_getbool("depth",&depth)) depth=false;
    /* y: depth migration, n: time migration */
 
    if (inv) { /* modelling */
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nx)) nx=1;
	if (!sf_histint(in,"n3",&ny)) ny=1;
	
	if (!sf_getint ("nh",&nh)) sf_error("Need nh=");
	/* number of offsets */

	if (!sf_getfloat ("dh",&dh)) sf_error("Need dh=");
	/* offset sampling */

	dh = 1./(sf_npfar(2*(nh-1))*dh);

	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"d3",&dy)) dy=dx;

	sf_putint(out,"n2",nh); sf_putfloat(out,"d2",dh);
	sf_putint(out,"n3",nx); sf_putfloat(out,"d3",dx);
	sf_putint(out,"n4",ny); sf_putfloat(out,"d4",dy);
	sf_putfloat(out,"o4",0.);
    } else { /* migration */
	if (!sf_getbool ("stack",&stack)) stack=true;
	/* if y: stack migrated image */

	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nh)) nh=1;
	if (!sf_histint(in,"n3",&nx)) nx=1;
	if (!sf_histint(in,"n4",&ny)) ny=1;

	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");
	if (!sf_histfloat(in,"d4",&dy)) dy=dx;

	if (stack) {
	    sf_putint(out,"n2",nx); sf_putfloat(out,"d2",dx);
	    sf_putint(out,"n3",ny); sf_putfloat(out,"d3",dy);
	    sf_putint(out,"n4",1);
	}
    }
    
    if (!sf_getfloat ("vel",&vel)) sf_error("Need vel=");
    /* constant velocity */

    if (!sf_histfloat(in,"o1",&w0)) w0=0.; 
    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");

    if (!sf_getint ("pad",&nw)) nw=nt;
    /* padding on the time axis */
    nw=sf_npfar(2*(nw-1));

    cosft_init(nw /* , w0, dw */);
    dw = SF_PI/(nw*dw);
    dh *= SF_PI;
    dx *= SF_PI;
    dy *= SF_PI;

    if (depth) {
	if (inv) {
	    dw *= vel * 0.5;
	} else {
	    dw *= 2./vel;
	}
    } else {
	dh *= vel * 0.5;
	dx *= vel * 0.5;
	dy *= vel * 0.5;
    }

    trace2 = sf_floatalloc(nw);
    trace = sf_floatalloc(nw);
    str = sf_floatalloc(nw);

    if (stack) keep = sf_floatalloc(nw);

/*  w2 = (/ (((iw-1)*dw)**2, iw = 1, nw) /) */

    if (!sf_getint("nf",&nf)) nf=4;
    /* [2,4,6,8] Interpolation accuracy */
    prefilter_init (nf, nw, 3*nw);

    for (iy = 0; iy < ny; iy++) {
	y = iy*dy;
	y *= y;
	for (ix = 0; ix < nx; ix++) {
	    x = ix*dx;
	    x *= x;

	    if (inv) {
		sf_read(trace,sizeof(float),nt,in);
		for (it=nt; it < nw; it++) { /* pad */
		    trace[it]=0.;
		}
		
		cosft_frw (trace,0,1);
	    } else if (stack) {
		for (it=0; it < nw; it++) {
		    keep[it] = 0.;
		}
	    }

	    for (ih = 0; ih < nh; ih++) {
		h = ih * dh;
		h *= h;
		xh = x*h;
		h += x + y;

		for (iw = 0; iw < nw; iw++) {
		    w2 = iw*dw;
		    w2 *= w2;

		    if (inv) { /* modeling */
			if (xh == 0.) {
			    str[iw] = sqrtf (w2 + h);
			} else { 
			    if (w2 == 0.) {
				str[iw] = -2.*dw;
			    } else {
				str[iw] = sqrtf (w2 + h + xh/w2);
			    }
			}
		    } else { /* migration */
			wh = w2-h;
			sq = wh*wh - 4.*xh;

			if (wh > 0. && sq > 0.) {
			    str[iw] = sqrtf(0.5*(wh + sqrtf (sq)));
			} else {
			    str[iw] = - 2.*dw;
			}
		    }
		}

		int1_init (str, 0., dw, nw, spline_int, nf, nw);

		if (!inv) {
		    sf_read(trace,sizeof(float),nt,in);
		    for (it=nt; it < nw; it++) { /* pad */
			trace[it]=0.;
		    }		
		    cosft_frw (trace,0,1);
		}

		prefilter_apply (nw, trace);
		int1_lop (false,false,nw,nw,trace,trace2);   

		if (inv || !stack) {
		    cosft_inv (trace2,0,1);
		    sf_write(trace2,sizeof(float),nt,out);
		} else {
		    for (iw=0; iw < nw; iw++) {
			keep[iw] += trace2[iw];
		    }
		}
	    } /* h */
	    if (!inv && stack) {
		cosft_inv (keep,0,1);
		sf_write(keep,sizeof(float),nt,out);
	    }
	} /* x */
    } /* y */

    exit (0);
}



