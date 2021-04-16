/* Inverse 1-D warping. 

September 2012 program of the month:
http://ahay.org/blog/2012/09/03/program-of-the-month-sfiwarp/
*/
/*
  Copyright (C) 2004 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_map4 mo;
    bool inv, each=true;
    int i, nt, n1, i2, n2, nw;
    float o1, d1, t0, dt, eps;
    sf_complex *ctrace, *ctrace2;
    float *trace, *str, *trace2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

    if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (inv) {
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	
	if (!sf_getint("n1",&n1)) n1=nt; /* output samples - for inv=y */
	if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
	/*( d1=1 output sampling - for inv=y )*/
	if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
	/*( o1=0 output origin - for inv=y )*/ 

	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
    } else {
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (!sf_histint(warp,"n1",&nt)) sf_error("No n1= in warp");
	if (!sf_histfloat(warp,"d1",&dt)) dt=d1;
	if (!sf_histfloat(warp,"o1",&t0)) t0=o1;
	
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
    }

    n2 = sf_leftsize(in,1);
    nw = sf_leftsize(warp,1);
    if (1 == nw) {
	each = false;
    } else if (n2 != nw) {
	sf_error("Need %d traces in warp, got %d",n2,nw);
    } 
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    trace2 = sf_floatalloc(n1);

    mo = sf_stretch4_init (n1, o1, d1, nt, eps);

    if (SF_COMPLEX == sf_gettype(in)) {
	ctrace = sf_complexalloc(nt);
	ctrace2 = sf_complexalloc(n1);
    } else {
	ctrace = ctrace2 = NULL;
    }

    for (i2=0; i2 < n2; i2++) {
	if (each || 0==i2) {
	    sf_floatread(str,nt,warp);
	    sf_stretch4_define (mo,str,false);
	}

	if (inv) {
	    if (SF_COMPLEX == sf_gettype(in)) {
		sf_complexread(ctrace,nt,in);
		for (i=0; i < nt; i++) {
		    trace[i] = crealf(ctrace[i]);
		}
		sf_stretch4_apply (false,mo,trace,trace2);
		for (i=0; i < n1; i++) {
		    ctrace2[i] = sf_cmplx(trace2[i],0.0f);
		}
		for (i=0; i < nt; i++) {
		    trace[i] = cimagf(ctrace[i]);
		}
		sf_stretch4_apply (false,mo,trace,trace2);
		for (i=0; i < n1; i++) {
#ifdef SF_HAS_COMPLEX_H
		    ctrace2[i] += sf_cmplx(0.0f,trace2[i]);
#else
		    ctrace2[i] = sf_cadd(ctrace2[i],sf_cmplx(0.0f,trace2[i]));
#endif
		}
		sf_complexwrite (ctrace2,n1,out);
	    } else {
		sf_floatread(trace,nt,in);
		sf_stretch4_apply (false,mo,trace,trace2);
		sf_floatwrite (trace2,n1,out);
	    }
	} else {
	    if (SF_COMPLEX == sf_gettype(in)) {
		sf_complexread (ctrace2,n1,out);
		for (i=0; i < n1; i++) {
		    trace2[i] = crealf(ctrace2[i]);
		}
		sf_stretch4_invert (false,mo,trace,trace2);
		for (i=0; i < nt; i++) {
		    ctrace[i] = sf_cmplx(trace[i],0.0f);
		}
		for (i=0; i < n1; i++) {
		    trace2[i] = cimagf(ctrace2[i]);
		}
		sf_stretch4_invert (false,mo,trace,trace2);
		for (i=0; i < nt; i++) {
#ifdef SF_HAS_COMPLEX_H
		    ctrace[i] += sf_cmplx(0.0f,trace[i]);
#else
		    ctrace[i] = sf_cadd(ctrace[i],sf_cmplx(0.0f,trace[i]));
#endif
		}
		sf_complexwrite(ctrace,nt,in);
	    } else {
		sf_floatread(trace2,n1,in);
		sf_stretch4_invert (false,mo,trace,trace2);
		sf_floatwrite (trace,nt,out);
	    }
	}
    }

    exit(0);
}
