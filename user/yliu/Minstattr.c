/* Estimate of instantaneous attributes. */
/*
  Copyright (C) 2011 Jilin University
  
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
#include <math.h>

#include <rsf.h>


int main (int argc, char* argv[])
{
    int nh, n1, n2, i1,i2, i, dim, n[SF_MAX_DIM];
    float *trace, *hilb, *dtrace, *dhilb, *attr, *dert=NULL, c, d1;
    char *type;
    bool hertz, hires, der2, verb;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n2 = 1;
    for (i=0; i < dim; i++) {
	n2 *= n[i];
    }
    n2 /= n1;

    if (NULL == (type=sf_getstring("type"))) type="amplitude";
    /* [amplitude,phase,frequency] attribute type, the default is amplitude */

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getint("order",&nh)) nh=10;
    /* Hilbert transformer order */

    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    if (!sf_getbool("hertz",&hertz)) hertz=true;
    /* if y, convert output to Hertz */

    if (!sf_getbool("hires",&hires)) hires=false;
    /* if y, compute high resolution instantaneous attributes */

    if (hires) {
	if (!sf_getbool("der2",&der2)) der2=false;
	/* if y, compute 2nd-order derivative to use with hires=y */
    }

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (hires) dert = sf_floatalloc(n1);
    trace = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);
    dtrace = sf_floatalloc(n1);
    dhilb = sf_floatalloc(n1);
    attr = sf_floatalloc(n1);

    sf_hilbert_init(n1, nh, c);
    if ('f'==type[0] || hires) sf_deriv_init(n1, nh, c);
    
    d1 = 1./(2.*SF_PI*d1);
    for (i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("trace %d of %d;",i2+1,n2);

	sf_floatread(trace,n1,in);

	for (i1=0; i1 < n1; i1++) {
	    hilb[i1] = 0.;
	    dtrace[i1] = 0.;
	    dhilb[i1] = 0.;
	    if (hires) dert[i1] = 0.;
	}
	if (hires) {
	    sf_deriv(trace,dert);
	    if (der2) {
		sf_deriv(dert,trace);
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] = dert[i1];
		}
	    }
	}
	sf_hilbert(trace,hilb);
	switch(type[0]) {
	    case 'a':
		for (i1=0; i1 < n1; i1++) {
		    attr[i1] = hypotf(trace[i1],hilb[i1]);
		}
		break;

	    case 'p':
		for (i1=0; i1 < n1; i1++) {
		    attr[i1] = atanf(hilb[i1]/
				     (trace[i1]+FLT_EPSILON))*90./SF_PI;
		    if (hires) {
			if (der2) {
			    attr[i1] *= -1.;
			} 
		    }
		}

		break;

	    case 'f':
		sf_deriv(trace,dtrace);
		sf_deriv(hilb,dhilb);
		for (i1=0; i1 < n1; i1++) {
		    attr[i1] = (trace[i1]*dhilb[i1]-dtrace[i1]*hilb[i1])/
			(trace[i1]*trace[i1]+hilb[i1]*hilb[i1]+FLT_EPSILON);
		}

		if (hertz) {
		    /* convert to Hertz */    
		    for (i1=0; i1 < n1; i1++) {
			attr[i1] *= d1;
		    }
		}
		break;

	    default:
		sf_error("Unknown attribute type=%c",type[0]);
		break;
	}
	sf_floatwrite(attr,n1,out);
    }
    if (verb) sf_warning(".");   
    exit(0);
}

/* 	$Id$	 */
