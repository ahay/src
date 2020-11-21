/* Local frequency estimation 

Takes: rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
*/
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

int main (int argc, char* argv[])
{
    int n1,n2, i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *num, *den, *phase, mean, d1;
    char key[6];
    bool hertz;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	n12 *= n[i];
    }
    n2 = n12/n1;

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

    trace = sf_floatalloc(n1);

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    phase = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("hertz",&hertz)) hertz=false;
    /* if y, convert output to Hertz */

    mean=0.;
    for (i=0, i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);

    num[0] = 0.;
    den[0] = 0.;
    i++;

	for (i1=1; i1 < n1-1; i1++) {
	    num[i] = 0.5*(trace[i1-1]+trace[i1+1]);
	    den[i] = trace[i1];
	    mean += den[i]*den[i];
	    i++;
	}
    num[n1-1] = 0.;
    den[n1-1] = 0.;
    i++;
    }
    mean = sqrtf(n12/mean);
    
    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }

    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, phase);

    if (hertz) {
	/* convert to Hertz */    
	d1 = 1./(2.*SF_PI*d1);
	for (i=0; i < n12; i++) {
	    if (phase[i] < -1.) {
		phase[i]=SF_PI*d1;
	    } else if (phase[i] > 1.) {
		phase[i] = 0.;
	    } else {
		phase[i] = acosf(phase[i])*d1;
	    }
	}
    }

    sf_floatwrite(phase,n12,out);

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
