/* Automatic gain control by shaping regularization. */
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

#include <math.h>
#include <rsf.h>

int main (int argc, char* argv[]) 
{
    bool verb;
    int n[SF_MAX_DIM], rect[SF_MAX_DIM];
    int n1, n2, i2, niter, dim, dim1, i1, i;
    char key[7];
    float *data, *gain, *igain, *one, eps;
    sf_file in, out, fgain;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,7,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]= i? 1: 125;
        /*( rect#=(125,1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i;
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);
    gain = sf_floatalloc (n1);
    igain = sf_floatalloc (n1);
    one = sf_floatalloc (n1);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */

    if (NULL != sf_getstring("gain")) {
	/* output gain file (optional) */
	fgain = sf_output("gain");
    } else {
	fgain = NULL;
    }

    sf_divn_init(dim1+1,n1,n,rect,niter,verb);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	for (i1=0; i1 < n1; i1++) {
	    gain[i1] = fabsf(data[i1]);
	    one[i1] = 1.0f;
	}

	sf_divne(one,gain,igain,eps);

	for (i1=0; i1 < n1; i1++) {
	    data[i1] *= igain[i1];
	}
	
	if (NULL != fgain) 
	    sf_floatwrite(igain,n1,fgain);

	sf_floatwrite(data,n1,out);
    }

    exit(0);
}
