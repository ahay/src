/* Convection filter. */
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
    bool verb, **mm;
    int niter, order, n1, n2, n12, i1, i2, iw, nw, rect[2], m[2];
    float **data, **rhs, ***conv, ***lhs, dif, mean;
    sf_file in, flt, mask;

    sf_init(argc, argv);
    in = sf_input("in");
    flt = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2= in input");
    n12 = n1*n2;
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */
    nw = 2*order;

    sf_shiftdim(in, flt, 2);
    sf_putint(flt, "n2", nw);

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* smoothing */

    if (!sf_getint("niter",&niter)) niter=100*order;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (NULL != sf_getstring("mask")) {
	mm = sf_boolalloc2(n12,2);
	mask = sf_input("mask");
    } else {
	mm = (bool**) sf_alloc(2,sizeof(bool*));
	mm[0] = mm[1] = NULL;
	mask = NULL;
    }

    data = sf_floatalloc2(n1,n2);
    rhs = sf_floatalloc2(n1,n2);
  
    conv = sf_floatalloc3(n1,n2,nw);
    lhs = sf_floatalloc3(n1,n2,nw);

    m[0] = n1;
    m[1] = n2;

    sf_multidivn_init(nw, 2, n12, m, rect, **lhs, NULL, verb);

    if (NULL != mask) {
	sf_floatread(data[0],n12,mask);
	sf_mask32 (false, order, 1, 1, n1, n2, 1, data[0], mm);
    }

    sf_floatread(data[0],n12,in);

    mean = 0.;
    for (i2=0; i2 < n2-1; i2++) {
	for (i1=0; i1 < order; i1++) {
	    rhs[i2][i1] = 0.0f;
	    for (iw=0; iw < nw; iw++) {
		lhs[iw][i2][i1] = 0.0f;
	    }
	}
	for (i1=order; i1 < n1-order; i1++) {
	    dif = data[i2+1][i1] - data[i2][i1];
	    rhs[i2][i1] = dif;
	    for (iw=1; iw <= order; iw++) {
		lhs[2*iw-2][i2][i1] = dif - data[i2+1][i1+iw] + data[i2][i1-iw];
		lhs[2*iw-1][i2][i1] = dif - data[i2+1][i1-iw] + data[i2][i1+iw];
	    }
	    for (iw=0; iw < nw; iw++) {
		mean += lhs[iw][i2][i1]*lhs[iw][i2][i1];
	    }
	}
	for (i1=n1-order; i1 < n1; i1++) {
	    rhs[i2][i1] = 0.0f;
	    for (iw=0; iw < nw; iw++) {
		lhs[iw][i2][i1] = 0.0f;
	    }
	}
    }
    
    for (i1=0; i1 < n1; i1++) {
	rhs[n2-1][i1] = 0.0f;
	for (iw=0; iw < nw; iw++) {
	    lhs[iw][n2-1][i1] = 0.0f;
	}
    }
    
    mean = sqrtf (mean/(n12*2));
    for (i2=0; i2 < n2-1; i2++) {
	for (i1=order; i1 < n1-order; i1++) {
	    if (NULL != mask && mm[0][i1+n1*i2]) {
		rhs[i2][i1] = 0.0f;
		for (iw=0; iw < nw; iw++) {
		    lhs[iw][i2][i1] = 0.0f;
		}
	    } else {
		rhs[i2][i1] /= mean;
		for (iw=0; iw < nw; iw++) {
		    lhs[iw][i2][i1] /= mean;
		}
	    }
	}
    }

    sf_multidivn (*rhs,**conv,niter);

    /* transpose */
    for (i2=0; i2 < n2; i2++) {
	for (iw=0; iw < nw; iw++) {
	    sf_floatwrite(conv[iw][i2],n1,flt);
	}
    }
  
    exit(0);
}
  

  
