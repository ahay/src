/* Multicomponent data registration analysis. */
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


static float *coord, ***out, *rat2, *num, *den, g0, dg, o1, d1, o2, d2, g1, ***dataW;
static int n2g, ntr, n1, n2, ng, order;
static bool shift;
static sf_bands spl, spl2;

void warpscan_init(int m1     /* input trace length */, 
		   float o11  /* input origin */,
		   float d11  /* float increment */,
		   int m2     /* output trace length */,
		   float o21  /* output origin */,
		   float d21  /* output increment */,
		   int ng1    /* number of scanned "gamma" values */,
		   float g01  /* first gamma */,
		   float dg1  /* gamma increment */,
		   int ntr1   /* number of traces */, 
		   int order1 /* interpolation accuracy */, 
		   int dim    /* dimensionality */, 
		   int *m     /* data dimensions [dim] */, 
		   int *rect  /* smoothing radius [dim] */, 
		   int niter  /* number of iterations */,
		   bool shift1 /* shift instead of strech */,
		   bool verb  /* verbosity */)
/*< initialize >*/
{
    n1 = m1;
    o1 = o11;
    d1 = d11;
    o2 = o21;
    d2 = d21;
    n2 = m2;
    ntr = ntr1;
    ng = ng1;
    g0 = g01;
    dg = dg1;
    g1 = g0+(ng-1)*dg;
    n2g = n2*ng*ntr;
    order = order1;
    shift = shift1;

    coord = sf_floatalloc (n2); 
    out =   sf_floatalloc3 (n2,ng,ntr);

    out =   sf_floatalloc3 (n2,ng,ntr);
    dataW =   sf_floatalloc3 (n2,ng,ntr);

    rat2 = sf_floatalloc (n2g);
    num = sf_floatalloc (n2g);
    den = sf_floatalloc (n2g);

    spl = sf_spline_init (order, n1);
    spl2 = sf_spline_init (order, n1);

    sf_divn_init(dim, n2g, m, rect, niter, verb);
}

void warpscanw(bool cheb     /* use Chebyshev scan */,
	       bool sign     /* use signed similarity */,
	       bool renr      /* weight by reference energy */,
	       bool dend      /* weight by data energy */,
	       float** inp   /* input data [ntr][n1] */, 
	       float** oth   /* target data [ntr][n2] */,
	       float** rergy /* reference energy [ntr][n2]*/, 
	       float** dergy /* data energy [ntr][n1]*/,
	       float* rat1)
/*< scan >*/
{
    float doth, dout, g;
    int i1, i2, ig, i;

    doth = 0.;
    dout = 0.;
    for (i2=0; i2 < ntr; i2++) {
	sf_banded_solve (spl, inp[i2]);
	
	if (dend) {
	  sf_banded_solve (spl2, dergy[i2]);
	}


	for (i1=0; i1 < n2; i1++) {
	    doth += oth[i2][i1]*oth[i2][i1];
	}
	
	for (ig=0; ig < ng; ig++) {
	    if (cheb) {
		g = 0.5*(g0+g1)+0.5*(g1-g0)*cosf(ig*SF_PI/(ng-1));
	    } else {
		g = g0 + ig*dg;
	    }

	    for (i1=0; i1 < n2; i1++) {
		coord[i1] = shift? o2+i1*d2+g: (o2+i1*d2)*g;
	    }

	    sf_int1_init (coord, o1, d1, n1, sf_spline_int, order, n2, 0.0);
	    sf_int1_lop (false,false,n1,n2,inp[i2],out[i2][ig]);

	    if (dend) {
	      sf_int1_init (coord, o1, d1, n1, sf_spline_int, order, n2, 0.0);
	      sf_int1_lop (false,false,n1,n2,dergy[i2],dataW[i2][ig]);
	    }

	    for (i1=0; i1 < n2; i1++) {
	      dout += out[i2][ig][i1]*out[i2][ig][i1];
	    }
	}
    }
    doth = sqrtf(ntr*n2/doth);
    dout = sqrtf(n2g/dout);

    for (i2=0; i2 < ntr; i2++) {
	for (ig=0; ig < ng; ig++) {
	    for (i1=0; i1 < n2; i1++) {
		i = (i2*ng + ig)*n2+i1;
		  den[i] = out[i2][ig][i1]*dout;
		  num[i] = oth[i2][i1]*dout;
	    }
	}
    }

    sf_divn(num,den,rat1);
	
    for (i2=0; i2 < ntr; i2++) {
	for (ig=0; ig < ng; ig++) {
	    for (i1=0; i1 < n2; i1++) {
		i = (i2*ng+ig)*n2+i1;
		  num[i] = out[i2][ig][i1]*doth;
		  den[i] = oth[i2][i1]*doth;
	    }
	}
    }
    sf_divn(num,den,rat2);
    
    if (sign) {
	sf_divn_combine_sign(rat1,rat2,rat1);
    } else {
	sf_divn_combine(rat1,rat2,rat1);
    }


    if (renr && dend) {
      for (i2=0; i2 < ntr; i2++) {
	for (ig=0; ig < ng; ig++) {
	  for (i1=0; i1 < n2; i1++) {
	    i = (i2*ng+ig)*n2+i1;
	    rat1[i] *= dataW[i2][ig][i1]*rergy[i2][i1];
	  }
	}
      }
    }

    if (renr) {
      for (i2=0; i2 < ntr; i2++) {
	for (ig=0; ig < ng; ig++) {
	  for (i1=0; i1 < n2; i1++) {
	    i = (i2*ng+ig)*n2+i1;
	    rat1[i] *= rergy[i2][i1];
	  }
	}
      }
    }


}

/* 	$Id: Mwarpscan.c 744 2004-08-17 18:46:07Z fomels $	 */
