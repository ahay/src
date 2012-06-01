/* Phase scan */
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
#include "phasescan.h"


static float *coord, ***out, *rat2, *num, *den, o1, d1, o2, d2, angle0, dangle;
static int n2g, ntr, n1, n2, na, order;
static sf_bands spl;

void phasescan_init(int m1     /* input trace length */, 
		   float o11  /* input origin */,
		   float d11  /* float increment */,
		   int m2     /* output trace length */,
		   float o21  /* output origin */,
		   float d21  /* output increment */,
                   int na1    /* number of scanned angle */,
		   float a01  /* first angle*/,
                   float da1  /* angle increment */,
		   int ntr1   /* number of traces */, 
		   int order1 /* interpolation accuracy */, 
		   int dim    /* dimensionality */, 
		   int *m     /* data dimensions [dim] */, 
		   int *rect  /* smoothing radius [dim] */, 
		   int niter  /* number of iterations */,
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
    na = na1; /*angle number*/
    angle0 = a01; /*first angle*/
    dangle = da1;/* angle increment*/
    n2g = n2*na*ntr;
    order = order1;

    coord = sf_floatalloc (n2); 
    out =   sf_floatalloc3 (n2,na,ntr);

    rat2 = sf_floatalloc (n2g);
    num = sf_floatalloc (n2g);
    den = sf_floatalloc (n2g);

    spl = sf_spline_init (order, n1);     
    sf_divn_init(dim, n2g, m, rect, niter, verb);
}

void phasescan(float** inp /* input data [ntr][n1] */, 
	      float** oth /* target data [ntr][n2] */,
	      float* rat1)
/*< scan >*/
{
    float doth, dout; /* doth and dout are for normalization*/
    int i1, i2, ia, i;/*i1 is time index, i2 is trace index*/
    float *htrace;
    htrace =  sf_floatalloc(n2);
    float theta,expc,exps,tmp;

    doth = 0.;
    dout = 0.;
    angle0 *= 1/180.0*SF_PI; 
    dangle *= 1/180.0*SF_PI; 
    sf_hilbert_init(n1,n2,1.0);
    for (i2=0; i2 < ntr; i2++) {
//	sf_banded_solve (spl, inp[i2]);

	for (i1=0; i1 < n2; i1++) {
	    doth += oth[i2][i1]*oth[i2][i1];
	}
        
	/*Hilbert Transform inp[i2] */
        sf_hilbert(inp[i2],htrace);
        
        /*rotate to angle theta */
        /*multiply by exp(i*theta*pi/180)*/
        for (ia=0; ia < na; ia++) {
            theta = angle0+ia*dangle;
            //sf_warning("theta=%g",theta);
            expc = cosf(theta);
            exps = sinf(theta);
            for (i1=0; i1 < n1; i1++){
                tmp = inp[i2][i1]*expc-htrace[i1]*exps;
                if(i1<n2) out[i2][ia][i1] = tmp;
            }
            for (i1=n1; i1 < n2; i1++) out[i2][ia][i1] = 0.0; 
            
            for (i1=0; i1 < n2; i1++) {
	        dout += out[i2][ia][i1]*out[i2][ia][i1];
    	    }
        }
     
      
	
    }

    doth = sqrtf(ntr*n2/doth);
    dout = sqrtf(n2g/dout);
    /* calculate local correlation*/
    for (i2=0; i2 < ntr; i2++) {
	for (ia=0; ia < na; ia++) {
	    for (i1=0; i1 < n2; i1++) {
		i = (i2*na + ia)*n2+i1;
		den[i] = out[i2][ia][i1]*dout;
		num[i] = oth[i2][i1]*dout;
	    }
	}
    }

    sf_divn(num,den,rat1);
	
    for (i2=0; i2 < ntr; i2++) {
	for (ia=0; ia < na; ia++) {
	    for (i1=0; i1 < n2; i1++) {
		i = (i2*na+ia)*n2+i1;
		num[i] = out[i2][ia][i1]*doth;
		den[i] = oth[i2][i1]*doth;
	    }
	}
    }
    sf_divn(num,den,rat2);
    
    sf_divn_combine(rat1,rat2,rat1);
}

/* 	$Id: Mwarpscan.c 744 2004-08-17 18:46:07Z fomels $	 */
