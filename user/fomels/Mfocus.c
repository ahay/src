/* Focusing indicator. 

Takes: rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
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
#include <math.h>

#include <rsf.h>

#include "divn.h"
#include "hilbert.h"

int main (int argc, char* argv[])
{
  int i3, n3, dim1, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *num, *den, *rat1, *rat2, *dat, mean;
    char key[6];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    if (!sf_getint("dim",&dim1)) {
      /* dimensionality */
      dim1=dim;
      n3=1;
    } else {
      n3 = sf_leftsize(in,dim1);
    }

    n12 = 1;
    for (i=0; i < dim1; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	n12 *= n[i];
    }

    dat = sf_floatalloc(n12);
    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    rat1 = sf_floatalloc(n12);
    rat2 = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    divn_init(dim1, n12, n, rect, niter);

    for (i3=0; i3 < n3; i3++) {
      sf_floatread(dat,n12,in);

      mean=0.;
      for (i=0; i < n12; i++) {
	dat[i] *= dat[i];
	mean += dat[i]*dat[i];
      }
      mean = sqrtf(n12/mean);

      for (i=0; i < n12; i++) {
	num[i] = mean;
	den[i] = dat[i]*mean;
      }
      
      divn (num, den, rat1);
    
      for (i=0; i < n12; i++) {
	den[i] = 1.;
	num[i] = dat[i];
      }
    
      divn (num, den, rat2);

      for (i=0; i < n12; i++) {
	rat1[i] *= rat2[i];
      }
    
      /*
	
      for (i=0; i < n12; i++) {
      if (rat1[i] > 0.) {
      if (rat2[i] > 0. || -rat2[i] < rat1[i]) {
      rat1[i] = fabsf(rat1[i]*rat2[i]);
      } else {
      rat1[i] = -fabsf(rat1[i]*rat2[i]);
      }
      } else {
      if (rat2[i] < 0. || rat2[i] < -rat1[i]) {
      rat1[i] = -fabsf(rat1[i]*rat2[i]);
      } else {
      rat1[i] = fabsf(rat1[i]*rat2[i]);
      }
      }
      }
      
      */

      sf_floatwrite(rat1,n12,out);
    }

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
