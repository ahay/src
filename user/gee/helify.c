/* Helix derivative by Kolmogoroff factorization */
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

#include "kolmog.h"

#define N12 4096

float helify(float a, float b  /* filter coefficients */, 
	     int k2            /* filter breakpoint */,
	     int nf            /* filter size */, 
	     float* f          /* filter */)
/*< Factor an isotropic Laplacian on a helix >*/
{
    const int n1 = 64, n2 = 64;
    const float gamma = 0.631974;
    float cy[N12], middle, z1, z2, scale;
    int i1;

    middle = 4.*(gamma + (1.-gamma)*0.5);
    for (i1=0; i1 < N12; i1++) {
	z1 = cosf (2*SF_PI*i1/(n1*n2));
	z2 = cosf (2*SF_PI*i1/(   n2));
	cy[i1] = a + b*(middle -2*gamma*(z1+z2) - 2*(1-gamma)*z1*z2); 
    }

    kolmog_init(N12,0,0);
    kolmog2(cy);

    scale = 1./cy[0];
    for (i1=0; i1 < k2; i1++) {
	f[i1]  = cy[i1+1]*scale;
    } 
    for (i1=k2; i1 < nf; i1++) {
	f[i1] = cy[n1+2-nf+i1]*scale;
    }

    return (scale*scale);
}
