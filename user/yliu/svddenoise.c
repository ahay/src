/* SVD denoising */
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

#include "svd.h"
#include "svddenoise.h"

static float *a, *u, *v, *mid, *pxt;

void svddenoise_lop(int m, int n, float pclip, float *x, float *y)
/*< SVD denoise operator >*/
{
    int i, j, max, ka;
    int nclip;
    const double eps = 1.0e-5;

    a = sf_floatalloc(m*n);
    u = sf_floatalloc(m*m);
    v = sf_floatalloc(n*n);
    mid = sf_floatalloc(m*n);
    pxt = sf_floatalloc(m*n);

    max = m > n ? m : n;
    ka = max +1;
    nclip = (int)(m*pclip/100.+0.5);

    for (i=0;i<(m*n);i++) {
	a[i] = x[i];
    }

    svdinit(m,n,ka,eps);
    svduav(a,u,v);
    for(i=0;i<m;i++) {
	for(j=0;j<n;j++)  {
	    if((i>=nclip)&&(i==j))  {   		 
		a[i*n+j] = 0.;
	    }
	} 
    }
    brmul(u,a,m,m,n,mid);                     
    brmul(mid,v,m,n,n,pxt);   

    for (i=0;i<(m*n);i++) {
	y[i] = pxt[i];
    }

    free(a);
    free(u);
    free(v);
    free(mid);
    free(pxt);
}

/* 	$Id$	 */
