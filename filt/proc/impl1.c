/* Anisotropic diffusion, 1-D */
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
/*^*/

#include "impl1.h"

static float t, *w, *d1, *w1;
static int nstep, n, nclip;
static bool up;
static sf_tris slv;

void impl1_init (float r     /* radius */, 
		 int n1      /* data size */, 
		 float tau   /* duration */, 
		 float pclip /* percentage clip */, 
		 bool up_in  /* weight version */)
/*< initialize >*/
{
    float q;

    q = 0.25/SF_PI;
    t = r*r*q;

    nstep = t/tau;
    if (nstep > 1) {
	t /= nstep;
    } else {
	nstep = 1;
    }

    n = n1;

    w = sf_floatalloc(n);
    d1 = sf_floatalloc(n);
    w1 = sf_floatalloc(n);

    slv = sf_tridiagonal_init (n);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;

    sf_warning("%s: nstep=%d tau=%g nclip=%d",__FILE__,nstep,t,nclip);
}

void impl1_close (void)
/*< free allocated storage >*/
{
    free(w);
    free(d1);
    free(w1);

    sf_tridiagonal_close (slv);
}

void impl1_apply (float *x)
/*< apply diffusion >*/
{
    int istep, i;
    float a, xsum, wsum;

    /******** 
	      dx/dt = 1/w D' 1/w D u
	      x_{t+1} = (I + t 1/w D' 1/w D)^{-1} x_t
	      = (w + D' t/w D)^{-1} w x_t 
    *********/


    for (istep=0; istep < nstep; istep++) {
	sf_grad2(n,x,w);

	for (i=0; i < n; i++) {
	    w1[i] = w[i];
	}

	a = sf_quantile(nclip,n,w1);
	
	if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
			    __FILE__,nclip);

	wsum = xsum = 0.;
	for (i=0; i < n; i++) {
	    w[i] = sqrtf(1.+w[i]/a);
	    wsum += w[i];
	    xsum += x[i]*x[i];
	}

	sf_warning("step %d of %d, xsum=%g, wsum=%g, a=%g",  
		   istep, nstep, xsum, wsum, a);

	for (i=0; i < n; i++) {
	    if (up) w[i] = 1./w[i];
	    x[i] *= w[i];
	}

	for (i=0; i < n; i++) {
	    d1[i] = w[i];
	    if (up) {
		w1[i] = -t*d1[i];
	    } else {
		w1[i] = -t/d1[i];
	    }
	}

	d1[0] -= w1[0];
	for (i=1; i < n-1; i++) {
	    d1[i] -= w1[i] + w1[i-1];
	}
	d1[n-1] -= w1[n-2];

	sf_tridiagonal_define (slv, d1, w1); 
	sf_tridiagonal_solve (slv, x);
    }
}
