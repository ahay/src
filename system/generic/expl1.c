/* Anisotropic diffusion, 1-D, explicit cascade */
/*
  Copyright (C) 2019 University of Texas at Austin
  
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

static float *w, *y;
static int nr, n, nclip;

void expl1_init (int rect /* smoothing radius */,
		 int n1      /* data size */, 
		 float pclip /* percentage clip */)
/*< initialize >*/
{
    nr = rect;
    n = n1;

    w = sf_floatalloc(n);
    y = sf_floatalloc(n);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;
}

void expl1_close (void)
/*< free allocated storage >*/
{
    free(w);
    free(y);
}

void expl1_apply (float *x)
/*< apply diffusion >*/
{
    int ir, i;
    float a, r;

    /******** 
	      dx/dt = 1/w D' 1/w D u
	      x_{t+1} = (I + t 1/w D' 1/w D)^{-1} x_t
	      = (w + D' t/w D)^{-1} w x_t 
    *********/


    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = 1/(r*r);
	
	sf_grad2(n,x,w);

	if (1==ir) {
	    for (i=0; i < n; i++) {
		y[i] = w[i];
	    }
		
	    /* quantile destroys input */
	    a = sf_quantile(nclip,n,y);
	
	    if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
				__FILE__,nclip);
	}

	for (i=0; i < n; i++) {
	    w[i] = 1.0f/sqrtf(1.0f+w[i]/a);
	}

	y[0] = x[0]+ r*(x[1]-x[0])*w[0]*w[0];
	for (i=1; i < n-1; i++) {
	    y[i] = x[i]+ r*((x[i-1]-x[i])*w[i-1] + (x[i+1]-x[i])*w[i])*w[i];
	}
	y[n-1] = x[n-1]+ r*(x[n-2]-x[n-1])*w[n-2]*w[n-1];

        for (i=0; i < n; i++) {
	    x[i] = y[i];
	}
    }
}
