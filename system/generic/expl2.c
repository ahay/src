/* Anisotropic diffusion, 2-D, by box cascade */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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

static float **y, **w;
static int n1, n2, n, nclip, nr;

void expl2_init (int rect   /* radius */, 
		 int n1_in, int n2_in /* data size */, 
		 float pclip          /* percentage clip */)
/*< initialize >*/
{
    int i;

    nr = rect;
    n1 = n1_in;
    n2 = n2_in;
    n = n1*n2;
    
    y = sf_floatalloc2(n1,n2);
    w = sf_floatalloc2(n1,n2);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;

    for (i=0; i < n; i++) {
	w[0][i] = 1.;
    }
}

void expl2_close (void)
/*< free allocated storage >*/
{
    free(*y);
    free(y);
    free(*w);
    free(w);
}

void expl2_apply (float **x)
/*< apply diffusion >*/
{
    int i1, i2, ir, i;
    float a, r, t;

    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = 0.5/(r*r);
	
	sf_sobel2(n1,n2,x,w);

	if (1==ir) {
	    for (i=0; i < n; i++) {
		y[0][i] = w[0][i];
	    }

	    a = sf_quantile(nclip,n,y[0]);
	    if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
				__FILE__,nclip);
	}

	for (i=0; i < n; i++) {
	    w[0][i] = 1.0f/sqrtf(1.0f+w[0][i]/a); 	
	}
	
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {

		t = 0.0f;
	       	if (i1 < n1-1) {
		    t += (x[i2][i1+1]-x[i2][i1])*w[i2][i1];
		}
		if (i1 > 0) {
		    t += (x[i2][i1-1]-x[i2][i1])*w[i2][i1-1];
		}
		if (i2 < n2-1) {
		    t += (x[i2+1][i1]-x[i2][i1])*w[i2][i1];
		}
		if (i2 > 0) {
		    t += (x[i2-1][i1]-x[i2][i1])*w[i2-1][i1];
		}
		
		y[i2][i1] = x[i2][i1] + r*w[i2][i1]*t;
	    }
	}

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		x[i2][i1] = y[i2][i1];
	    }
	}
    }
}

