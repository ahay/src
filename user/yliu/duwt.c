/* Digital undecimated (stationary) wavelet transform by lifting scheme */
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
#include <rsf.h>
#include <math.h>
#include "duwt.h"

static int nt, scale;
static bool inv, unit;
static void (*transform)(bool);
static float *t, *w;

static void biorthogonal(bool adj) 
/* Lifting CDF 9/7 biorthogonal wavelet transform in place */
{
    int i,j;
    float a;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    if (inv) {                            /*reverse dwt9/7 transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]  /= a;
	        }
	        t[0] /= a;                         /*left boundary*/
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] *= a;
		    /* Undo Scale */
	        }
	        if (i+j < nt) t[i+j] *= a;         /*right boundary*/  

                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   += (t[i+j]+t[i-j])*a;
		    /* Undo Update 2 */
	        }
	        t[0] += 2*a*t[j];                  /*left boundary*/

	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])*a;
		    /* Undo Predict 2 */
	        }	 
	        if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
                    /* Undo Step 2 */

                a= 0.05298011854f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   += (t[i+j]+t[i-j])*a;
		    /* Undo Update 1 */
	        }
	        t[0] += 2*a*t[j];                  /*left boundary*/ 

	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])*a;
		    /* Undo Predict 1 */
	        }	 
	        if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
                    /* Undo Step 1 */


	    } else {                               /*adjoint transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]  *= a;
	        }
	        t[0] *= a;                         /*left boundary*/
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] /= a;
		    /* Undo Scale */
	        }
	        if (i+j < nt) t[i+j] /= a;         /*right boundary*/  

                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j]   -= t[i]*a;
                    t[i-j]   -= t[i]*a;
		    /* Undo Update 2 */
	        }
	        t[j] -= 2*a*t[0];                  /*left boundary*/

	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
                    t[i]     -= t[i+j]*a;
                    t[i+2*j] -= t[i+j]*a;
		    /* Undo Predict 2 */
	        }	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j];  /*right boundary*/  
                    /* Undo Step 2 */

                a= 0.05298011854f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j]   -= t[i]*a;
                    t[i-j]   -= t[i]*a;
		    /* Undo Update 1 */
	        }
	        t[j] -= 2*a*t[0];                  /*left boundary*/ 

	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
                    t[i]     -= t[i+j]*a;
                    t[i+2*j] -= t[i+j]*a;
		    /* Undo Predict 1 */
	        }	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j];  /*right boundary*/  
                    /* Undo Step 1 */

	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {        /*different scale*/
	    a = -1.586134342f;
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 1 */
	    }	 
	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= -0.05298011854f;
	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 1 */
	    }
                /* Step 1 */

	    a = 0.8829110762f;
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 2 */
	    }	 
	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= 0.4435068522f;
	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 2 */
	    }
                /* Step 2 */

            a= 1/(1.230174105f);
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] *= a;
	    }	 
	    if (i+j < nt) t[i+j] *= a;         /*right boundary*/  
	    t[0] /= a;                         /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]  /= a;
		/* Scale */
	    }
	}
    }

}

static void linear(bool adj) 
/* Lifting linear-interpolation transform in place */
{
    int i, j;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    if (inv) {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   -= (t[i+j]+t[i-j])/4;
		}
		t[0] -= t[j]/2;
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])/2;
		}	 
		if (i+j < nt) t[i+j] += t[i];
	    } else {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j] += t[i]/4;
		    t[i-j] += t[i]/4;
		}
		t[j] += t[0]/2;
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i]     -= t[i+j]/2;
		    t[i+2*j] -= t[i+j]/2;
		}	 
		if (i+j < nt) t[i] -= t[i+j];
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] -= (t[i]+t[i+2*j])/2;
		/* d = o - P e */
	    }	 
	    if (i+j < nt) t[i+j] -= t[i];    
	    t[0] += t[j]/2;
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])/4;
		/* s = e + U d */
	    }
	}
    }
}

static void haar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j, shift, k;
    float *temp;
    temp = sf_floatalloc(nt);


    if (adj) {
	shift = 0;
	for (j=(scale-1); j>= 1; j--) {
	    shift = pow(2,(j-1));

	    for (i=0; i< nt; i++) {
		temp[i]=t[j*nt+i];
	    }
	    if (inv) {
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] - t[(j-1)*nt+k]/2;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += t[j*nt+k];
		}

	    } else {
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] +=  t[j*nt+k]/2;
		}

		for (i=0; i< nt; i++) {
		    if ((i-shift)<0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] - t[(j-1)*nt+k];
		}
	    }
	}

    } else {
	shift = 0;
	for (j=1; j< scale; j++) {
	    shift = pow(2,(j-1));
	    for (i=0; i< nt; i++) {
		t[j*nt+i] = t[(j-1)*nt+i];
	    }
	    for (i=0; i< nt; i++) {
		if ((i+shift) > (nt-1)) {
		    k = i + shift - nt;
		} else {
		    k = i + shift;
		}
		t[(j-1)*nt+i] -= t[j*nt+k];
	    }
	    for (i=0; i< nt; i++) {
		if ((i-shift) < 0) {
		    k = i - shift + nt;
		} else {
		    k = i - shift;
		}
		temp[k] = t[j*nt+i] + t[(j-1)*nt+k]/2;
	    }
	    for (i=0; i< nt; i++) {
		t[j*nt+i]=temp[i];
	    }
	}
    }
}

void wavelet_init(int n /* data size */, bool inv1, bool unit1, char type, int max1) 
/*< allocate space >*/
{
    int i, j;
    float wi;

    inv = inv1;
    unit = unit1;
    scale = max1;

    for (nt=1; nt < n; nt *= 2) ;
    t = sf_floatalloc(scale*n);

    switch(type) {
	case 'h': 
	    transform = haar;
	    break;
	case 'l':
	    transform = linear;
	    sf_error("Type=%c is unavailable now!",type);
	    break;
	case 'b':
	    transform = biorthogonal;
	    sf_error("Type=%c is unavailable now!",type);
	    break;

	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }

    if (unit) {
	w = sf_floatalloc(scale*nt);
	for (j=0; j< nt; j++) {
	    w[0*nt+j] = sqrtf(0.5);
	}
	wi = 1.;
	for (j=1; j< scale; j++, wi *=2) {
	    for (i=0; i< nt; i++) {
		w[j*nt+i] = sqrtf(wi);
	    }
	}
    }
}

void wavelet_close(void) 
/*< deallocate space >*/
{
    free (t);
    if (unit) free(w);
}

void wavelet_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int it, i, j;

    sf_adjnull (adj,add,nx,ny,x,y);

    if (adj) {
	for (j=0; j< scale; j++) {
	    for(i=0; i< nx; i++) {
		t[j*nx+i] = y[j*nx+i];
	    }
	}
	if (unit) {
	    for (it=0; it < scale*nt; it++) {
		if (inv) {
		    t[it] /= w[it];
		} else {
		    t[it] *= w[it];
		}
	    }
	}


    } else {
	for (it=0; it < nx; it++) {
	    t[it]=x[it];
	}
	for (it=nx; it < nt; it++) {
	    t[it] = 0.;
	}
    }

    transform(adj);    

    if (adj) {
	for (it=0; it < nx; it++) {
	    x[it] += t[it];
	}

    } else {
	if (unit) {
	    for (it=0; it < scale*nt; it++) {
		t[it] *= w[it];
	    }
	}
	for (j=0; j< scale; j++) {
	    for(i=0; i< nx; i++) {
		y[j*nx+i] += t[j*nx+i];
	    }
	}

    }
}

/* 	$Id$	 */
