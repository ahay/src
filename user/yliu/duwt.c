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
/* Lifting CDF 9/7 biorthogonal undecimated wavelet transform in place */
{
    float a;
    int i, j, shift, k, shiftt;
    float *temp, *tempt;
    temp = sf_floatalloc(nt);
    tempt = sf_floatalloc(nt);

    if (adj) {
	shift = 0;
	shiftt = 1;
	for (j=(scale-1); j>= 1; j--) {
	    shift = pow(2,(j-1));
	    shiftt = pow(2,j);
	    if (inv) {                            /*reverse UWT 9/7 transform*/
		a= 1.230174105f;
		for (i=0; i< nt; i++) {
		    t[(j-1)*nt+i] *=a;
		    t[j*nt+i] /=a;
		}
		/* Undo Scale */
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    temp[i] = t[j*nt+k];
		}
		a= -0.4435068522f;		
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] + tempt[k];
		}
		a = -0.8829110762f;
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += tempt[k];
		}
		/* Undo Step2 */

		for (i=0; i< nt; i++) {
		    temp[i] = t[j*nt+i];
		}
		a= 0.05298011854f;		
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[i] + tempt[k];
		}
		a = 1.586134342f;
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += tempt[k];
		}
		/* Undo Step1 */

	    } else {                               /*adjoint UWT 9/7 transform*/
		a= 1.230174105f;
		for (i=0; i< nt; i++) {
		    t[(j-1)*nt+i] /=a;
		    t[j*nt+i] *=a;
		}
		/* Adjoint Scale */
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    temp[i] = t[j*nt+k];
		}

		a= -0.4435068522f;	
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += tempt[k];
		}

		a = -0.8829110762f;
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] + tempt[k];
		}
		/* Adjoint Step2 */

		for (i=0; i< nt; i++) {
		    temp[i] = t[j*nt+i];
		}

		a= 0.05298011854f;		
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += tempt[k];
		}

		a = 1.586134342f;
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[i] + tempt[k];
		}
		/* Adjoint Step1 */
	    }
	}
    } else {
	shift = 0;
	shiftt = 1;
	for (j=1; j< scale; j++) {
	    shift = pow(2,(j-1));
	    shiftt = pow(2,j);
	    for (i=0; i< nt; i++) {
		t[j*nt+i] = t[(j-1)*nt+i];
	    }
	    a = -1.586134342f;
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
	    }
	    for (i=0; i< nt; i++) {
		if ((i+shift) > (nt-1)) {
		    k = i + shift - nt;
		} else {
		    k = i + shift;
		}
		t[(j-1)*nt+i] += tempt[k];
	    }
	    a = -0.05298011854f;
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
	    }
	    for (i=0; i< nt; i++) {
		if ((i-shift) < 0) {
		    k = i - shift + nt;
		} else {
		    k = i - shift;
		}
		t[j*nt+i] = t[j*nt+i] + tempt[k];
	    }
              /* Step 1 */

	    a = 0.8829110762f;
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[j*nt+i] + t[j*nt+k])*a;
	    }
	    for (i=0; i< nt; i++) {
		if ((i+shift) > (nt-1)) {
		    k = i + shift - nt;
		} else {
		    k = i + shift;
		}
		t[(j-1)*nt+i] += tempt[k];
	    }
            a= 0.4435068522f;
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])*a;
	    }
	    for (i=0; i< nt; i++) {
		if ((i-shift) < 0) {
		    k = i - shift + nt;
		} else {
		    k = i - shift;
		}
		temp[k] = t[j*nt+i] + tempt[k];
	    }
	    for (i=0; i< nt; i++) {
		if ((i-shift) < 0) {
		    k = i - shift + nt;
		} else {
		    k = i - shift;
		}
		t[j*nt+k]=temp[i];
	    }
             /* Step 2 */

            a= 1/(1.230174105f);
	    for (i=0; i< nt; i++) {
		t[(j-1)*nt+i] *=a;
		t[j*nt+i] /=a;
	    }
		/* Scale */
	}
    }
}

static void linear(bool adj) 
/* Lifting linear-interpolation transform in place */
{
    int i, j, shift, k, shiftt;
    float *temp, *tempt;
    temp = sf_floatalloc(nt);
    tempt = sf_floatalloc(nt);

    if (adj) {
	shift = 0;
	shiftt = 1;
	for (j=(scale-1); j>= 1; j--) {
	    shift = pow(2,(j-1));
	    shiftt = pow(2,j);

	    for (i=0; i< nt; i++) {
		temp[i]=t[j*nt+i];
	    }
	    if (inv) {

		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])/2;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift) < 0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] - tempt[k]/2;
		}
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])/2;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] = tempt[k] + t[(j-1)*nt+i];
		}
	    } else {
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[j*nt+i] + t[j*nt+k])/2;
		}
		for (i=0; i< nt; i++) {
		    if ((i+shift) > (nt-1)) {
			k = i + shift - nt;
		    } else {
			k = i + shift;
		    }
		    t[(j-1)*nt+i] += tempt[k]/2;
		}
		for (i=0; i<nt; i++) {
		    if ((i-shiftt) < 0) {
			k = i -shiftt + nt;
		    } else {
			k = i - shiftt;
		    }
		    tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])/2;
		}
		for (i=0; i< nt; i++) {
		    if ((i-shift)<0) {
			k = i - shift + nt;
		    } else {
			k = i - shift;
		    }
		    t[j*nt+i] = temp[k] - tempt[k];
		}
	    }
	}

    } else {
	shift = 0;
	shiftt = 1;
	for (j=1; j< scale; j++) {
	    shift = pow(2,(j-1));
	    shiftt = pow(2,j);
	    for (i=0; i< nt; i++) {
		t[j*nt+i] = t[(j-1)*nt+i];
	    }
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[j*nt+i] + t[j*nt+k])/2;
	    }
	    for (i=0; i< nt; i++) {
		if ((i+shift) > (nt-1)) {
		    k = i + shift - nt;
		} else {
		    k = i + shift;
		}
		t[(j-1)*nt+i] -= tempt[k];
	    }
	    for (i=0; i<nt; i++) {
		if ((i-shiftt) < 0) {
		    k = i -shiftt + nt;
		} else {
		    k = i - shiftt;
		}
		tempt[i] = (t[(j-1)*nt+i] + t[(j-1)*nt+k])/2;
	    }
	    for (i=0; i< nt; i++) {
		if ((i-shift) < 0) {
		    k = i - shift + nt;
		} else {
		    k = i - shift;
		}
		temp[k] = t[j*nt+i] + tempt[k]/2;
	    }
	    for (i=0; i< nt; i++) {
		t[j*nt+i]=temp[i];
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
	    break;
	case 'b':
	    transform = biorthogonal;
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
