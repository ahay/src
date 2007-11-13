/* Digital freqlet transform */
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

#include "freqlet.h"

static int nt;
static bool inv, unit;
static float *w;
static sf_complex *t, *z;
static void (*transform)(bool);

static void haar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j;
    sf_complex z0;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    z0 = z[j];
	    for (i=0; i < nt-j; i += 2*j) {
		if (inv) {
		    t[i]   -= t[i+j]/(2*z0);
		    t[i+j] += t[i]*z0;
		} else {
		    t[i+j] += t[i]*z0/2;
		    t[i]   -= t[i+j]/z0;
		}
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    z0 = z[j];
	    for (i=0; i < nt-j; i += 2*j) {
		t[i+j] -= t[i]*z0;
		t[i]   += t[i+j]/(2*z0);
	    }	    
	}
    }
}


static void linear(bool adj) 
/* Lifting linear-interpolation transform in place */
{
    int i, j;
    sf_complex z0;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    z0 = z[j];
	    if (inv) {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   -= (t[i+j]/z0+t[i-j]*z0)/4;
		}
		t[0] -= t[j]/(2*z0);
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]*z0+t[i+2*j]/z0)/2;
		}	 
		if (i+j < nt) t[i+j] += t[i]*z0;
	    } else {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j] += t[i]*z0/4;
		    t[i-j] += t[i]/(4*z0);
		}
		t[j] += t[0]*z0/2;
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i]     -= t[i+j]/(2*z0);
		    t[i+2*j] -= t[i+j]*z0/2;
		}	 
		if (i+j < nt) t[i] -= t[i+j]/z0;
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    z0 = z[j];

	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] -= (t[i]*z0+t[i+2*j]/z0)/2;
		/* d = o - P e */
	    }	 
	    if (i+j < nt) t[i+j] -= t[i]*z0;    
	    t[0] += t[j]/(2*z0);
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]/z0+t[i-j]*z0)/4;
		/* s = e + U d */
	    }
	}
    }
}

void freqlet_init(int n /* data size */, bool inv1, bool unit1, char type) 
/*< allocate space >*/
{
    int i, j;
    float wi;

    inv = inv1;
    unit = unit1;

    for (nt=1; nt < n; nt *= 2) ;
    t = sf_complexalloc(nt);
    z = sf_complexalloc(nt);
    
    for (j=1; j <= nt/2; j *= 2) {
	z[j] = 1.; 
    }
    
    switch(type) {
	case 'h': 
	    transform = haar;
	    break;
	case 'l':
	    transform = linear;
	    break;
	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }

    if (unit) {
	w = sf_floatalloc(nt);

	w[0] = sqrtf((float) nt);
	wi = 0.5;	
	for (j=1; j <= nt/2; j *= 2, wi *= 2) {
	    
	    for (i=0; i < nt-j; i += 2*j) {
		w[i+j] = sqrtf(wi);
	    }
	}
    }
}

void freqlet_set(float w0)
/*< set frequency >*/
{
    int j;

    for (j=1; j <= nt/2; j *= 2) {
	z[j] = sf_cmplx(cosf(w0*j),sinf(w0*j));
    }
}

void freqlet_close(void) 
/*< deallocate space >*/
{
    free (t);
    free (z);
    if (unit) free(w);
}

void freqlet_lop(bool adj, bool add, int nx, int ny, 
		 sf_complex *x, sf_complex *y)
/*< linear operator >*/
{
    int it, i, j;

    sf_cadjnull (adj,add,nx,ny,x,y);

    if (adj) {
	t[0] = y[0];
	it = 1;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		if (it < ny) {
		    t[i+j]=y[it];
		    it++;
		} else {
		    t[i+j]=0.;
		}
	    }	    	    
	}

	if (unit) {
	    for (it=0; it < nt; it++) {
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
	    for (it=0; it < nt; it++) {
		t[it] *= w[it];
	    }
	}

	y[0] += t[0];
	it = 1;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		y[it] += t[i+j];
		it++;
		if (it >= ny) return;
	    }	    	    
	}
    }
}



