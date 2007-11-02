/* Digital wavelet transform */
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

#include "wavelet.h"

static int nt;
static bool inv;
static void (*transform)(bool);
static float *t, *tt, s;

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

static void ulinear(bool adj) 
/* Unitary linear-interpolation transform in place */
{
    int i, j;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		tt[i+j] = t[i+j]*s;
		tt[i] = t[i]*s;
		t[i+j] = 0.;
		t[i] = 0.;
	    }
	    t[0] += tt[0];
	    t[j] += tt[0];
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i-2*j] -= tt[i]/2;
		t[i-j] += tt[i];
		t[i] +=6.0f* tt[i]/2;
		t[i+j] += tt[i];
		t[i+2*j] -= tt[i]/2;
	    }
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+2*j] -= tt[i+j];
		t[i+j] += 2*tt[i+j];
		t[i] -= tt[i+j];
	    }	 
	    if (i+j < nt) {
		t[i+j] += tt[i+j];
		t[i] -= tt[i+j];
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-2*j; i += 2*j) {
		tt[i+j] = 2*t[i+j] - (t[i]+t[i+2*j]);
		/* d = o - P e */
	    }	 
	    if (i+j < nt) tt[i+j] = t[i+j]-t[i];    
	    tt[0] = t[0]+t[j];
	    for (i=2*j; i < nt-j; i += 2*j) {
		tt[i] = t[i-j]+t[i+j] + (6.0f*t[i]-t[i-2*j]-t[i+2*j])/2;
		/* s = e + U d */
	    }
	    for (i=0; i < nt-j; i += 2*j) {
		t[i+j] = tt[i+j]*s;
		t[i] = tt[i]*s;
	    }
	}
    }
}

static void haar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		if (inv) {
		    t[i]   -= t[i+j]/2;
		    t[i+j] += t[i];
		} else {
		    t[i+j] += t[i]/2;
		    t[i]   -= t[i+j];
		}
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		t[i+j] -= t[i];
		t[i]   += t[i+j]/2;
	    }	    
	}
    }
}

static void uhaar(bool adj) 
/* Unitary Haar transform */
{
    int i, j;
    float a, b;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		a = t[i]+t[i+j];
		b = t[i]-t[i+j];
		t[i+j] = a*s;
		t[i] = b*s;
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		a = t[i+j]-t[i];
		b = t[i+j]+t[i];
		t[i+j] = a*s;
		t[i] = b*s;
	    }	    
	}
    }
}

void wavelet_init(int n /* data size */, bool inv1, bool unit, char type) 
/*< allocate space >*/
{
    inv = inv1;

    for (nt=1; nt < n; nt *= 2) ;
    t = sf_floatalloc(nt);

    if (unit) s = sqrtf(0.5);

    switch(type) {
	case 'h': 
	    transform = unit? uhaar: haar;
	    break;
	case 'l':
	    transform = unit? ulinear: linear;
	    break;
	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }

    if (transform == ulinear) tt = sf_floatalloc(nt);
}

void wavelet_close(void) 
/*< deallocate space >*/
{
    free (t);
    if (transform == ulinear) free(tt);
}

void wavelet_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int it, i, j;

    sf_adjnull (adj,add,nx,ny,x,y);

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
