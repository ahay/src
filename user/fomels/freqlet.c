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
#ifdef SF_HAS_COMPLEX_H
		    t[i]   -= t[i+j]/(2*z0);
		    t[i+j] += t[i]*z0;
#else
		    t[i] = sf_cadd(t[i],
				   sf_crmul(
				       sf_cmul(t[i+j],sf_conjf(z0)),-0.5));
		    t[i+j] = sf_cadd(t[i+j],sf_cmul(t[i],z0));		    
#endif
		} else {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += t[i]*z0/2;
		    t[i]   -= t[i+j]/z0;
#else
		    t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_cmul(t[i],z0),0.5));
		    t[i]   = sf_cadd(t[i],sf_cmul(t[i+j],sf_neg(sf_conjf(z0))));
#endif
		}
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    z0 = z[j];
	    for (i=0; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i+j] -= t[i]*z0;
		t[i]   += t[i+j]/(2*z0);
#else
		t[i+j] = sf_cadd(t[i+j],sf_cmul(t[i],sf_neg(z0)));
		t[i]   = sf_cadd(t[i],
				 sf_crmul(sf_cmul(t[i+j],sf_conjf(z0)),0.5));
#endif
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
#ifdef SF_HAS_COMPLEX_H
		    t[i]   -= (t[i+j]/z0+t[i-j]*z0)/4;
#else
		    t[i] = sf_cadd(t[i],
				   sf_crmul(
				       sf_cadd(
					   sf_cmul(t[i+j],sf_conjf(z0)),
					   sf_cmul(t[i-j],z0)),
				       -0.25))
#endif
		}
#ifdef SF_HAS_COMPLEX_H
		t[0] -= t[j]/(2*z0);
#else
		t[0] = sf_cadd(t[0],sf_crmul(sf_cmul(t[j],sf_conjf(z0)),-0.5));
#endif
		for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t[i]*z0+t[i+2*j]/z0)/2;
#else
		    t[i+j] = sf_cadd(t[i+j],
				     sf_crmul(
					 sf_cadd(sf_cmul(t[i],z0),
						 sf_cmul(t[i+2*j],sf_conjf(z0))),
					 0.5));
#endif
		}
#ifdef SF_HAS_COMPLEX_H		
		if (i+j < nt) t[i+j] += t[i]*z0;
#else
		if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_cmul(t[i],z0));
#endif
	    } else {
		for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i+j] += t[i]*z0/4;
		    t[i-j] += t[i]/(4*z0);
#else
		    t[i+j] = sf_cadd(t[i+j],
				    sf_crmul(sf_mul(t[i],z0),0.25));
		    t[i-j] = sf_cadd(t[i-j],
				    sf_crmul(sf_mul(t[i],sf_conjf(z0)),0.25));
#endif
		}
#ifdef SF_HAS_COMPLEX_H	
		t[j] += t[0]*z0/2;
#else
		t[j] = sf_cadd(t[j],sf_crmul(sf_cmul(t[0],z0),0.5));
#endif
		for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H	
		    t[i]     -= t[i+j]/(2*z0);
		    t[i+2*j] -= t[i+j]*z0/2;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cmul(t[i+j],sf_conjf(z0)),
						 -0.5));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(sf_cmul(t[i+j],z0),
							 -0.5));
#endif
		}	 
#ifdef SF_HAS_COMPLEX_H	
		if (i+j < nt) t[i] -= t[i+j]/z0;
#else
		if (i+j < nt) t[i] = sf_cadd(t[i],
					     sf_neg(
						 sf_cmul(t[i+j],sf_conjf(z0))));
#endif
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {
	    z0 = z[j];

	    for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H	
		t[i+j] -= (t[i]*z0+t[i+2*j]/z0)/2;
		/* d = o - P e */
#else
		t[i+j] = sf_cadd(t[i+j],
				 sf_crmul(
				     sf_cadd(sf_cmul(t[i],z0),
					     sf_cmul(t[i+2*j],sf_conjf(z0))),
				     -0.5));
#endif
	    }	 
#ifdef SF_HAS_COMPLEX_H
	    if (i+j < nt) t[i+j] -= t[i]*z0;    
	    t[0] += t[j]/(2*z0);
#else
	    if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_cneg(sf_cmul(t[i],z0)));
	    t[0] = sf_cadd(t[0],sf_crmul(sf_cmul(t[j],sf_conjf(z0)),0.5));
#endif
	    for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t[i+j]/z0+t[i-j]*z0)/4;
		/* s = e + U d */
#else
		t[i] = sf_cadd(t[i],sf_crmul(sf_cadd(
						 sf_cmul(t[i+j],sf_conjf(z0)),
						 sf_cmul(t[i-j],z0)),0.25));
#endif
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
	z[j] = sf_cmplx(1.,0.); 
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
		    t[i+j]=sf_cmplx(0.,0.);
		}
	    }	    	    
	}

	if (unit) {
	    for (it=0; it < nt; it++) {
		if (inv) {
#ifdef SF_HAS_COMPLEX_H
		    t[it] /= w[it];
#else
		    t[it] = sf_crmul(t[it],1.0f/w[it]);
#endif
		} else {
#ifdef SF_HAS_COMPLEX_H
		    t[it] *= w[it];
#else
		    t[it] = sf_crmul(t[it],w[it]);
#endif
		}
	    }
	}
    } else {
	for (it=0; it < nx; it++) {
	    t[it]=x[it];
	}
	for (it=nx; it < nt; it++) {
	    t[it] = sf_cmplx(0.,0.);
	}
    }

    transform(adj);    

    if (adj) {
	for (it=0; it < nx; it++) {
#ifdef SF_HAS_COMPLEX_H
	    x[it] += t[it];
#else
	    x[it] = sf_add(x[it],t[it]);
#endif
	}
    } else {
	if (unit) {
	    for (it=0; it < nt; it++) {
#ifdef SF_HAS_COMPLEX_H
		t[it] *= w[it];
#else
		t[it] = sf_crmul(t[it],w[it]);
#endif
	    }
	}

#ifdef SF_HAS_COMPLEX_H
	y[0] += t[0];
#else
	y[0] = sf_cadd(y[0],t[0]);
#endif
	it = 1;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		y[it] += t[i+j];
#else
		y[it] = sf_cadd(y[it],t[i+j]);
#endif
		it++;
		if (it >= ny) return;
	    }	    	    
	}
    }
}



