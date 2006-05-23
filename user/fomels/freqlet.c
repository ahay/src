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

static int nt, nk;
static bool inv;
static float *t, **a;

static void linear(bool adj) 
/* Lifting linear-interpolation transform in place */
{
    int i, j, k;

    if (adj) {
	for (j=nt/2, k=nk-1; j >= 1; j /= 2, k--) {
	    if (inv) {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   -= (t[i+j]+t[i-j])/(4*a[k][i]);
		}
		t[0] -= t[j]/(2*a[k][0]);
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])/(2*a[k][i+j]);
		}	 
		if (i+j < nt) t[i+j] += t[i]/a[k][i+j];
	    } else {
		for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j] += t[i]/(4*a[k][i]);
		    t[i-j] += t[i]/(4*a[k][i]);
		}
		t[j] += t[0]/(2*a[k][0]);
		for (i=0; i < nt-2*j; i += 2*j) {
		    t[i]     -= t[i+j]/(2*a[k][i+j]);
		    t[i+2*j] -= t[i+j]/(2*a[k][i+j]);
		}	 
		if (i+j < nt) t[i] -= t[i+j]/a[k][i+j];
	    }
	}
    } else {
	for (j=1, k=0; j <= nt/2; j *= 2, k++) {
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] -= (t[i]+t[i+2*j])/(2*a[k][i+j]);
		/* d = o - P e */
	    }	 
	    if (i+j < nt) t[i+j] -= t[i]/a[k][i+j];    
	    t[0] += t[j]/(2*a[k][0]);
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])/(4*a[k][i]);
		/* s = e + U d */
	    }
	}
    }
}

void freqlet_init(int n, bool inv1) 
/*< allocate space >*/
{
    int j;

    inv = inv1;

    for (nt=1; nt < n; nt *= 2) ;
    for (j=1, nk=0; j <= nt/2; j *= 2, nk++);

    t = sf_floatalloc(nt);
    a = sf_floatalloc2(nk,nt);
}

void freqlet_set(const float* aa)
/*< set local frequency >*/
{
    int i,k;
    float c;

    for (i=0; i < nt; i++) {
	c = aa[i];
	a[0][i] = c;
	for (k=1; k < nk; k++) {
	    c = 2.*c*c-1.;
	    a[k][i] = c;
	}
    }
}

void freqlet_close(void) 
/*< deallocate space >*/
{
    free (t);
    free (*a);
    free (a);
}

void freqlet_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
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

    linear(adj);

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
