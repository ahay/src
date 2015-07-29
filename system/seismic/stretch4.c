/* Inverse spline interpolation */
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

#include "stretch4.h"

#ifndef _stretch4_h

typedef struct Map4 *map4;
/* abstract data type */
/*^*/

#endif

struct Map4 {
    int nt, nd, ib, ie;
    float t0,dt, eps;
    int *x; 
    bool *m;
    float **w, *diag, *offd[3];
    sf_bands slv;
    sf_tris tslv;
};

map4 stretch4_init (int n1, float o1, float d1 /* regular axis */, 
		    int nd                     /* data samples */, 
		    float eps                  /* regularization */)
/*< initialize >*/
{
    int i;
    map4 str;
    
    str = (map4) sf_alloc (1, sizeof(*str));

    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->eps = eps;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc2 (4,nd);
    str->diag = sf_floatalloc (str->nt);
    
    for (i = 0; i < 3; i++) {
	str->offd[i] = sf_floatalloc (str->nt-1-i);
    }
  
    str->slv = sf_banded_init (str->nt,3);
    str->tslv = sf_spline4_init(str->nt);

    return str;
}

void stretch4_define (map4 str, const float* coord /* [nd] */)
/*< set coordinates >*/
{
    int id, ix, i1, n1, i, j, i2;
    float rx, d, o[3], *w;
    
    n1 = str->nt;

    d = str->eps*2./3.;
    o[0] = -str->eps/8.;
    o[1] = -str->eps/5.;
    o[2] = -str->eps/120.;

    for (i1 = 0; i1 < n1; i1++) {
	/* regularization */
	str->diag[i1] = d;
	for (j=0; j < 3; j++) {
	    if (i1 < n1-1-j) str->offd[j][i1] = o[j];
	}
    }
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt - 1.; 
	ix = floorf(rx);
	rx -= ix;

	if (ix <= -4 || ix >= n1) {
	    str->m[id] = true; 
	    continue;
	}

	str->x[id] = ix; 
	str->m[id] = false; 
	w = str->w[id];

	sf_spline4_int(rx,w);
	
	i1 = SF_MAX(0,-ix);
	i2 = SF_MIN(4,n1-ix);

	for (i = i1; i < i2; i++) { 
	    str->diag[ix+i] += w[i] * w[i];
	    for (j=0; j < i2-i-1; j++) {
		str->offd[j][ix+i] += w[i] * w[i+j+1];
	    }
	}
    }

    sf_banded_define (str->slv, str->diag, str->offd);
    
    str->ib = -1;
    for (i1 = 0; i1 < n1; i1++) {
	if (str->diag[i1] != d) {
	    str->ib = i1-1; 
	    break;
	}
    }
    
    str->ie = n1+3;
    for (i1 = n1-1; i1 >= 0; i1--) {
	if (str->diag[i1] != d) {
	    str->ie = i1+4;
	    break;
	}
    }
}

void cstretch4_apply (map4 str, 
		      const sf_complex* ord /* [nd] */, 
		      sf_complex* mod       /* [n1] */)
/*< complex transform ordinates to model >*/
{    
    int id, it;
    float *real, *imag, *ford;
    
    real = sf_floatalloc(str->nt);
    imag = sf_floatalloc(str->nt);
    ford = sf_floatalloc(str->nd);

    for (id=0; id < str->nd; id++) {
	ford[id] = crealf(ord[id]);
    }

    stretch4_apply (str,ford,real);

    for (id=0; id < str->nd; id++) {
	ford[id] = cimagf(ord[id]);
    }

    stretch4_apply (str,ford,imag);

    for (it=0; it < str->nt; it++) {
	mod[it] = sf_cmplx(real[it],imag[it]);
    }

    free(real);
    free(imag);
    free(ford);
}


void stretch4_apply (map4 str, 
		     const float* ord /* [nd] */, 
		     float* mod       /* [n1] */)
/*< transform ordinates to model >*/
{
    int id, it, i, nt, i1, i2;
    float *w, *mm;
    
    mm = str->diag;
    nt = str->nt;

    for (it = 0; it < nt; it++) {
	mm[it] = 0.;
    }
    
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	it = str->x[id]; 
	w = str->w[id]; 
	
	i1 = SF_MAX(0,-it);
	i2 = SF_MIN(4,nt-it);

	for (i=i1; i < i2; i++) {
	    mm[it+i] += w[i]*ord[id];
	}
    }    

    sf_banded_solve (str->slv, mm);

    for (it = 0; it <= str->ib; it++) {
	mm[it] = 0.;
    }
    
    for (it = str->ie; it < nt; it++) {
	mm[it] = 0.;
    }

    sf_spline4_post(nt,0,nt,mm,mod);

    for (it = 0; it <= str->ib; it++) {
	mod[it] = 0.;
    }
    
    for (it = str->ie; it < nt; it++) {
	mod[it] = 0.;
    }
}

void cstretch4_invert (map4 str, 
		       sf_complex* ord       /* [nd] */, 
		       const sf_complex* mod /* [n1] */)
/*< convert model to ordinates by spline interpolation >*/
{    
    int id, it;
    float *real, *imag, *fmod;
    
    real = sf_floatalloc(str->nd);
    imag = sf_floatalloc(str->nd);
    fmod = sf_floatalloc(str->nt);

    for (it=0; it < str->nt; it++) {
	fmod[it] = crealf(mod[it]);
    }

    stretch4_invert (str,real,fmod);

    for (it=0; it < str->nt; it++) {
	fmod[it] = cimagf(mod[it]);
    }

    stretch4_invert (str,imag,fmod);

    for (id=0; id < str->nd; id++) {
	ord[id] = sf_cmplx(real[id],imag[id]);
    }

    free(real);
    free(imag);
    free(fmod);
}

void stretch4_invert (map4 str, 
		      float* ord       /* [nd] */, 
		      const float* mod /* [n1] */)
/*< convert model to ordinates by spline interpolation >*/
{
    int id, it, i, nt, i1, i2;
    float *w, *mm;

    mm = str->diag;
    nt = str->nt;

    for (it = 0; it < nt; it++) {
	mm[it] = mod[it];
    }

    sf_tridiagonal_solve(str->tslv, mm);

    for (id = 0; id < str->nd; id++) {
	ord[id] = 0.;
	if (str->m[id]) continue;
	
	it = str->x[id]; 
	w = str->w[id]; 
	
	i1 = SF_MAX(0,-it);
	i2 = SF_MIN(4,nt-it);

	for (i=i1; i < i2; i++) {
	    ord[id] += w[i]*mm[it+i];
	}
    } 
}

void stretch4_close (map4 str)
/*< free allocated storage >*/
{
    int i;

    free (str->x);
    free (str->m);
    free (str->w[0]);
    free (str->w);
    free (str->diag);

    for (i = 0; i < 3; i++) {
	free (str->offd[i]);
    }
    
    sf_banded_close (str->slv);
    free (str);
}

/* 	$Id: stretch4.c 11336 2013-11-13 21:21:13Z sfomel $	 */
