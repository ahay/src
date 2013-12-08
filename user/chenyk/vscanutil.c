/* Velocity analysis unitlity */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include<rsf.h>
#include"vscanutil.h"

#ifndef _vscanutil_h

typedef struct Vint1 *vint1;
typedef struct Fint1 *fint1;
/* abstract data types */
/*^*/

typedef float (*mapfunc)(float,int);
/*^*/

#endif

struct Vint1 {
    float **spl, w[4];
    int n1, nw, dim;
    sf_tris slv;
};

struct Fint1 {
    float *spl, *t, w[4];
    int n1, nw, nt, ir;
    sf_tris slv;
};

float* fint1_coeff (fint1 fnt, int n)
/*< extract n-th spline coefficient >*/
{
    float* coeff;
    coeff = &(fnt->spl[n+fnt->nw]);

    return coeff;
}

float* vint1_coeff (vint1 fnt, int n, int dim)
/*< extract n-th spline coefficient for dimension dim >*/
{
    float* coeff;
    coeff = &(fnt->spl[dim][n+fnt->nw]);
    
    return coeff;
}

vint1 vint1_init (int nw  /* trace extension */, 
		  int n1  /* trace length */, 
		  int dim /* number of functions */)
/*< initialize multi-function interpolation >*/
{
    vint1 fnt;

    fnt = (vint1) sf_alloc (1, sizeof (*fnt));

    fnt->nw = nw; 
    fnt->n1 = n1; 
    fnt->dim = dim;
    fnt->spl = sf_floatalloc2 (n1+2*nw,dim);
    fnt->slv = sf_spline4_init (n1+2*nw);
    
    return fnt;
}

fint1 fint1_init (int nw /* trace extension */, 
		  int n1 /* trace length */,
		  int mute /* taper length */)
/*< intialize single-function interpolation >*/
{
    fint1 fnt;
    float t;
    int it;
    
    fnt = (fint1) sf_alloc (1, sizeof (*fnt));
    
    fnt->nw = nw; 
    fnt->n1 = n1; 
    fnt->spl = sf_floatalloc (n1+2*nw);
    fnt->slv = sf_spline4_init (n1+2*nw);
    fnt->nt = mute;
    if (mute > 0) {
	fnt->t = sf_floatalloc(mute);
	
	for (it=0; it < mute; it++) {
	    t = sinf(0.5*SF_PI*(it+1.)/(mute+1.));
	    fnt->t[it]= t*t; 
	}
    }
    fnt->ir = -1;
    
    return fnt;
}

void vint1_set (vint1 fnt, float** dat /* [dim][n1] */)
/*< set multi-function grid >*/
{
    int i;
    for (i = 0; i < fnt->dim; i++) {
	extend (fnt->nw,fnt->n1,dat[i],fnt->spl[i]);

	fnt->spl[i][                  0] *= (5./6.);
	fnt->spl[i][fnt->n1+2*fnt->nw-1] *= (5./6.);

	sf_tridiagonal_solve (fnt->slv,fnt->spl[i]);
    }
}

void fint1_set (fint1 fnt, float* dat)
/*< set single-function grid >*/
{
    extend (fnt->nw,fnt->n1,dat,fnt->spl);
    fnt->spl[0] *= (5./6.);
    fnt->spl[fnt->n1+2*fnt->nw-1] *= (5./6.);
    sf_tridiagonal_solve (fnt->slv,fnt->spl);
}

void fint1_close (fint1 fnt)
/*< free allocated storage >*/
{
    free (fnt->spl);
    if (fnt->nt > 0) free (fnt->t);
    sf_tridiagonal_close (fnt->slv);
    free (fnt);
}

void vint1_close (vint1 fnt)
/*< free allocated storage >*/
{
    free (fnt->spl[0]);
    free (fnt->spl);
    sf_tridiagonal_close (fnt->slv);
    free (fnt);
}

float fint1_apply (fint1 fnt /* interpolation object */, 
		   int i     /* grid point */, 
		   float x   /* offset */, 
		   bool der  /* to compute derivative */) 
/*< interpolate >*/
{
    float f;
    int j, k;

    if (der) {
	sf_spline4_der (x,fnt->w);
    } else {
	sf_spline4_int (x,fnt->w);
    }
    
    f = 0.;
    for (j = 0; j < 4; j++) {
	k = i+fnt->nw/2+j+1;
	f += fnt->w[j]*fnt->spl[k];
    }
    return f;
}

void vint1_apply (vint1 fnt /* interpolation object */, 
		  int i     /* grid point */, 
		  float x   /* offset */, 
		  bool der  /* to compute derivative */, 
		  float* f  /* output value [dim] */) 
/*< interpolate >*/
{
    int j, k, n;
    
    if (der) {
	sf_spline4_der (x,fnt->w);
    } else {
	sf_spline4_int (x,fnt->w);
    }
    
    for (n=0; n < fnt->dim; n++) {
	f[n] = 0.;
	for (j = 0; j < 4; j++) {
	    k = i+fnt->nw/2+j+1;
	    f[n] += fnt->w[j]*fnt->spl[n][k];
	}
    }
}

void stretch(fint1 str                  /* interpolation object */, 
	     mapfunc map                /* mapping function */,
	     int n1, float d1, float o1 /* old sampling */,
	     int n2, float d2, float o2 /* new sampling */,
	     float *trace               /* new trace [n2] */,
	     float maxstr               /* maximum stretch */)
/*< trace interpolation >*/
{
    int i2, it, im, ip, i;
    float t, tp;

    tp = -1.;
    ip = -1;
    im = str->nt;
    for (i2=0; i2 < n2; i2++) {
	t = o2+i2*d2;
	t = map(t,i2);
	t = (t-o1)/d1;
	it = floorf(t);
	if (it < 0 || it >= n1 || 
	    (tp > 0. && fabsf(t-tp) < maxstr)) { /* too much stretch */
	    trace[i2]=0.;
	    if (ip < 0 || ip != i2-1) {
		for (i=i2-1, im=0; i >=0 && im < str->nt; i--, im++) {
		    trace[i] *= str->t[im];
		}
	    }
	    ip = i2;
	    im=0;
	} else {
	    trace[i2] = fint1_apply(str,it,t-it,false);
	    if (im < str->nt) {
		trace[i2] *= str->t[im];
		im++;
	    }
	}
	tp = t;
    }
}

static const int nw = 3;
static const float a[] = {7./3., -5./3., 1./3.};

void extend (int ne     /* padding */, 
	     int nd     /* data length */, 
	     float *dat /* data [nd] */, 
	     float *ext /* extension [nd+2*ne] */)
/*< 1-D extension >*/
{
    int i, j;
    float s;

    for (i=0; i < nd; i++) {
	ext[ne+i] = dat[i];
    }
    for (i=ne-1; i >= 0; i--) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i+j+1];
	}
	ext[i] = s;
    }
    for (i=nd+ne; i < nd+2*ne; i++) {
	for (s=0., j=0; j < nw; j++) {
	    s += a[j]*ext[i-j-1];
	}
	ext[i] = s;
    }
}

void extend2 (int ne         /* padding */, 
	      int n1, int n2 /* data size */, 
	      float** dat    /* data [n2][n1] */, 
	      float** ext    /* extension [n2+2*ne][n1+2*ne] */, 
	      float* tmp1    /* temporary storage [n2] */, 
	      float* tmp2    /* temporary storage [n2+2*ne] */)
/*< 2-D extension >*/
{
    int i1, i2;
    for (i2=0; i2 < n2; i2++) {
	extend (ne,n1,dat[i2],ext[i2+ne]);
    }
    for (i1=0; i1 < n1+2*ne; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    tmp1[i2] = ext[i2+ne][i1];
	}
	extend (ne,n2,tmp1,tmp2);
	for (i2=0; i2 < n2+2*ne; i2++) {
	    ext[i2][i1] = tmp2[i2];
	} 
    }
}

