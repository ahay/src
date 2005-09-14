/* Dynamical programming */
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
#include <float.h>
#include <stdlib.h>

#include <rsf.h>

static float *prev, *next, **what, *dist, *prob;
static int n1, n2, gt;

static float find_minimum(int ic, int nc, int jc, float c, float *pick);
static float interpolate(float fc, int i1);

void dynprog_init(int nz    /* vertical */, 
		  int nx    /* horizontal */, 
		  int gate /* picking gate */, 
		  float an  /* anisotropy */)
/*< Initialize >*/
{
    int i2;

    n1=nz;
    n2=nx;
    gt=gate;

    prev = sf_floatalloc(n2);
    next = sf_floatalloc(n2);
    dist = sf_floatalloc(n2);
    prob = sf_floatalloc(2*gate-1);
    what = sf_floatalloc2(n2,n1);

    for (i2=0; i2 < n2; i2++) {
	dist[i2] = hypotf(i2,an);
    }
}

void dynprog_close(void)
/*< deallocate >*/
{
    free(prev);
    free(next);
    free(dist); 
    free(prob);
    free(*what); 
    free(what);
}

static float find_minimum(int ic, int nc, int jc, float c, float *pick)
{
    float fm, f0, fp, a, b;
    
    if (0==ic) {
	ic++;
	fm=c;
	f0=prob[ic];
	fp=prob[ic+1];
    } else if (nc-1==ic) {
	ic--;
	fm=prob[ic-1];
	f0=prob[ic];
	fp=c;
    } else {
	fm=prob[ic-1];
	f0=c;
	fp=prob[ic+1];
    }
    ic += jc;
    a = fm+fp-2.*f0;
    if (a <= 0.) {
	if (fm < f0 && fm < fp) {
	    *pick = ic-1;
	    return fm;
	} 
	if (fp < f0 && fp < fm) {
	    *pick = ic+1;
	    return fp;
	} 
	*pick = ic;
	return f0;
    }
    b = 0.5*(fm-fp);
    a = b/a;
    if (a > 1.) {
	*pick = ic+1;
	return fp;
    }
    if (a < -1.) {
	*pick = ic-1;
	return fm;
    }
    if (f0 < 0.5*b*a) {
      *pick = ic;
      return f0;
    }
    f0 -= 0.5*b*a;
    *pick=ic+a;
    return f0;
}

void dynprog(int i0 /* starting velocity */,
	     float** weight /* [n1][n2] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, i2, i, ic, ib, ie, it;
    
    for (i2=0; i2 < n2; i2++) {
	w = 0.5*(weight[1][i2]+weight[0][i0]);
	prev[i2] = dist[SF_ABS(i2-i0)]*w;
	what[1][i2] = i0;
    }

    for (i1=2; i1 < n1; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    w = weight[i1][i2];
	    ib = SF_MAX(i2-gt,-1);
	    ie = SF_MIN(i2+gt,n2);
	    c = FLT_MAX;
	    ic = -1;
	    for (i=ib+1; i < ie; i++) {
		w2 = 0.5*(w+weight[i1-1][i]);
		d = dist[SF_ABS(i2-i)]*w2+prev[i];
		it = i-ib-1;
		if (d < c) {
		    c = d;
		    ic = it;
		}
		prob[it]=d;
	    }
		
	    next[i2]=find_minimum(ic,ie-ib-1,ib+1,c,&what[i1][i2]);
	}
	if (i1==n1-1) return;
	for (i2=0; i2 < n2; i2++) {
	    prev[i2]=next[i2];
	}
    }
}

void dynprog_traj(float *traj /* [n1] */)
/*< extract trajectory >*/
{
    int i2, i1;
    float c, d, fc;

    c = FLT_MAX;
    fc = 0;
    for (i2=0; i2 < n2; i2++) {
	d = next[i2];
	if (d < c) {
	    c = d;
	    fc = i2;
	}
    }
    for (i1=n1-1; i1 >= 0; i1--) {
	traj[i1]=fc;
	fc = interpolate(fc,i1);
    }
}

static float interpolate(float fc, int i1)
{
    int ic;

    ic = floorf(fc);
    fc -= ic;
    if (n2-1 <= ic) return what[i1][n2-1];
    if (0 > ic) return what[i1][0];

    fc = what[i1][ic]*(1.-fc)+what[i1][ic+1]*fc;
    return fc;
}
