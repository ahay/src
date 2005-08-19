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

#include <rsf.h>

static float **prev, **next, ***what, *dist, *time;
static int n1, n2, gt;

static float find_minimum(int ic, int nc, int jc, float c, float *pick);
static float interpolate(float fc, int i1, int i);

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

    prev = sf_floatalloc2(n2,n2);
    next = sf_floatalloc2(n2,n2);
    dist = sf_floatalloc(n2);
    time = sf_floatalloc(2*gate-1);
    what = sf_floatalloc3(n2,n2,n1);

    for (i2=0; i2 < n2; i2++) {
	dist[i2] = hypotf(i2,an);
    }
}

void dynprog_close(void)
/*< deallocate >*/
{
    free(*prev); free(prev);
    free(*next); free(next);
    free(dist); 
    free(time);
    free(**what); free(*what); free(what);
}

static float find_minimum(int ic, int nc, int jc, float c, float *pick)
{
    float fm, f0, fp, a, b;
    
    if (0==ic) {
	ic++;
	fm=c;
	f0=time[ic];
	fp=time[ic+1];
    } else if (nc-1==ic) {
	ic--;
	fm=time[ic-1];
	f0=time[ic];
	fp=c;
    } else {
	fm=time[ic-1];
	f0=c;
	fp=time[ic+1];
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
    f0 -= 0.5*b*a;
    *pick=ic+a;
    return f0;
}

void dynprog(float** weight /* [n2][n1] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, is, i2, i, ic, ib, ie, it;
    
    for (i2=0; i2 < n2; i2++) {
	w = weight[i2][1];
	for (is=0; is < n2; is++) {
	    w2 = 0.5*(w+weight[0][is]);
	    prev[i2][is] = dist[SF_ABS(i2-is)]*w2;
	    what[1][i2][is] = is;
	}
    }

    for (i1=2; i1 < n1; i1++) {
	sf_warning("step %d of %d",i1,n1);
	for (i2=0; i2 < n2; i2++) {
	    w = weight[i2][i1];
	    ib = SF_MAX(i2-gt,-1);
	    ie = SF_MIN(i2+gt,n2);
	    for (is=0; is < n2; is++) {
		c = FLT_MAX;
		ic = -1;
		for (i=ib+1; i < ie; i++) {
		    w2 = 0.5*(w+weight[i][i1-1]);
		    d = dist[SF_ABS(i2-i)]*w2+prev[i][is];
		    it = i-ib-1;
		    if (d < c) {
			c = d;
			ic = it;
		    }
		    time[it]=d;
		}
		
		next[i2][is]=find_minimum(ic,ie-ib-1,ib+1,c,&what[i1][i2][is]);
	    }
	}
	if (i1==n1-1) return;
	for (i2=0; i2 < n2; i2++) {
	    for (is=0; is < n2; is++) {
		prev[i2][is]=next[i2][is];
	    }
	}
    }
}

void dynprog_traj(float o2, float d2, float *traj /* [n1] */)
/*< extract trajectory >*/
{
    int ic, i, i2, is, i1;
    float c, d, fc;

    c = FLT_MAX;
    ic = i = -1;
    for (i2=0; i2 < n2; i2++) {
	for (is=0; is < n2; is++) {
	    d = next[i2][is];
	    if (d < c) {
		c = d;
		fc = i2;
		i = is;
	    }
	}
    }
    for (i1=n1-1; i1 >= 0; i1--) {
	traj[i1]=o2+fc*d2;
	fc = interpolate(fc,i1,i);
    }
}

static float interpolate(float fc, int i1, int i)
{
    int ic;

    ic = floorf(fc);
    fc -= ic;
    if (n2-1 <= ic) return what[i1][n2-1][i];
    if (0 > ic) return what[i1][0][i];

    fc = what[i1][ic][i]*(1.-fc)+what[i1][ic+1][i]*fc;
    return fc;
}
