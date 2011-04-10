/* Triangle smoothing along slope */
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
/*^*/

#include "trisl.h"

static int n1, n2, rect;
static float **p, **tmp, *trace1, *trace2, amp;
static sf_triangle tr;

static void forw(int i2, const float* t1, float* t2);
static void back(int i2, float* t1, const float* t2);

void trisl_init(int m1, int m2       /* data dimensions */, 
		int rect1, int rect2 /* triangle radius */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    rect = rect2;

    tmp = sf_floatalloc2(n1,n2+2*rect);

    trace1 = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n1);

    tr = sf_triangle_init (rect1,n1);
    amp = 1./rect;
}

void trisl_set(float** p1 /* slope [m2][m1] */)
/*< set slopes >*/
{
    p = p1;
}

void trisl_close(void)
/*< free allocated storage >*/
{
    free(*tmp);
    free(tmp);
    free(trace1);
    free(trace2);
    sf_triangle_close(tr);
}

static void forw(int i2, const float* t1, float* t2)
{
    int i1, it;
    float t;

    if (i2 < rect-1) {
	back(2*(rect-1)-i2,t2,t1);
    } else if (i2 > n2+rect-1) {
	back(2*(n2+rect-1)-i2,t2,t1);
    } else if (i2 == rect-1 || i2 == n2+rect-1) {
	for (i1=0; i1 < n1; i1++) {
	    t2[i1] += t1[i1];
	}
    } else {
	for (i1=0; i1 < n1; i1++) {
	    t = i1 + p[i2-rect][i1];
	    if (t < 0.) continue;
	    it = t;
	    if (it >=0 && it < n1) {
		t -= it;
		if (it < n1-1) {
		    t2[it]   += t1[i1]*(1.-t);
		    t2[it+1] += t1[i1]*t;
		} else if (t==0.) {
		    t2[it]   += t1[i1];
		}
	    }
	}
    }
}

static void back(int i2, float* t1, const float* t2)
{
    int i1, it;
    float t;

    if (i2 < rect-1) {
	forw(2*(rect-1)-i2,t2,t1);
    } else if (i2 > n2+rect-1) {
	forw(2*(n2+rect-1)-i2,t2,t1);
    } else if (i2 == rect-1 || i2 == n2+rect-1) {
	for (i1=0; i1 < n1; i1++) {
	    t1[i1] += t2[i1];
	}	
    } else {
	for (i1=0; i1 < n1; i1++) {	    
	    t = i1 + p[i2-rect][i1];
	    if (t < 0.) continue;
	    it = t;
	    if (it >=0 && it < n1) {
		t -= it;
		if (it < n1-1) {
		    t1[i1] += t2[it]*(1.-t) + t2[it+1]*t;
		} else if (t==0.) {
		    t1[i1] += t2[it];
		}
	    }
	}
    }
}

static void roll(bool adj, float** in)
{
    int i2;

    if (adj) {
	for (i2=n2+2*rect-2; i2 >=0; i2--) {	    
	    back(i2,in[i2],in[i2+1]);
	}
    } else {
	for (i2=0; i2 < n2+2*rect-1; i2++) {
	    forw(i2,in[i2],in[i2+1]);
	}
    }
}

static void shift(bool adj, float** in)
{
    int i1, i2, ir;

    if (adj) {
	for (i2=0; i2 < n2+rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		trace1[i1] = -in[i2+rect][i1];
	    }
	    
	    for (ir=i2+rect-1; ir >= i2; ir--) {
		for (i1=0; i1 < n1; i1++) {
		    trace2[i1] = trace1[i1];
		    trace1[i1]=0.;
		}
		back(ir,trace1,trace2);
	    }
	    
	    for (i1=0; i1 < n1; i1++) {
		in[i2][i1] += trace1[i1];
	    }
	}
    } else {
	for (i2=n2+rect-1; i2 >=0; i2--) {
	    for (i1=0; i1 < n1; i1++) {
		trace1[i1] = in[i2][i1];
	    }
	    
	    for (ir=i2; ir < i2+rect; ir++) {
		forw(ir,trace1,trace2);
		for (i1=0; i1 < n1; i1++) {
		    trace1[i1] = trace2[i1];
		    trace2[i1]=0.;
		}
	    }
	    
	    for (i1=0; i1 < n1; i1++) {
		in[i2+rect][i1] -= trace1[i1];
	    }
	}
    }
}


void trisl_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i1, i2;
    float *trace;

    if (nx != n1*n2 || ny != nx) 
	sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj, add, nx, ny, x, y);

    for (i2=0; i2 < n2; i2++) {
	trace = tmp[i2+rect];
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = adj? y[i2*n1+i1]: x[i2*n1+i1];
	}
	if (!adj) sf_smooth (tr,0,1,false,false,trace);
    }

    if (adj) {
	for (i2=0; i2 < rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[n2+rect+i2][i1] = 0.;
		tmp[rect-1-i2][i1]  = 0.;
	    }
	}
    } else {
	for (i2=0; i2 < rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[n2+rect+i2][i1] = tmp[n2+rect-1-i2][i1];
		tmp[rect-1-i2][i1]  = tmp[rect+i2][i1];
	    }
	}
    }

    for (i2=0; i2 < n2+2*rect; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    tmp[i2][i1] *= amp;
	}
    }

    roll(false,tmp);
    shift(false,tmp);
    shift(true,tmp);
    roll(true,tmp);
    
    for (i2=0; i2 < n2+2*rect; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    tmp[i2][i1] *= amp;
	}
    }
    
    if (adj) {
	for (i2=0; i2 < rect; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tmp[n2+rect-1-i2][i1] += tmp[n2+rect+i2][i1];
		tmp[rect+i2][i1]      += tmp[rect-1-i2][i1];
	    }
	}
    }

    for (i2=0; i2 < n2; i2++) {
	trace = tmp[i2+rect];
	if (adj) sf_smooth2 (tr,0,1,false,false,trace);
	for (i1=0; i1 < n1; i1++) {
	    if (adj) {
		x[i2*n1+i1] += trace[i1];
	    } else {
		y[i2*n1+i1] += trace[i1];
		
	    }
	}
    }
}
