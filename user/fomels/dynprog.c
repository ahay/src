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

static float **prev, **next, *dist;
static int n1, n2, ***what;

void dynprog_init(int nz   /* vertical */, 
		  int nx   /* horizontal */, 
		  float an /* anisotropy */)
/*< Initialize >*/
{
    int i2;

    n1=nz;
    n2=nx;

    prev = sf_floatalloc2(n2,n2);
    next = sf_floatalloc2(n2,n2);
    dist = sf_floatalloc(n2);
    what = sf_intalloc3(n2,n2,n1);

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
    free(**what); free(*what); free(what);
}

void dynprog(int gate       /* picking gate */, 
	     float** weight /* [n1][n2] */,
	     int *traj      /* [n1] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, is, i2, i, ic, ib, ie;
    
    for (i2=0; i2 < n2; i2++) {
	w = weight[1][i2];
	for (is=0; is < n2; is++) {
	    w2 = 0.5*(w+weight[0][is]);
	    prev[i2][is] = dist[SF_ABS(i2-is)]*w2;
	    what[1][i2][is] = is;
	}
    }

    for (i1=2; i1 < n1; i1++) {
	sf_warning("step %d of %d",i1,n1);
	for (i2=0; i2 < n2; i2++) {
	    w = weight[i1][i2];
	    ib = SF_MAX(i2-gate,-1);
	    ie = SF_MIN(i2+gate,n2);
	    for (is=0; is < n2; is++) {
		c = FLT_MAX;
		ic = -1;
		for (i=ib+1; i < ie; i++) {
		    w2 = 0.5*(w+weight[i1-1][i]);
		    d = dist[SF_ABS(i2-i)]*w2+prev[i][is];
		    if (d < c) {
			c = d;
			ic = i;
		    }
		}
		next[i2][is]=c;
		what[i1][i2][is]=ic;
	    }
	}
	if (i1==n1-1) break;
	for (i2=0; i2 < n2; i2++) {
	    for (is=0; is < n2; is++) {
		prev[i2][is]=next[i2][is];
	    }
	}
    }

    c = FLT_MAX;
    ic = i = -1;
    for (i2=0; i2 < n2; i2++) {
	for (is=0; is < n2; is++) {
	    d = next[i2][is];
	    if (d < c) {
		c = d;
		ic = i2;
		i = is;
	    }
	}
    }
    for (i1=n1-1; i1 >= 0; i1--) {
	traj[i1]=ic;
	ic = what[i1][ic][i];
    }
}
