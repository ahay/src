/* Dynamical programming in 3-D */
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

static float **prev, **next, **dist;
static int ***what;
static int n1, n2, n3, gt2, gt3;

void dynprog2_init(int nz                   /* vertical */, 
		   int nx, int ny           /* horizontal */, 
		   int xgate, int ygate     /* picking gate */, 
		   float anx, float any     /* anisotropy */)
/*< Initialize >*/
{
    int i2, i3;
    float d;

    n1=nz;
    n2=nx;
    n3=ny;
    gt2=xgate;
    gt3=ygate;

    prev = sf_floatalloc2(n2,n3);
    next = sf_floatalloc2(n2,n3);
    dist = sf_floatalloc2(n2,n3);
    what = sf_intalloc3(n2,n3,n1);

    for (i3=0; i3 < n3; i3++) {
	d = hypotf(any*i3,anx);
	for (i2=0; i2 < n2; i2++) {
	    dist[i3][i2] = hypotf(i2,d);
	}
    }
}

void dynprog2_close(void)
/*< deallocate >*/
{
    free(*prev); free(prev);
    free(*next); free(next);
    free(*dist); free(dist); 
    free(**what); free(*what); free(what);
}

void dynprog2(int i0 /* starting velocity */,
	      float*** weight /* [n1][n3][n2] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, i2, i3, j2, j3, ic, ib3, ie3, ib2, ie2;
    
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    w = 0.5*(weight[1][i3][i2]+weight[0][i3][i0]);
	    prev[i3][i2] = dist[0][SF_ABS(i2-i0)]*w;
	    what[1][i3][i2] = i0;
	}
    }

    for (i1=2; i1 < n1; i1++) {
	sf_warning("depth %d of %d",i1+1,n1);

	for (i3=0; i3 < n3; i3++) {
	    ib3 = SF_MAX(i3-gt3,-1);
	    ie3 = SF_MIN(i3+gt3,n3);
	    for (i2=0; i2 < n2; i2++) {
		w = weight[i1][i3][i2];
		ib2 = SF_MAX(i2-gt2,-1);
		ie2 = SF_MIN(i2+gt2,n2);
		c = FLT_MAX;
		ic = -1;

		for (j3=ib3+1; j3 < ie3; j3++) {
		    for (j2=ib2+1; j2 < ie2; j2++) {
			w2 = 0.5*(w+weight[i1-1][j3][j2]);
			d = dist[SF_ABS(i3-j3)][SF_ABS(i2-j2)]*w2+prev[j3][j2];
			
			if (d < c) {
			    c = d;
			    ic = j2-ib2-1;
			}
		    }
		}
		
		what[i1][i3][i2] = ic;
		next[i3][i2]=c;
	    }
	}
	if (i1==n1-1) return;
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		prev[i3][i2]=next[i3][i2];
	    }
	}
    }
}

void dynprog2_traj(int **traj /* [n3][n1] */)
/*< extract trajectory >*/
{
    int i3, i2, i1, ic;
    float c, d;

    for (i3=0; i3 < n3; i3++) {
	c = FLT_MAX;
	ic = 0;
	for (i2=0; i2 < n2; i2++) {
	    d = next[i3][i2];
	    if (d < c) {
		c = d;
		ic = i2;
	    }
	}
	for (i1=n1-1; i1 >= 0; i1--) {
	    traj[i3][i1]=ic;
	    ic = what[i1][i3][ic];
	}
    }
}
