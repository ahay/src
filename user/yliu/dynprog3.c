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
#include <stdlib.h>

#include "findmin2.h"

#include <rsf.h>

static float **prev, **next, ****what, **dist, **dist0, **prob, ***ttime;
static int n1, n2, n3, gt2, gt3;

static void interpolate(float f2, float f3, int i1, float *w);

float *** dynprog3_init(int nz                /* vertical */, 
			int nx1, int nx2      /* horizontal */, 
			int gate1, int gate2  /* picking gate */, 
			float an1, float an2  /* anisotropy */,
			bool savetime         /* save traveltime */)
/*< Initialize >*/
{
    int i2, i3;
    float hyp;

    n1=nz;
    n2=nx1;
    n3=nx2;
    gt2=gate1;
    gt3=gate2;

    prev = sf_floatalloc2(n2,n3);
    next = sf_floatalloc2(n2,n3);
    dist = sf_floatalloc2(n2,n3);
    dist0 = sf_floatalloc2(n2,n3);

    prob = sf_floatalloc2(2*gate1-1,2*gate2-1);
    what = sf_floatalloc4(2,n2,n3,n1);

    ttime = savetime? sf_floatalloc3(n2,n3,n1): NULL;

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    hyp = (i2*i2)/(an1*an1) + (i3*i3)/(an2*an2);
	    dist0[i3][i2] = sqrtf(hyp);
	    dist[i3][i2] = sqrtf(1. + hyp);
	}
    }

    return ttime;
}

void dynprog3_close(void)
/*< deallocate >*/
{
    if (NULL != ttime) {
	free(**ttime);
	free(*ttime);
	free(ttime);
    }
    free(*prev);
    free(prev);
    free(*next);
    free(next);
    free(*dist0); 
    free(dist0);
    free(*dist); 
    free(dist); 
    free(*prob);
    free(prob);
    free(***what);
    free(**what);
    free(*what); 
    free(what);
}


void dynprog3(int q2, int q3  /* starting velocity */,
	      float*** weight /* [n1][n3][n2] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, i2, i3, k2, k3, ib2, ie2, ib3, ie3, m2, m3, l2, l3;
    
    /* first layer */
    if (NULL != ttime) {
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		w = 0.5*(weight[0][i3][i2]+weight[0][q3][q2]);
		ttime[0][i3][i2] = dist0[SF_ABS(i3-q3)][SF_ABS(i2-q2)];
	    }
	}
    }

    /* second layer */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    w = 0.5*(weight[1][i3][i2]+weight[0][q3][q2]);
	    prev[i3][i2] = dist[SF_ABS(i3-q3)][SF_ABS(i2-q2)]*w;
	    what[1][i3][i2][0] = q2;
	    what[1][i3][i2][1] = q3;

	    if (NULL != ttime) ttime[1][i3][i2]=prev[i3][i2];
	}
    }

    for (i1=2; i1 < n1; i1++) {
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		w = weight[i1][i3][i2];
		
		ib2 = SF_MAX(i2-gt2,-1);
		ie2 = SF_MIN(i2+gt2,n2);
		
		ib3 = SF_MAX(i3-gt3,-1);
		ie3 = SF_MIN(i3+gt3,n3);
		
		c = SF_HUGE;
		l2 = -1;
		l3 = -1;
		
		for (k3=ib3+1; k3 < ie3; k3++) {
		    for (k2=ib2+1; k2 < ie2; k2++) {
			
			w2 = 0.5*(w+weight[i1-1][k3][k2]);
			
			/* time = distance*slowness */
			d = dist[SF_ABS(i3-k3)][SF_ABS(i2-k2)]*w2+prev[k3][k2];

			m2 = k2-ib2-1;
			m3 = k3-ib3-1;

			if (d < c) {
			    c = d;
			    l2 = m2;
			    l3 = m3;
			}

			prob[m3][m2]=d;
		    }
		}
		
		next[i3][i2]=find_minimum(l2,ie2-ib2-1,ib2+1,
					  l3,ie3-ib3-1,ib3+1,
					  c,prob,what[i1][i3][i2]);
	    }
	}

	if (i1==n1-1) return;

	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		prev[i3][i2]=next[i3][i2];
		if (NULL != ttime) ttime[i1][i3][i2]=prev[i3][i2];
	    }
	}
    } /* i1 */
}

void dynprog3_traj(float **traj /* [2][n1] */)
/*< extract trajectory >*/
{
    int i3, i2, i1;
    float c, d, f[2];

    c = SF_HUGE;
    f[0] = 0;
    f[1] = 0;

    /* minimum at the bottom */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    d = next[i3][i2];

	    if (d < c) {
		c = d;
		f[0] = i2;
		f[1] = i3;
	    }	    
	}
    }

    /* coming up */
    for (i1=n1-1; i1 >= 0; i1--) {
	traj[0][i1]=f[0];
	traj[1][i1]=f[1];
	
	interpolate(f[0],f[1],i1,f);
    }
}

static void interpolate(float f2, float f3, int i1, float *w)
{
    int i2, i3, i;

    i2 = floorf(f2);
    f2 -= i2;

    if (i2 >= n2-1) {
	i2 = n2-2;
	f2 = 0.;
    } else if (i2 < 0) {
	i2 = 0;
	f2 = 0.;
    }

    i3 = floorf(f3);
    f3 -= i3;

    if (i3 >= n3-1) {
	i3 = n3-2;
	f3 = 0.;
    } else if (i3 < 0) {
	i3 = 0;
	f3 = 0.;
    }

    for (i=0; i < 2; i++) {
	w[i] = 
	    what[i1][i3  ][i2  ][i]*(1.-f2)*(1.-f3) +
	    what[i1][i3+1][i2  ][i]*(1.-f2)*f3 +
	    what[i1][i3  ][i2+1][i]*f2*(1.-f3) +
	    what[i1][i3+1][i2+1][i]*f2*f3;
    }
}
