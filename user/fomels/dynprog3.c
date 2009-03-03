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

#include <rsf.h>

static float **prev, **next, ****what, **dist, **prob;
static int n1, n2, n3, gt2, gt3;

static float find_minimum(int ic1, int nc1, int jc1,
			  int ic2, int nc2, int jc2,
			  float c, float *pick);
static void interpolate(float f2, float f3, int i1, float *w);

void dynprog3_init(int nz                /* vertical */, 
		   int nx1, int nx2      /* horizontal */, 
		   int gate1, int gate2  /* picking gate */, 
		   float an1, float an2  /* anisotropy */)
/*< Initialize >*/
{
    int i2, i3;

    n1=nz;
    n2=nx1;
    n3=nx2;
    gt2=gate1;
    gt3=gate2;

    prev = sf_floatalloc2(n2,n3);
    next = sf_floatalloc2(n2,n3);
    dist = sf_floatalloc2(n2,n3);
    prob = sf_floatalloc2(2*gate1-1,2*gate2-1);
    what = sf_floatalloc4(2,n2,n3,n1);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    dist[i3][i2] = sqrtf(1. + (i2*i2)/(an1*an1) + (i3*i3)/(an2*an2));
	}
    }
}

void dynprog3_close(void)
/*< deallocate >*/
{
    free(*prev);
    free(prev);
    free(*next);
    free(next);
    free(*dist); 
    free(dist); 
    free(*prob);
    free(prob);
    free(***what);
    free(**what);
    free(*what); 
    free(what);
}

static float find_minimum(int ic1, int nc1, int jc1,
			  int ic2, int nc2, int jc2,
			  float fc, float *pick)
{
    float fm0, f0m, f00, f0p, fp0, fpp, a, b, c, d, e, den, x1, x2, df;
    
    if (0==ic1) {
	ic1++;
    } else if (nc1-1==ic1) {
	ic1--;
    }

    if (0==ic2) {
	ic2++;
    } else if (nc2-1==ic2) {
	ic2--;
    }
    
    f00=prob[ic2][ic1];
    f0p=prob[ic2][ic1+1];
    f0m=prob[ic2][ic1-1];
    fp0=prob[ic2+1][ic1];
    fm0=prob[ic2-1][ic1];

    /* later: select one of four corners */
    fpp=prob[ic2+1][ic1+1];
    
    ic1 += jc1;
    ic2 += jc2;

    a = (fm0 + fp0 - 2*f00)/2.;
    b = (f0m + f0p - 2*f00)/2.;
    c = f00 - f0p - fp0 + fpp;
    d = (fp0 - fm0)/2.;
    e = (f0p - f0m)/2.;

    den = 4*a*b - c*c;
    if (den <= 0. || a <= 0. || b <= 0.) {
	/* no minimum */

	fpp = SF_HUGE;

	if (fpp > f00) {
	    fpp = f00;
	    pick[0] = ic1;
	    pick[1] = ic2;
	}
	if (fpp > f0p) {
	    fpp = f0p;
	    pick[0] = ic1+1;
	    pick[1] = ic2;
	}
	if (fpp > f0m) {
	    fpp = f0m;
	    pick[0] = ic1-1;
	    pick[1] = ic2;
	}
	if (fpp > fp0) {
	    fpp = fp0;
	    pick[0] = ic1;
	    pick[1] = ic2+1;
	}
	if (fpp > fm0) {
	    fpp = fm0;
	    pick[0] = ic1;
	    pick[1] = ic2-1;
	}
	
	return fpp;
    }
	 
    x1 = (c*e-2*b*d)/den;
    x2 = (c*d-2*a*e)/den;

    if (x1 > 1.) {
	pick[0] = ic1+1;
	if (x2 > 1.) {
	    pick[1] = ic2+1;
	    return prob[ic2-jc2+1][ic1-jc1+1];
	} 
	if (x2 < -1.) {
	    pick[1] = ic2-1;
	    return prob[ic2-jc2-1][ic1-jc1+1];
	} 
	pick[1] = ic2;
	return f0p;
    }
    if (x1 < -1.) {
	pick[0] = ic1-1;
	if (x2 > 1.) {
	    pick[1] = ic2+1;
	    return prob[ic2-jc2+1][ic1-jc1-1];
	} 
	if (x2 < -1.) {
	    pick[1] = ic2-1;
	    return prob[ic2-jc2-1][ic1-jc1-1];
	} 
	pick[1] = ic2;
	return f0m;
    }
    if (x2 > 1.) {
	pick[0] = ic1;
	pick[1] = ic2+1;
	return fp0;
    }
    if (x2 < -1.) {
	pick[0] = ic1;
	pick[1] = ic2-1;
	return fm0;
    }

    df = (a*e*e + b*d*d - c*d*e)/den;

    if (f00 < df) {
      pick[0] = ic1;
      pick[1] = ic2;
      return f00;
    }

    f00 -= df;

    pick[0]=ic1+x1;
    pick[1]=ic2+x2;

    return f00;
}

void dynprog3(int i0          /* starting velocity */,
	      float*** weight /* [n1][n2][n3] */)
/*< apply >*/
{
    float d, c, w, w2;
    int i1, i2, i3, k2, k3, ib2, ie2, ib3, ie3, m2, m3, l2, l3;
    
    /* first layer */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    w = 0.5*(weight[1][i3][i2]+weight[0][0][i0]);
	    prev[i3][i2] = dist[i3][SF_ABS(i2-i0)]*w;
	    what[1][i3][i2][0] = i0;
	    what[1][i3][i2][1] = 0;
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
					  c,what[i1][i3][i2]);
	    }
	}

	if (i1==n1-1) return;

	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		prev[i3][i2]=next[i3][i2];
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
