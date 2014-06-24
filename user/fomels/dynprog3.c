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

static float **prev, **next, ****what, **dist, **dist0, **prob, ***ttime;
static int n1, n2, n3, gt2, gt3;

static float find_minimum(int ic1, int nc1, int jc1,
			  int ic2, int nc2, int jc2,
			  float c, float *pick);
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

static float find_minimum(int ic1, int nc1, int jc1,
			  int ic2, int nc2, int jc2,
			  float fc, float *pick)
{
    float f0m, f00, f0p, fpm, fp0, fpp, fmm, fm0, fmp, fmin;
    float den, x1, x2, df;
    
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
    fp0=prob[ic2][ic1+1];
    fm0=prob[ic2][ic1-1];

    f0p=prob[ic2+1][ic1];
    fpp=prob[ic2+1][ic1+1];
    fmp=prob[ic2+1][ic1-1];

    f0m=prob[ic2-1][ic1];
    fpm=prob[ic2-1][ic1+1];
    fmm=prob[ic2-1][ic1-1];
    
    ic1 += jc1;
    ic2 += jc2;

    den = 64*f00*f00 + 32*f00*f0m - 32*f0m*f0m -32*f0p*f0p + 32*f0p*(f00 - 2*f0m) - 32*fm0*fm0 + 
	16*fm0*(2*f00 +5*f0m + 5*f0p) + 7*fmm*fmm - 16*fmm*(4*f00 + f0m + f0p + fm0) + 7*fmp*fmp - 
	2*fmp*(32*f00 + 8*f0m + 8*f0p + 8*fm0 - 25*fmm) - 32*fp0*fp0 + 
	16*fp0*(2*f00 + 5*f0m + 5*f0p - 4*fm0 - fmm - fmp) + 7*fpm*fpm - 
	2*fpm*(32*f00 + 8*f0m + 8*f0p + 8*fm0 - 25*fmm - 7*fmp + 8*fp0) + 7*fpp*fpp - 
	2*fpp*(32*f00 + 8*f0m + 8*f0p + 8*fm0 - 7*fmm - 25*fmp + 8*fp0 - 25*fpm);

    if (den <= 0.) {
	/* no minimum */

	fmin = SF_HUGE;

	if (fmin > f00) {
	    fmin = f00;
	    pick[0] = ic1;
	    pick[1] = ic2;
	}
	if (fmin > f0p) {
	    fmin = f0p;
	    pick[0] = ic1+1;
	    pick[1] = ic2;
	}
	if (fmin > f0m) {
	    fmin = f0m;
	    pick[0] = ic1-1;
	    pick[1] = ic2;
	}

	if (fmin > fp0) {
	    fmin = fp0;
	    pick[0] = ic1;
	    pick[1] = ic2+1;
	}
	if (fmin > fpp) {
	    fmin = fpp;
	    pick[0] = ic1+1;
	    pick[1] = ic2+1;
	}
	if (fmin > fpm) {
	    fmin = fpm;
	    pick[0] = ic1-1;
	    pick[1] = ic2+1;
	}

	if (fmin > fm0) {
	    fmin = fm0;
	    pick[0] = ic1;
	    pick[1] = ic2-1;
	}
	if (fmin > fmp) {
	    fmin = fmp;
	    pick[0] = ic1+1;
	    pick[1] = ic2-1;
	}
	if (fmin > fmm) {
	    fmin = fmm;
	    pick[0] = ic1-1;
	    pick[1] = ic2-1;
	}
	
	return fmin;
    }
	 
    x1 = -2*(8*fm0*fm0 + 4*fm0*(2*f00 - f0m - f0p) - fmm*fmm +
	     fmm*(8*f00 - f0m - 7*f0p + 4*fm0) - fmp*fmp + 
	     fmp*(8*f00 - 7*f0m - f0p + 4*fm0 - 14*fmm) - 8*fp0*fp0 - 
	     4*fp0*(2*f00 - f0m - f0p - 3*fmm - 3*fmp) + fpm*fpm - 
	     fpm*(8*f00 - f0m - 7*f0p + 12*fm0 + 4*fp0) + fpp*fpp - 
	     fpp*(8*f00 - 7*f0m - f0p + 12*fm0 + 4*fp0 - 14*fpm))/den;
    x2 = 2*(-8*f00*f0m + 8*f00*f0p - 8*f0m*f0m + 8*f0p*f0p +
	    4*fm0*(f0m - f0p) + fmm*fmm - 
	    fmm*(8*f00 + 4*f0m + 12*f0p - fm0) - fmp*fmp + 
	    fmp*(8*f00 + 12*f0m + 4*f0p - fm0) + 
	    fp0*(4*f0m - 4*f0p + 7*fmm - 7*fmp) + fpm*fpm - 
	    fpm*(8*f00 + 4*f0m + 12*f0p - 7*fm0 - 14*fmm - fp0) - 
	    fpp*fpp + 
	    fpp*(8*f00 + 12*f0m + 4*f0p - 7*fm0 - 14*fmp - fp0))/den;

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

    df = (256*f00*f00*f00/9 - 68*f00*f0m*f0m/3 + 52*f0m*f0m*f0m/9 +
	  52*f0p*f0p*f0p/9 - 68*f0p*f0p*(f00 - f0m)/3 - 
	  4*f0p*(30*f00*f0m - 17*f0m*f0m)/3 + 52*fm0*fm0*fm0/9 - 
	  2*fm0*fm0*(34*f00 + 15*f0m + 15*f0p)/3 + 
	  2*fm0*(32*f00*f0m - 15*f0m*f0m - 15*f0p*f0p + 
		 2*f0p*(16*f00 - 17*f0m))/3 + 10*fmm*fmm*fmm/9 - 
	  fmm*fmm*(20*f00 + 11*f0m + 13*f0p + 11*fm0)/3 + 
	  fmm*(-64*f00*f00 + 24*f00*f0m - 6*f0m*f0m + 10*f0p*f0p + 
	       4*f0p*(10*f00 - f0m) - 6*fm0*fm0 + 
	       fm0*(24*f00 + 53*f0m + 51*f0p))/3 + 10*fmp*fmp*fmp/9 - 
	  fmp*fmp*(20*f00 + 13*f0m + 11*f0p + 11*fm0 - 26*fmm)/3 + 
	  fmp*(-64*f00*f00 + 40*f00*f0m + 10*f0m*f0m - 6*f0p*f0p + 
	       4*f0p*(6*f00 - f0m) - 6*fm0*fm0 + 
	       fm0*(24*f00 + 51*f0m + 53*f0p) + 26*fmm*fmm + 
	       2*fmm*(12*f00 - 16*f0m - 16*f0p - 21*fm0))/3 + 
	  52*fp0*fp0*fp0/9 - 2*fp0*fp0*(34*f00 + 15*f0m + 15*f0p - 
					34*fm0 - 5*fmm - 5*fmp)/3 +
	  fp0*(64*f00*f0m - 30*f0m*f0m - 30*f0p*f0p + 
	       4*f0p*(16*f00 - 17*f0m) + 68*fm0*fm0 - 
	       4*fm0*(30*f00 + 17*f0m + 17*f0p) - 13*fmm*fmm + 
	       fmm*(40*f00 + 51*f0m + 37*f0p - 4*fm0) - 13*fmp*fmp
	       + fmp*(40*f00 + 37*f0m + 51*f0p - 4*fm0 - 70*fmm))/3 + 
	  10*fpm*fpm*fpm/9 - 
	  fpm*fpm*(20*f00 + 11*f0m + 13*f0p + 13*fm0 - 26*fmm - 
		   6*fmp + 11*fp0)/3 + 
	  fpm*(-64*f00*f00 + 24*f00*f0m - 6*f0m*f0m +
	       10*f0p*f0p + 4*f0p*(10*f00 - f0m) + 10*fm0*fm0 + 
	       fm0*(40*f00 + 51*f0m + 37*f0p) + 26*fmm*fmm + 
	       2*fmm*(12*f00 - 21*f0m - 35*f0p - 16*fm0) + 6*fmp*fmp 
	       - 8*fmp*(f00 + 2*f0m + 2*f0p + 2*fm0 - 3*fmm) - 
	       6*fp0*fp0 + fp0*(24*f00 + 53*f0m + 51*f0p - 4*fm0 - 
				32*fmm - 16*fmp))/3 + 
	  10*fpp*fpp*fpp/9 - fpp*fpp*(20*f00 + 13*f0m + 11*f0p +
				      13*fm0 - 6*fmm - 26*fmp + 
				      11*fp0 - 26*fpm)/3 + 
	  fpp*(-64*f00*f00 + 40*f00*f0m + 10*f0m*f0m - 6*f0p*f0p + 
	       4*f0p*(6*f00 - f0m) + 10*fm0*fm0 + 
	       fm0*(40*f00 + 37*f0m + 51*f0p) + 6*fmm*fmm -
	       8*fmm*(f00 + 2*f0m + 2*f0p + 2*fm0) + 26*fmp*fmp + 
	       2*fmp*(12*f00 - 35*f0m - 21*f0p - 16*fm0 + 12*fmm) - 
	       6*fp0*fp0 + fp0*(24*f00 + 51*f0m + 53*f0p - 4*fm0 - 
				16*fmm - 32*fmp) + 26*fpm*fpm +
	       2*fpm*(12*f00 - 16*f0m - 16*f0p - 35*fm0 + 12*fmm + 
		      12*fmp - 21*fp0))/3)/den;

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
					  c,what[i1][i3][i2]);
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
