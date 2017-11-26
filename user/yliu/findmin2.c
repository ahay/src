/* Find minimum in 2-D */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

float find_minimum(int ic1, int nc1, int jc1,
		   int ic2, int nc2, int jc2,
		   float fc, float** prob, float *pick)
/*< find minimum >*/
{
    float f0m, f00, f0p, fpm, fp0, fpp, fmm, fm0, fmp;
    float den, den2, x1, x2, df;
    const float eps=0.01;
    
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

    den2 = den*den+eps;
    den = den/den2;
	 
    x1 = -2*(8*fm0*fm0 + 4*fm0*(2*f00 - f0m - f0p) - fmm*fmm +
	     fmm*(8*f00 - f0m - 7*f0p + 4*fm0) - fmp*fmp + 
	     fmp*(8*f00 - 7*f0m - f0p + 4*fm0 - 14*fmm) - 8*fp0*fp0 - 
	     4*fp0*(2*f00 - f0m - f0p - 3*fmm - 3*fmp) + fpm*fpm - 
	     fpm*(8*f00 - f0m - 7*f0p + 12*fm0 + 4*fp0) + fpp*fpp - 
	     fpp*(8*f00 - 7*f0m - f0p + 12*fm0 + 4*fp0 - 14*fpm))*den;
    x2 = 2*(-8*f00*f0m + 8*f00*f0p - 8*f0m*f0m + 8*f0p*f0p +
	    4*fm0*(f0m - f0p) + fmm*fmm - 
	    fmm*(8*f00 + 4*f0m + 12*f0p - fm0) - fmp*fmp + 
	    fmp*(8*f00 + 12*f0m + 4*f0p - fm0) + 
	    fp0*(4*f0m - 4*f0p + 7*fmm - 7*fmp) + fpm*fpm - 
	    fpm*(8*f00 + 4*f0m + 12*f0p - 7*fm0 - 14*fmm - fp0) - 
	    fpp*fpp + 
	    fpp*(8*f00 + 12*f0m + 4*f0p - 7*fm0 - 14*fmp - fp0))*den;

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
		      12*fmp - 21*fp0))/3)*den;

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
