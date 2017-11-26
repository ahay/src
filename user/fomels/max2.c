/* Maximum in 2-D */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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
#include <rsf.h>

bool not_max(float **dat, int i2, int i1)
/*< check if local maximum >*/
{
    float d;
    int k1, k2;

    d = dat[i2][i1];
    for (k2=-1; k2 <= 1; k2++) {
	for (k1=-1; k1 <= 1; k1++) {
	    if ((k1 || k2) && dat[i2+k2][i1+k1] >= d) return true;
	}
    }
    return false;
}


float max2(float **dat, int i2, int i1, float *xx, float *yy)
/*< find local maximum >*/
{
    float f0,f1,f2, f11,f12,f22;
    float f00,f0m,f0p,fm0,fmm,fmp,fp0,fpm,fpp;
    float det, x, y, df;

    f00=dat[i2][i1];
    fm0=dat[i2][i1-1];
    fp0=dat[i2][i1+1];
    f0m=dat[i2-1][i1];
    fmm=dat[i2-1][i1-1];
    fpm=dat[i2-1][i1+1];
    f0p=dat[i2+1][i1];
    fmp=dat[i2+1][i1-1];
    fpp=dat[i2+1][i1+1];

    f0=(5*f00 + 2*f0m + 2*f0p + 2*fm0 - fmm - fmp + 2*fp0 - fpm - fpp)/9.;
    f1=(-fm0 - fmm - fmp + fp0 + fpm + fpp)/6.;
    f2=(-f0m + f0p - fmm + fmp - fpm + fpp)/6.;
    f11=(-2*f00 - 2*f0m - 2*f0p + fm0 + fmm + fmp + fp0 + fpm + fpp)/3.;
    f12=(fmm - fmp - fpm + fpp)/4.;
    f22=(-2*f00 + f0m + f0p - 2*fm0 + fmm + fmp - 2*fp0 + fpm + fpp)/3.;
    
    det=f12*f12 - f11*f22;
    x=(f1*f22-f2*f12)/det;
    y=(f2*f11-f1*f12)/det;

    df=(f11*f2*f2+f22*f1*f1-2*f1*f12*f2)/(2.*det);
    
    if (x > -1.0f && x < 1.0f && y > -1.0f && y < 1.0f && df > 0) {
	*xx=x;
	*yy=y;
	f0 += df;
    } else {
	*xx=0.0f;
	*yy=0.0f;
    }

    return f0;
}
