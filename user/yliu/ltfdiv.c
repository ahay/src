/* Streaming local time-frequency decomposition. */
/*
  Copyright (C) 2025 Jilin University
  
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
#include <stdlib.h>
#include <string.h>

static int n1=1, n2=1, nw=1;
static float *eps, lam=1, tol=10;
static float **arr1, **arr2;
static bool center=false, smooth=false;
static sf_triangle *tr;


void ltfdiv_init(int ndat        /* input trace length */, 
		         int nbase       /* base length [nbase=2*nw]*/,
                 int *rects      /* nonstationary smoothing radius [nw] */,
                 float *eps0     /* nonstationary localization scalar [nw] */,
                 float lam0     /* regularization */,
                 float tol0     /* tolerance of ratio correction */,
                 bool iscenter   /* bidirectional flag */,
		         bool issmooth   /* smooth flag (anti-leackage) */)
/*< initialize >*/
{
    int i2;
    n1 = ndat;
    n2 = nbase;
    nw = n2/2;
    eps = eps0;
    lam = lam0;
    tol = tol0;
    center = iscenter;
    smooth = issmooth;
    if (!rects) smooth = false;
    if (smooth){
        tr = sf_alloc(nw, sizeof(*tr));
        for(i2=0;i2<nw;i2++){
            tr[i2] = sf_triangle_init(rects[i2],n1,false);
        }
    }
    arr1 = sf_floatalloc2(n1,n2);
    arr2 = sf_floatalloc2(n1,n2);
    memset(*arr1, 0, sizeof(float)*n1*n2);
    memset(*arr2, 0, sizeof(float)*n1*n2);

}


void ltfdiv_close(void)
/*< close >*/
{
    int i2;
    if (smooth){
        for(i2=0;i2<nw;i2++){
            sf_triangle_close(tr[i2]);
        }
    }
    if (center){
        free(arr2);
    }
    free(arr1);
}


void ltfdiv(float *num, /* input data [ndat] */
            float **den, /* bases [ndat,nbase] */
            float **rat    /* coeffs [ndat,nbase] */)
/*< Streaming local time-frequency decomposition. >*/
{
    int i1, i2;
    float den0 = 0, num0= 0, num1= 0;
    float tmp0 =0, tmp1=0;

    memset(*rat, 0, sizeof(float)*n1*n2);

    for (i2=0; i2<nw; i2++){
        den0 += den[i2][0]*den[i2][0];
        den0 += den[i2+nw][0]*den[i2+nw][0];

    }
    den0 += lam;

    for (i1=0; i1<n1; i1++){
        num0 = 0;
        num1 = 0;
        tmp0 = num[i1]/den0;
        if (center) tmp1 = num[n1-1-i1]/den0;
        if (i1>0)
	    for (i2=0; i2<n2; i2++){
		num0 += eps[i2%nw]*den[i2][i1]*arr1[i2][i1-1]/den0;
		if (center)
		    num1 += eps[i2%nw]*den[i2][n1-1-i1]*arr2[i2][n1-i1]/den0;
	    }
        for (i2=0; i2<n2; i2++){
	    arr1[i2][i1] = den[i2][i1]*(tmp0-num0)+eps[i2%nw]*arr1[i2][i1-1];
	    
            if(center){
                arr2[i2][n1-1-i1] = den[i2][n1-1-i1]*(tmp1-num1)+eps[i2%nw]*arr2[i2][n1-i1];
                rat[i2][n1-1-i1] += 0.5*arr2[i2][n1-1-i1];
                rat[i2][i1] += 0.5*arr1[i2][i1];
            } else rat[i2][i1] += arr1[i2][i1];
        }
    }
    
    if(smooth){
        for (i2=0; i2<nw; i2++){
            sf_smooth(tr[i2], 0, 1, false, rat[i2]);
            sf_smooth(tr[i2], 0, 1, false, rat[i2+nw]);
        }
    }
    for (i1=0; i1<n1; i1++){
        num0 =0;
        for (i2=0; i2<nw; i2++){
            num0 += den[i2][i1]*rat[i2][i1]+den[i2+nw][i1]*rat[i2+nw][i1];
        }
        for (i2=0; i2<n2; i2++) {
            if (num0 != 0 && SF_ABS(num[i1]/num0)<tol) {
                rat[i2][i1] *= num[i1]/num0;
            } 
        }
    }
}
