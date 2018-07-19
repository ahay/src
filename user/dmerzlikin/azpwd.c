/* Azimuthal Plane-Wave Destruction Interface  */
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
#include "allp3.h"

static int n123;
static float *pwdx, *pwdy, *pwddata, *az;
static allpass ap, aq;

void azpwd_init(int n1, int n2, int n3  /* data size */, 
		 float d1, float d2, float d3 /* sampling */,
                 float o1, float o2, float o3 /* origins */,
		 int nw /* [1,2,3] accuracy order */,
		 float *dip1, float *dip2 /* dip distribution */,
		 float *faz /* azimuth distribution */,
                 int nj1, int nj2 /* antialiasing */)
/*< Initialize AzPWD  >*/
{

    n123 = n1*n2*n3;

    /* iline */
    ap = allpass_init(nw, nj1, n1,n2,n3, dip1);

    /* xline */
    aq = allpass_init(nw, nj2, n1,n2,n3, dip2);

    /* iline + xline */
    allpass3_init(ap,aq);

    pwdx = sf_floatalloc(n123);
    pwdy = sf_floatalloc(n123);
    pwddata = sf_floatalloc(n123);

    az = faz;

}

void azpwd_lop (bool adj, bool add, int fnx, int fny, float* x, float* y)
/*< Apply AzPWD >*/
{

    float conv;
    int i;

    if((n123 != fnx) || (n123 != fny)) sf_error("Wrong dimensions");

    sf_adjnull(adj,add,fnx,fny,x,y);

    for (i=0; i < n123; i++){

	pwdx[i] = 0.0;
	pwdy[i] = 0.0;
	pwddata[i] = 0.0;

    }

    conv = 3.14159265 / 180;

    if (adj == false){

	allpass32_lop(adj,false,n123,n123,x,pwdx,x,pwdy);

    	for (i=0; i < n123; i++){

		pwdx[i] = pwdx[i]*(cosf(conv*az[i]));

		pwdy[i] = pwdy[i]*(sinf(conv*az[i]));			

		y[i] += pwdx[i] + pwdy[i];

	}
	
    } else {

	for (i=0; i < n123; i++){

		pwddata[i] = y[i];

		y[i] = y[i]*cosf(conv*az[i]);		

		pwddata[i] = pwddata[i]*sinf(conv*az[i]);

	}

	allpass32_lop(adj,false,n123,n123,pwdx,y,pwdy,pwddata);


	for (i=0; i < n123; i++){

		x[i] += pwdx[i] + pwdy[i];

	}

    }/* AzPWD */

}

void azpwd_close(void)
/*< free allocated storage >*/
{
    free(pwdx);
    free(pwdy);
    free(pwddata);
    sf_warning("azpwd: memory free");
}
