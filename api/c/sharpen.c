/* sharpening */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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

#include "_bool.h"
#include "_defs.h"
#include "adjnull.h"
#include "komplex.h"
#include "alloc.h"
#include "error.h"
#include "quantile.h"
#include "weight.h"

#include "sharpen.h"

static int np=0, n=1;
static float *ww, f;

void sf_sharpen_init(int n1     /* data size */,
		     float perc /* quantile percentage */,
		     float fact /* multiplicative factor */) 
/*< initialize >*/
{
    int i;

    n = n1;
    np = n*perc*0.01;
    if (np < 0) np=0;
    if (np >= n-1) np=n-2;

    ww = sf_floatalloc(n);
    sf_weight_init(ww);

    for (i=0; i < n; i++) {
	ww[i] = 1.0;
    }

    f = fact;
}

void sf_sharpen_close(void)
/*< free allocated storage >*/
{
    free(ww);
}

float sf_sharpen(const float *pp) 
/*< compute weight for sharpening regularization >*/
{
    int i, n1;
    float wp, wmax, wmin, wi;

    wmax = 0.;
    for (i=0; i < n; i++) {
	ww[i] = fabsf(pp[i]);
	if (ww[i] > wmax) wmax=ww[i];
    }
    wmin = FLT_EPSILON*wmax;

    wp = 0.;
    for (n1=np; n1 < n-1; n1++) {
	wp = sf_quantile(n1,n,ww);
	if (wp > wmin) break;
    }
  
    for (i=0; i < n; i++) {
	wi = fabsf(pp[i])+wmin;
	ww[i] = expf(-f*wp*wp/(wi*wi));
    }

    return wp;
}

void sf_csharpen(const sf_complex *pp) 
/*< compute weight for sharpening regularization >*/
{
    int i, n1;
    float wp, wmax, wmin, wi;

    wmax = 0.;
    for (i=0; i < n; i++) {
	ww[i] = cabsf(pp[i]);
	if (ww[i] > wmax) wmax=ww[i];
    }
    wmin = FLT_EPSILON*wmax;

    wp = 0.;
    for (n1=np; n1 < n-1; n1++) {
	wp = sf_quantile(n1,n,ww);
	if (wp > wmin) break;
    }
  
    for (i=0; i < n; i++) {
	wi = cabsf(pp[i])+wmin;
	ww[i] = expf(-f*wp*wp/(wi*wi));
    }
}


