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

#include <rsf.h>
/*^*/

static int np=0, n=1;

void sharpen_init(int n1,float perc) 
/*< initialize >*/
{
    n = n1;
    np = n*perc*0.01;
    if (np < 0) np=0;
    if (np >= n) np=n-1;
}

void csharpen(const sf_complex *pp, float *ww) 
/*< compute weight for sharpening regularization >*/
{
    int i, n1;
    float wi, wp;

    for (i=0; i < n; i++) {
	wi = crealf(conjf(pp[i])*pp[i]);
	ww[i]=wi*wi;
    }
    for (n1=np; n1 < n; n1++) {
	wp = sf_quantile(n1,n,ww);
	if (wp > 0.) break;
    }
    for (i=0; i < n; i++) {
	wi = crealf(conjf(pp[i])*pp[i]);
	ww[i] = expf(-wp/wi);
    }
}
