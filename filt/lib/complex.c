/* Complex number operations */
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

#include "complex.h"
#include "c99.h"

#include "kiss_fft.h"
/*^*/

kiss_fft_cpx sf_csqrtf (kiss_fft_cpx c)
/*< replacement for csqrtf >*/
{
    float d, r, s;
    kiss_fft_cpx v;

    if (c.i == 0) {
      if (c.r < 0) {
	  v.r = 0.;
	  v.i = copysignf (sqrtf (-c.r), c.i);
      } else {
	  v.r =  fabsf (sqrtf (c.r));
	  v.i =  copysignf (0, c.i);
      }
    } else if (c.r == 0) {
	r = sqrtf (0.5 * fabsf (c.i));
	v.r = r;
	v.i = copysignf (r, c.i);
    } else {
	d = hypotf (c.r, c.i);
	/* Use the identity   2  Re res  Im res = Im x
	   to avoid cancellation error in  d +/- Re x.  */
	if (c.r > 0) {
	    r = sqrtf (0.5 * d + 0.5 * c.r);
	    s = (0.5 * c.i) / r;
        } else {
	    s = sqrtf (0.5 * d - 0.5 * c.r);
	    r = fabsf ((0.5 * c.i) / s);
        }
	v.r = r;
	v.i = copysignf (s, c.i);
    }
    return v;
}

kiss_fft_cpx sf_cdiv(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex division >*/
{
    kiss_fft_cpx c;
    float r,den;
    if (fabsf(b.r)>=fabsf(b.i)) {
	r = b.i/b.r;
	den = b.r+r*b.i;
	c.r = (a.r+r*a.i)/den;
	c.i = (a.i-r*a.r)/den;
    } else {
	r = b.r/b.i;
	den = b.i+r*b.r;
	c.r = (a.r*r+a.i)/den;
	c.i = (a.i*r-a.r)/den;
    }
    return c;
}

#define SF_CMUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
/*^*/

#define SF_CCMUL(m,a,b) \
    do{ (m).r = (a).r*(b).r + (a).i*(b).i;\
        (m).i = (a).i*(b).r - (a).r*(b).i; }while(0)
/*^*/


