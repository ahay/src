/* Pseudo-random numbers: uniform and normally distributed */
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
#include <math.h>

#include <rsf.h>
/*^*/

#include "randn.h"

static const double s = 0.449871, tt = -0.386595;
static const double a = 0.19600, b = 0.25472;
static const double r1 = 0.27597, r2 = 0.27846;
static const double two_pi = 2*3.14159265358979323846264338327950288;

static int    iset = 0;
static float  vset = 0.0;

float randn_one_bm (void)
/*< return a random number (normally distributed, Box-Muller method) >*/
{
    double x1, x2, y1, y2, z1, z2;
    
    if (iset ==0)
    {
        x1 = genrand_real1 ();
        x2 = genrand_real1 ();
        
        z1 = sqrt(-2.0*log(x1));
        z2 = two_pi*x2;
        
        y1 = z1*cos(z2);
        y2 = z1*sin(z2);
        
        iset = 1;
        vset = y1;
        return ((float) y2);
    }
    else
    {
        iset = 0;
        return (vset);
    }
}

float randn_one (void)
/*< return a random number (normally distributed) >*/
{
    double u, v, x, y, q;
    
    do {
	u = genrand_real1 ();
	v = genrand_real1 ();
	v = 1.7156 * (v - 0.5);
	/*     Evaluate the quadratic form */
	x = u - s;
	y = fabsf(v) - tt;
	q = x*x + y*(a*y - b*x);
	/*     Accept P if inside inner ellipse */
	if (q < r1) break;
	/*     Reject P if outside outer ellipse */
	if (q > r2) continue;
	/*     Reject P if outside acceptance region */
    } while (v*v >= -4.0*logf(u)*u*u);
    return ((float) (v/u));
}

void randn (int nr, float *r /* [nr] */)
/*< fill an array with normally distributed numbers >*/
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = randn_one ();
    }
}

void random0 (int nr, float *r /* [nr] */)
/*< fill an array with uniformly distributed numbers >*/
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = (float) genrand_real1 ();
    }
}

/* 	$Id$	 */
