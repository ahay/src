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
#include <math.h>

#include "randn.h"
#include "_defs.h"
#include "_bool.h"
#include "mt19937ar.h"

static bool   iset = true;
static float  vset = 0.0;

float sf_randn_one_bm (void)
/*< return a random number (normally distributed, Box-Muller method) >*/
{
    double x1, x2, y1, y2, z1, z2;
    
    if (iset) {
        x1 = genrand_real1 ();
        x2 = genrand_real1 ();
        
        z1 = sqrt(-2.0*log(x1));
        z2 = 2.0*SF_PI*x2;
        
        y1 = z1*cos(z2);
        y2 = z1*sin(z2);
        
        iset = false;
        vset = y1;
        return ((float) y2);
    } else {
        iset = true;
        return vset;
    }
}

void sf_randn (int nr, float *r /* [nr] */)
/*< fill an array with normally distributed numbers >*/
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = sf_randn_one_bm ();
    }
}

void sf_random (int nr, float *r /* [nr] */)
/*< fill an array with uniformly distributed numbers >*/
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = (float) genrand_real1 ();
    }
}

/* 	$Id: randn.c 1768 2006-03-20 09:41:07Z fomels $	 */
