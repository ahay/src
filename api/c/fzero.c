/* Zero finder. */
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

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "_bool.h"
/*^*/

#include "fzero.h"
#include "_defs.h"
#include "error.h"
#include "c99.h"

float sf_zero (float (*func)(float) /* function f(x) */, 
	       float a, float b     /* interval */, 
	       float fa, float fb   /* f(a) and f(b) */,
	       float toler          /* tolerance */, 
	       bool verb            /* verbosity flag */)
/*< Return c such that f(c)=0 (within tolerance). 
  fa and fb should have different signs. >*/
{
    float c, fc, m, s, p, q, r, e, d;
    char method[256];

    if (0. == fa) return a;
    if (0. == fb) return b;

    if (SF_SIG(fa) == SF_SIG(fb)) 
	sf_error("%s: need different sign for zero finding, "
		 "got %f and %f",__FILE__,fa,fb);

    c = b;  
    e = d = b - a;  

    fc = fb;
    /* Main loop, exit from middle of the loop */
    while (fb != 0.) {
	/* Insure that b is the best result so far, a is the previous
	   value of b, and c is on the opposite of the zero from b. */
	if (SF_SIG(fb) == SF_SIG(fc)) {
	    c = a;  fc = fa;
	    e = d = b - a;  
	}
	if (fabsf(fc) < fabsf(fb)) {
	    a = b;    b = c;    c = a;
	    fa = fb;  fb = fc;  fc = fa;
	}
	
	/* Convergence test and possible exit */
	m = 0.5*(c - b);
	if ((fabsf(m) <= toler) || (fb == 0.0)) return b;
   
	/* Choose bisection or interpolation */
	if ((fabsf(e) < toler) || (fabsf(fa) <= fabsf(fb))) {
	    /* Bisection */
	    e = d = m;
	    if (verb) strcpy(method,"bisection");
	} else {
	    /* Interpolation */
	    s = fb/fa;
	    if (a == c) {
		/* Linear interpolation */
		p = 2.0*m*s;
		q = 1.0 - s;
		if (verb) strcpy(method,"linear interpolation");
	    } else {
		/* Inverse quadratic interpolation */
		q = fa/fc;
		r = fb/fc;
		p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
		q = (q - 1.0)*(r - 1.0)*(s - 1.0);
		if (verb) strcpy(method,"inverse quadratic interpolation");
	    }
	    if (p > 0) {
		q = -q;
	    } else { 
		p = -p;
	    }
	    /* Is interpolated point acceptable */
	    if ((2.0*p < 3.0*m*q - fabsf(toler*q)) && (p < fabsf(0.5f*e*q))) {
		e = d;  d = p/q;
	    } else {
		if (verb) strcpy(method,"interpolation not accepted");
		e = d = m; 
	    }
	} /* Interpolation */
	
	/* Next point */
	a = b;
	fa = fb;
	if (fabsf(d) > toler) {
	    b += d;
	} else if (b > c) {
	    b -= toler;
	} else { 
	    b += toler;
	}
	fb = func(b);
	
	if (verb) fprintf(stderr,"(%g,%g,%g) fb=%e method=%s\n",a,b,c,fb,method);
    } /* Main loop */
    
    return b;
}

/* 	$Id$	 */
