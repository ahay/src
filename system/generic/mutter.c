/* Muting. */
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

#include "mutter.h"

static float dt, t0;
static int nt;
static bool abs0, inner, hyper;              

void mutter_init (int n1 , float o1, float d1 /* time axis */, 
		  bool abs1                   /* if twosided */,
		  bool inner1                 /* inner mute */, 
		  bool hyper1                 /* hyperbolic mute */)
/*< Initialize >*/
{
    nt = n1;
    t0 = o1;
    dt = d1;
    abs0 = abs1;
    inner = inner1;
    hyper = hyper1;
}

void mutter (float tp     /* time step */, 
	     float slope0 /* first slope */, 
	     float slopep /* second slope */, 
	     float x      /* offset */, 
	     float *data  /* trace */,
	     bool  nan   /* nan instaed of zeros*/)
/*< Mute >*/
{
    int it;
    float wt, t;

    if (abs0) x = fabsf(x);

    for (it=0; it < nt; it++) {
	t = t0+it*dt;
	if (hyper) t *= t;
	wt = t - x * slope0;
	if ((inner && wt > 0.) || (!inner && wt < 0.)) {
	    if (nan)
#ifdef NAN
	        data[it]= NAN;
#else
	    	data[it]= 0.0 /0.0;
#endif
	    else
	    	data[it] = 0.;
	    	
	} else {
	    wt = t - tp - x * slopep;
	    if ((inner && wt >=0.) || (!inner && wt <= 0.)) {
		wt = sinf(0.5 * SF_PI * 
			  (t-x*slope0)/(tp+x*(slopep-slope0)));
		data[it] *= (wt*wt);
	    } 
	}
    }
}

/* 	$Id: mutter.c 7107 2011-04-10 02:04:14Z ivlad $	 */

