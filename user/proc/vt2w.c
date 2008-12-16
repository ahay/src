/* Compute the Stolt stretch parameter */
/*
  Copyright (C) 2007 University of Texas at Austin
   
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

float vt2w (int nt      /* time axis length */, 
	    float *vt   /* [nt] velocity in time */,
	    float *str  /* [nt] stretch (output) */)
/*< Evaluate Stolt stretch >*/
{
    int it, it1;
    float v2,vrms2,s2,s, sum1,sum2,sum3,sum4, wi, wbar;

    sum1=sum2=sum3=sum4=0.;
    for (it=0; it < nt; it++) {
	it1 = it+1;
	v2=vt[it]*vt[it];
	sum1 += v2;     /* sum1 = \int v^2 = t v_rms^2 */ 
	vrms2 = sum1/it1;			     
	sum2 += v2*v2;  /* sum2 = \int v^4 */
	sum3 += sum1; 	/* sum3 = \int t v_rms^2 */
	s2=sum3*2.;
      
	str[it] = sqrtf(s2); /* s*vrms */

	if (sum1 > FLT_EPSILON) {
	    s = sum2*it1/(sum1*sum1); /* heterogeneity */
	    wi = 1.-(s2/(sum1*it1))*(v2/vrms2-s);
	} else {
	    s = 1.;
	    wi = 1.;
	}

	sum4 += wi;      /* sum4 = \int w */
    }
    wbar = sum4/nt;
    
    return wbar;
}
