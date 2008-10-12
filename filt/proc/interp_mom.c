/* B-spline interpolation. */
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
#include <rsf.h>

#include "interp_mom.h"

static void mom8_int (float x, float* w);
static void mom6_int (float x, float* w);
static void mom4_int (float x, float* w);

void mom_int (float x, int n, float* w)
/*< interpolation function >*/
{
    switch (n) {    
	case 8: 
	    mom8_int (x, w); 
	    break;
	case 6: 
	    mom6_int (x, w); 
	    break;
	case 4: 
	    mom4_int (x, w); 
	    break;
	default: 
	    sf_error ("%s: mom_int length %d not implemented",
		      __FILE__,n); 
	    break;
    }
}

static void mom4_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = (3. + x*((7. - 2.*x)*x - 8.))/16.;
    w[1] = (10. + x2*(6.*x - 13.))/16.;
    w[2] = (3.  + x*(8. + (5. - 6.*x)*x))/16.;
    w[3] = x2*(1. + 2.*x)/16.;
}

static void mom6_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = (3. + x*(x*(26. + x*((11. - 2.*x)*x - 24.)) - 14.))/288.;
    w[1] = (64. + x*(x*(40. + x*(48. + x*(10.*x - 43.))) - 116.))/288.;
    w[2] = (77. + x2*((31. - 10.*x)*x2 - 66.))/144.;
    w[3] = (32. + x*(58. + x*(20. + x*(x*(10.*x - 19.) - 24.))))/144.;
    w[4] = (3. + x*(14. + x*(26. + x*(24. + (7. - 10.*x)*x))))/288.;
    w[5] = x2*x2*(1. + 2.*x)/288.;
}

static void mom8_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = (3 + x*(x*(57. + x*(x*(85. + x*((15. - 2.*x)*x - 48.)) - 90.)) 
		   - 20.))/11520.;
    w[1] = (298. + x*(x*(1158. + x*(x*(x*(192. + x*(14.*x - 89.)) - 30.) 
				    - 600.)) - 940.))/11520.;
    w[2] = (2741. + x*(3.*x*(205. + x*(490. + x*(x*((73. - 14.*x)*x - 80.) 
						 - 215.))) - 3820.))/11520.;
    w[3] = (5436. + 5.*x2*(x2*(236. + x2*(14.*x - 57.)) - 732.))/11520.;
    w[4] = (2741. + 5.*x*(764. + x*(123. + x*(x*(x*(48. + (41. - 14.*x)*x) 
						 - 129.) - 294.))))/11520.;
    w[5] = (298. + x*(940. + 
		      3.*x*(386. + x*(200. + x*(x*(x*(14.*x - 25.) - 64.) 
						- 10.)))))/11520.;
    w[6] = (3. + x*(20. + x*(57. + 
			     x*(90. + x*(85. + x*(48. + (9. - 
							 14.*x)*x))))))/11520.;
    w[7] = x2*x2*x2*(1. + 2.*x)/11520.;
}

/* 	$Id: interp_mom.c 999 2005-02-12 14:00:52Z fomels $	 */
