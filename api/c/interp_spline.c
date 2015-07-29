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
#include "interp_spline.h"
#include "error.h"

static void spline8_int (float x, float* w);
static void spline6_int (float x, float* w);
static void spline3_int (float x, float* w);
static void spline2_int (float x, float* w);

static void spline8_der (float x, float* w);
static void spline6_der (float x, float* w);
static void spline3_der (float x, float* w);
static void spline2_der (float x, float* w);

void sf_spline_int (float x, int n, float* w)
/*< interpolation function >*/
{
    switch (n) {    
	case 8: 
	    spline8_int (x, w); 
	    break;
	case 6: 
	    spline6_int (x, w); 
	    break;
	case 4: 
	    sf_spline4_int (x, w); 
	    break;
	case 3: 
	    spline3_int (x, w); 
	    break;
	case 2: 
	    spline2_int (x, w); 
	    break;
	default: 
	    sf_error ("%s: spline_int length %d not implemented",
		      __FILE__,n); 
	    break;
    }
}

void sf_spline_der (float x, int n, float* w)
/*< derivative computation >*/
{
    switch (n) {    
	case 8: 
	    spline8_der (x, w); 
	    break;
	case 6: 
	    spline6_der (x, w); 
	    break;
	case 4: 
	    sf_spline4_der (x, w); 
	    break;
	case 3: 
	    spline3_der (x, w); 
	    break;
	case 2: 
	    spline2_der (x, w); 
	    break;
	default: 
	    sf_error ("%s: spline_der length %d not implemented",
		      __FILE__,n); 
	    break;
    }
}

static void spline2_int (float x, float* w)
{
    w[0] = 1. - x; 
    w[1] = x; 
}

static void spline2_der (float x, float* w)
{
    w[0] = -1.; 
    w[1] = 1.; 
}

static void spline3_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = 0.5*(x2-x)+0.125;
    w[1] = 0.75-x2;
    w[2] = 0.5*(x2+x)+0.125;
}

static void spline3_der (float x, float* w)
{
    w[0] = x-0.5;
    w[1] = -2.*x;
    w[2] = x+0.5;
}

void sf_spline4_int (float x, float* w)
/*< Cubic spline interpolation >*/
{
    float x2;
    
    x2 = x*x;
    w[0] = (1. + x*((3. - x)*x-3.))/6.;
    w[1] = (4. + 3.*(x -2.)*x2)/6.;
    w[2] = (1. + 3.*x*(1. + (1. - x)*x))/6.;
    w[3] = x2*x/6.;
}

void sf_spline4_der (float x, float* w)
/*< Cubic spline derivative >*/
{
    float x1;
    
    x1 = 1.-x;
    w[0] = -0.5*x1*x1;
    w[1] = 0.5*x*(3.*x-4.);
    w[2] = 0.5*x1*(3.*x+1.);
    w[3] = 0.5*x*x;
}

static void spline6_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = (1. + x*(x*(10. + x*((5. - x)*x-10.))-5.))/120.;
    w[1] = (26. + 5.*x*(x*(4. + x*(4. + x*(x-4.)))-10.))/120.;
    w[2] = (66. + 10.*x2*((3. - x)*x2-6.))/120.;
    w[3] = (26. + 10.*x*(5. + x*(2. + x*(x*(x-2.)-2.))))/120.;
    w[4] = (1. + 5.*x*(1. + x*(2. + x*(2. + (1. - x)*x))))/120.;
    w[5] = x2*x2*x/120.;
}

static void spline6_der (float x, float* w)
{
    float x1, x2;
    
    x1 = 1.-x;
    x2 = x*x;
    
    w[0] = -x1*x1*x1*x1/24.;
    w[1] = (x*(8. + x*(12. + x*(5.*x - 16.))) - 10.)/24.;
    w[2] = (x2*(1. - 5.*x/12.) - 1.)*x;
    w[3] = (5. + x*(4. + x*(x*(5.*x - 8.) - 6.)))/12.;
    w[4] = (1. + x*(4. + x*(6. + x*(4. - 5.*x))))/24.;
    w[5] = x2*x2/24.;
}

static void spline8_int (float x, float* w)
{
    float x2;
    
    x2 = x*x;
    w[0] = (1. + x*(x*(21. + x*(x*(35. + x*((7. - x)*x-21.))-35.))-7.))/5040.;
    w[1] = (120. + 7.*x*(x*(72. + x*(x2*(12. + x*(x-6.))-40.))-56.))/5040.;
    w[2] = (1191. + 7.*x*(x*(45. + x*(95. + 3.*x*(x*((5. - x)*x-5.)
						  -15.)))-245.))/5040.;
    w[3] = (2416. + 35.*x2*(x2*(16. + x2*(x-4.))-48.))/5040.;
    w[4] = (1191. + 35.*x*(49. + x*(9. + x*(x*(x*(3. + (3. - x)*x)-9.)
					    -19.))))/5040.;
    w[5] = (120. + 7.*x*(56. + x*(72. + x*(40. + 3.*x2*(x*(x-2.)-4.)))))/5040.;
    w[6] = (1. + 7.*x*(1. + x*(3. + x*(5. + x*(5. + x*(3. + 
						       (1. - x)*x))))))/5040.;
    w[7] = x2*x2*x2*x/5040.;
}

static void spline8_der (float x, float* w)
{
    float x1, x2;
    
    x1 = 1.-x;
    x2 = x*x;
    
    w[0] = -x1*x1*x1*x1*x1*x1/720.;
    w[1] = (x*(144. + x*(x2*(60. + x*(7.*x - 36.)) - 120.)) - 56.)/720.;
    w[2] = (3.*x*(30. + x*(95. + x*(x*(x*(30. - 7.*x)- 25.)- 60.)))- 245.)/720.;
    w[3] = x*(x2*(64. + x2*(7.*x-24.)) - 96.)/144.;
    w[4] = (49. + x*(18. + x*(x*(x*(15. + x*(18. - 7.*x)) - 36.) - 57.)))/144.;
    w[5] = (56. + 3.*x*(48. + x*(40. + x2*(x*(7.*x-12.) - 20.))))/720.;
    w[6] = (1. + x*(6. + x*(15. + x*(20. + x*(15. + x*(6. - 7.*x))))))/720.;
    w[7] = x2*x2*x2/720.;
}

/* 	$Id: interp_spline.c 7107 2011-04-10 02:04:14Z ivlad $	 */
