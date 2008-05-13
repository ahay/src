/* All-pass plane-wave destruction filter coefficients */
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

#include "apfilt.h"

void passfilter (int nw   /* size */, 
		 float p  /* slope */, 
		 float* a /* output filter [2*nw+1] */)
/*< find filter coefficients >*/
{
    switch (nw) {
	case 1:   	
	    a[0] = (1.-p)*(2.-p)/12.;
	    a[1] = (2.+p)*(2.-p)/6.;	
	    a[2] = (1.+p)*(2.+p)/12.;
	    break;
	case 2:
	    a[0] = (1.-p)*(2.-p)*(3.-p)*(4.-p)/1680.;
	    a[1] = (4.-p)*(2.-p)*(3.-p)*(4.+p)/420.;
	    a[2] = (4.-p)*(3.-p)*(3.+p)*(4.+p)/280.;
	    a[3] = (4.-p)*(2.+p)*(3.+p)*(4.+p)/420.;
	    a[4] = (1.+p)*(2.+p)*(3.+p)*(4.+p)/1680.;
	    break;
	case 3:
	    a[0] = (1.-p)*(2.-p)*(3.-p)*(4.-p)*(5.-p)*(6.-p)/665280.;
	    a[1] = (6.-p)*(2.-p)*(3.-p)*(4.-p)*(5.-p)*(6.+p)/110880.;
	    a[2] = (6.-p)*(5.-p)*(3.-p)*(4.-p)*(5.+p)*(6.+p)/44352.;
	    a[3] = (6.-p)*(5.-p)*(4.-p)*(4.+p)*(5.+p)*(6.+p)/33264.;
	    a[4] = (6.-p)*(5.-p)*(3.+p)*(4.+p)*(5.+p)*(6.+p)/44352.;
	    a[5] = (6.-p)*(2.+p)*(3.+p)*(4.+p)*(5.+p)*(6.+p)/110880.;
	    a[6] = (1.+p)*(2.+p)*(3.+p)*(4.+p)*(5.+p)*(6.+p)/665280.;
	    break;
	default:
	    sf_error("%s: filter size %d is not implemented",__FILE__,nw);
	    break;
    }
}

void aderfilter (int nw   /* size */, 
		 float p  /* slope */, 
		 float* a /* output filter [2*nw+1] */)
/*< find coefficients for filter derivative >*/
{
    switch (nw) {
	case 1: 
	    a[0] = (3.-2.*p)/12.;
	    a[1] = p/3.;
	    a[2] = (-3.-2.*p)/12.;
	    break;
	case 2:
	    a[0] = (5.-2.*p)*(5.+p*(p-5.))/840.;	
	    a[1] = (80.+p*(p*(4.*p-15.)-20.))/420.;
	    a[2] = p*(25.-2.*p*p)/140.;
	    a[3] = (p*(p*(4.*p+15.)-20.)-80.)/420.;
	    a[4] = (-5.-2.*p)*(5.+p*(p+5.))/840.;
	    break;
	case 3:
	    a[0] = (1764. + 
		    p*(-3248. + p*(2205. + p*(-700. + 
					      (105. - 6.*p)*p))))/665280.;
	    a[1] = (2772. + p*(-2436. + p*(525. + 
					   p*(70. + p*(-35. + 3.*p)))))/55440.;
	    a[2] = (6300. + p*(-336. + p*(-1281. + 
					  p*(196. + (35. - 6.*p)*p))))/44352.;
	    a[3] = (p*(1876. + p*p*(-154. + 3.*p*p)))/16632.;
	    a[4] = (-6300. + p*(-336. + p*(1281. + 
					   p*(196. + (-35. - 6.*p)*p))))/44352.;
	    a[5] = (-2772. + p*(-2436. + p*(-525. + 
					    p*(70. + p*(35. + 3.*p)))))/55440.;
	    a[6] = (-1764. + p*(-3248. + p*(-2205. + 
					    p*(-700. + 
					       (-105. - 6.*p)*p))))/665280.;
	    break;
	default:
	    sf_error("%s: filter size %d not implemented",__FILE__,nw);
	    break;
    }
}

/* 	$Id$	 */
