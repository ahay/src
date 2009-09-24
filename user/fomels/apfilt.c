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
    int j, k, n;
    double ak;
    
    n = nw*2;
   
    for (k=0; k <= n; k++) {
	ak = 1.0;
	for (j=0; j < n-k; j++) {
	    ak *= (k+j+1)*(n-j-p)/(2*(2*j+1)*(j+1));
	}
	for (; j < n; j++) {
	    ak *= (p+j+1.0)/(2*(2*j+1));
	}
	a[k] = ak;
    }
}

void aderfilter (int nw   /* size */, 
		 float p  /* slope */, 
		 float* a /* output filter [2*nw+1] */)
/*< find coefficients for filter derivative >*/
{

    int i, j, k, n;
    double ak, aj;
    
    n = nw*2;
   
    for (k=0; k <= n; k++) {
	a[k] = 0.;
	for (i=0; i < n; i++) {
	    ak = -1.0;
	    for (j=0; j < n-k; j++) {
		aj = (k+j+1.0)/(2*(2*j+1)*(j+1));
		if (j != i) {
		    ak *= aj*(n-j-p);
		} else {
		    ak *= -aj;
		}
	    }
	    for (; j < n; j++) {
		aj = 1.0/(2*(2*j+1));
		if (j != i) {
		    ak *= aj*(p+j+1.0);
		} else {
		    ak *= aj;
		}
	    }
	    a[k] += ak;
	}
    }
}

/* 	$Id$	 */
