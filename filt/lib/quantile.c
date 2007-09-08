/* Computing quantiles by Hoare's algorithm */
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

#include "quantile.h"

float sf_quantile(int q    /* quantile */, 
		  int n    /* array length */, 
		  float* a /* array [n] */) 
/*< find quantile (caution: a is changed) >*/ 
{
    float *i, *j, ak, *low, *hi, buf, *k;

    low=a;
    hi=a+n-1;
    k=a+q; 
    while (low<hi) {
	ak = *k;
	i = low; j = hi;
	do {
	    while (*i < ak) i++;     
	    while (*j > ak) j--;     
	    if (i<=j) {
		buf = *i;
		*i++ = *j;
		*j-- = buf;
	    }
	} while (i<=j);
	if (j<k) low = i; 
	if (k<i) hi = j;
    }
    return (*k);
}
