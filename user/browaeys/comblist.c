/* Subroutines for combinations of k elements out of n */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include <float.h>
#include <rsf.h>

#include "comblist.h"

int binomial (int n, int k)
/*< Computes the binomial coefficient C(n,k) = n! / k!*(n-k)! >*/
/* This code is distributed under the GNU LGPL license. Modified from original author John Burkardt */
{
    int i, icnk, mn, mx;

    mn = SF_MIN(k,(n-k));

    if (mn < 0) {

	icnk = 0;

    } else if (mn == 0) {

	icnk = 1;

    } else {

	mx = SF_MAX(k,(n-k));
	icnk = mx + 1;

	for (i = 2; i < (mn+1); i++) {
	    icnk = (icnk*(mx + i))/i;
	}
    }

    return(icnk);
}

void comb_next (int n, int k, int *a, int *done)
/*< Computes combinations of k elements out of n, computed one at a time, in lexicographical order. >*/
/* Set done to FALSE before the first call.
   Use output value from previous call on subsequent calls.
   Output value will be TRUE as long as there are more
   combinations to compute, and FALSE when the list is exhausted. */
/* This code is distributed under the GNU LGPL license. Modified from original author John Burkardt */
{
    int i;
    int j;

    if (*done) {

	for (i = 0; i < k; i++) a[i] = i+1;

	if ( 1 < k ) {
	    *done = 0;
	} else {
	    *done = 1;
	}

    } else {

	if (a[k-1] < n) {
	    a[k-1] += 1;
	    return;
	}

	for (i = k; i > 1; i--) {
	    if (a[i-2] < (n-k+i-1)) {
		a[i-2] += 1;
		for (j = i; j < (k+1); j++) a[j-1] = a[i-2] + j - (i-1);
		return;
	    }
	}

	*done = 1;
    }

    return;
}
