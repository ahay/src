/* Priority queue (heap sorting) */
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
#include "alloc.h"
#include "pqueue.h"

#ifndef _sf_pqueue_h

enum {SF_IN, SF_FRONT, SF_OUT};
/*^*/

#endif

static float **x, **xn, **x1;

void sf_pqueue_init (int n)
/*< Initialize heap with the maximum size >*/
{
    x = (float **) sf_alloc ((n+1),sizeof (float *)); 
}

void sf_pqueue_start (void)
/*< Set starting values >*/
{
    xn = x;
    x1 = x+1;
}

void sf_pqueue_close (void)
/*< Free the allocated storage >*/
{
    free (x);
}

void sf_pqueue_insert (float* v)
/*< Insert an element (smallest first) >*/
{
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v > **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

void sf_pqueue_insert2 (float* v)
/*< Insert an element (largest first) >*/
{
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v < **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

float* sf_pqueue_extract (void)
/*< Extract the smallest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;
    
    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) n && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*t <= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

float* sf_pqueue_extract2 (void)
/*< Extract the largest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;
    
    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) n && **xc < **(xc+1)) {
	    c++; xc++;
	}
	if (*t >= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

void sf_pqueue_update (float **v)
/*< restore the heap: the value has been altered >*/
{
  unsigned int c;
  int n;
  float **xc, **xi;

  xi = v; 
  n = (int) (xn-x); c = (unsigned int) (xi-x);
  for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
      xc = x + c;
      if (c < (unsigned int) n && **xc > **(xc+1)) {
	  c++; xc++;
      }
      if (**v <= **xc) break;
      *xi = *xc; xi = xc;
  }
  xi = v; c = (unsigned int) (xi-x);
  for (c >>= 1; c > 0; c >>= 1) {
      xc = x + c;
      if (**v > **xc) break;
      *xi = *xc; xi = xc; 
  }
  *xi = *v; 
}

/* 	$Id$	 */

