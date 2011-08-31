/* Priority queue (OMP) */
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "pqueueomp.h"

static float ***x, ***xn, ***x1;

void pqueue_init (int n)
/*< Initialize heap with the maximum size >*/
{
    int its, mts;

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    for (its=0; its < mts; its++)
	x[its] = (float **) sf_alloc ((n+1),sizeof (float *)); 
}

void pqueue_start (void)
/*< Set starting values >*/
{
    int its;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    xn[its] = x[its];
    x1[its] = x[its]+1;
}

void pqueue_close (void)
/*< Free the allocated storage >*/
{
    free (x);
}

void pqueue_insert (float* v)
/*< Insert an element (smallest first) >*/
{
    int its;
    float **xi, **xq;
    unsigned int q;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    xi = ++xn[its];
    *xi = v;
    q = (unsigned int) (xn[its]-x[its]);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x[its] + q;
	if (*v > **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

void pqueue_insert2 (float* v)
/*< Insert an element (largest first) >*/
{
    int its;
    float **xi, **xq;
    unsigned int q;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    xi = ++xn[its];
    *xi = v;
    q = (unsigned int) (xn[its]-x[its]);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x[its] + q;
	if (*v < **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

float* pqueue_extract (void)
/*< Extract the smallest element >*/
{
    int its;
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;
    
#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    v = *(x1[its]);
    *(xi = x1[its]) = t = *(xn[its]--);
    n = (int) (xn[its]-x[its]);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) n && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*t <= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

float* pqueue_extract2 (void)
/*< Extract the largest element >*/
{
    int its;
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    v = *(x1[its]);
    *(xi = x1[its]) = t = *(xn[its]--);
    n = (int) (xn[its]-x[its]);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) n && **xc < **(xc+1)) {
	    c++; xc++;
	}
	if (*t >= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

void pqueue_update (float **v)
/*< restore the heap: the value has been altered >*/
{
    int its;
    unsigned int c;
    int n;
    float **xc, **xi;

#ifdef _OPENMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif
    
    xi = v; 
    n = (int) (xn[its]-x[its]); c = (unsigned int) (xi-x[its]);
    for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
	xc = x[its] + c;
	if (c < (unsigned int) n && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (**v <= **xc) break;
	*xi = *xc; xi = xc;
    }
    xi = v; c = (unsigned int) (xi-x[its]);
    for (c >>= 1; c > 0; c >>= 1) {
	xc = x[its] + c;
	if (**v > **xc) break;
	*xi = *xc; xi = xc; 
    }
    *xi = *v; 
}

/* 	$Id: pqueue.c 7537 2011-07-28 06:13:30Z sfomel $	 */
