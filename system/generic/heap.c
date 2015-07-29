/* Another version of priority queue. */
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

#include <stdlib.h>
#include "heap.h"

#ifndef _heap_h

typedef struct Point {
  double v;
  void* d;
  struct Point **h;
} Point; 
/*^*/

#endif

static Point **x, **xn, **x1;

void heap_init (int n)
/*< initialize with the maximum number of points >*/
{
  x = (Point **) malloc ((n+1)*sizeof (Point *)); 
  xn = x;
  x1 = x+1;
}

void heap_close (void)
/*< free the storage >*/
{
  free (x);
}

void heap_insert (Point* v)
/*< insert a point >*/
{
  Point **xi, **xp;
  unsigned int p;

  xi = ++xn;
  *xi = v;
  p = xn-x;
  for (p >>= 1; p > 0; p >>= 1) {
    xp = x + p;
    if (v->v < (*xp)->v) break;
    (*xp)->h = xi; *xi = *xp; xi = xp; 
  }
  v->h = xi; *xi = v; 
}

Point* heap_extract (void)
/*< extract the maximum >*/
{
  unsigned int c;
  int n;
  Point *v, *t;
  Point **xi, **xc;

  v = *x1;
  *(xi = x1) = t = *(xn--);
  n = xn-x;
  for (c = 2; c <= (unsigned int) n; c <<= 1) {
    xc = x + c;
    if (c < (unsigned int) n && (*xc)->v < (*(xc+1))->v) {
      c++; xc++;
    }
    if (t->v >= (*xc)->v) break;
    (*xc)->h = xi; *xi = *xc; xi = xc;
  }
  t->h = xi; *xi = t; 
  return v;
}

void heap_update (Point *v)
/*< restore the heap: the value has been altered >*/
{
  unsigned int c;
  int n;
  Point **xc, **xi;

  xi = v->h; *xi = v; 
  n = xn-x; c = xi-x;
  for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
    xc = x + c;
    if (c < (unsigned int) n && (*xc)->v < (*(xc+1))->v) {
      c++; xc++;
    }
    if (v->v >= (*xc)->v) break;
    (*xc)->h = xi; *xi = *xc; xi = xc;
  }
  v->h = xi; *xi = v; c = xi-x;
  for (c >>= 1; c > 0; c >>= 1) {
    xc = x + c;
    if (v->v < (*xc)->v) break;
    (*xc)->h = xi; *xi = *xc; xi = xc; 
  }
  v->h = xi; *xi = v; 
}
