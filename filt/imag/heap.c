/* 
 * Copyright (c) 1997 Stanford Exploration Project
 * All Rights Reserved
 *
 * File: heap.c 
 */

/*
 * This file provides the implementation for the
 * simple heap priority queue
 *
 * Reference: Sedgewick "Algorithms in C"
 * 
 * Author: Sergey Fomel
 */

#include <stdlib.h>
#include "heap.h"

static Point **x, **xn, **x1;

void heap_init (int n)
{
  x = (Point **) malloc ((n+1)*sizeof (Point *)); 
  xn = x;
  x1 = x+1;
}
  
void heap_close (void)
{
  free (x);
}

void heap_insert (Point* v)
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
{
  unsigned int c;
  int n;
  Point *v, *t;
  Point **xi, **xc;

  v = *x1;
  *(xi = x1) = t = *(xn--);
  n = xn-x;
  for (c = 2; c <= n; c <<= 1) {
    xc = x + c;
    if (c < n && (*xc)->v < (*(xc+1))->v) {
      c++; xc++;
    }
    if (t->v >= (*xc)->v) break;
    (*xc)->h = xi; *xi = *xc; xi = xc;
  }
  t->h = xi; *xi = t; 
  return v;
}

void heap_update (Point *v)
{
  unsigned int c;
  int n;
  Point **xc, **xi;

  xi = v->h; *xi = v; 
  n = xn-x; c = xi-x;
  for (c <<= 1; c <= n; c <<= 1) {
    xc = x + c;
    if (c < n && (*xc)->v < (*(xc+1))->v) {
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













