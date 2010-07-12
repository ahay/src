/* 
 * Copyright (c) 1997 Stanford Exploration Project
 * All Rights Reserved
 *
 * File: sqroot.c 
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
#include "pqueue.h"

#define ABS(x) ((x) < 0 ? -(x) : (x))

static float **x, **xn, **x1;

void heap_init (int n)
{
  x = (float **) malloc ((n+1)*sizeof (float *)); 
  xn = x;
  x1 = x+1;
}
  
void heap_close (void)
{
  free (x);
}

void heap_insert (float* v)
{
  float **xi, **xp;
  unsigned int p;

  xi = ++xn;
  *xi = v;
  p = xn-x;
  for (p >>= 1; p > 0; p >>= 1) {
    xp = x + p;
    if (*v > **xp) break;
    *xi = *xp; xi = xp;
  }
  *xi = v; 
}

float* heap_extract (void)
{
  unsigned int c;
  int n;
  float *v, *t;
  float **xi, **xc;

  v = *x1;
  *(xi = x1) = t = *(xn--);
  n = xn-x;
  for (c = 2; c <= n; c <<= 1) {
    xc = x + c;
    if (c < n && **xc > **(xc+1)) {
      c++; xc++;
    }
    if (*t <= **xc) break;
    *xi = *xc; xi = xc;
  }
  *xi = t;
  return v;
}

void heap_insertl (float* v)
{
  float **xi, **xp;
  unsigned int p;

  xi = ++xn;
  *xi = v;
  p = xn-x;
  for (p >>= 1; p > 0; p >>= 1) {
    xp = x + p;
    if (ABS(*v) > ABS(**xp)) break;
    *xi = *xp; xi = xp;
  }
  *xi = v; 
}

float* heap_extractl (void)
{
  unsigned int c;
  int n;
  float *v, *t;
  float **xi, **xc;

  v = *x1;
  *(xi = x1) = t = *(xn--);
  n = xn-x;
  for (c = 2; c <= n; c <<= 1) {
    xc = x + c;
    if (c < n && ABS(**xc) > ABS(**(xc+1))) {
      c++; xc++;
    }
    if (ABS(*t) <= ABS(**xc)) break;
    *xi = *xc; xi = xc;
  }
  *xi = t;
  return v;
}












