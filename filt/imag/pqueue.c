#include <rsf.h>

#include "pqueue.h"

static float **x, **xn, **x1;

/* Initialize heap with the maximum size */
void pqueue_init (int n)
{
    x = (float **) sf_alloc ((n+1),sizeof (float *)); 
}

/* Set starting values */
void pqueue_start (void)
{
    xn = x;
    x1 = x+1;
}
  
/* Free the allocated storage */
void pqueue_close (void)
{
    free (x);
}

/* Insert an element */
void pqueue_insert (float* v)
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

/* Extract the smallest element */
float* pqueue_extract (void)
{
  unsigned int c;
  int n;
  float *v, *t;
  float **xi, **xc;

  v = *(x1);
  *(xi = x1) = t = *(xn--);
  n = (int) (xn-x);
  if (n < 0) return NULL;
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

/* 	$Id: pqueue.c,v 1.2 2003/09/30 14:30:53 fomels Exp $	 */

