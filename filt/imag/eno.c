#include <math.h>
#include <float.h>

#include <rsf.h>

#include "eno.h"

/* concrete data type */
struct Eno {
    int order, n;
    float** diff;
};

static const float big_number = FLT_MAX;

#ifndef MAX
#define MAX(a,b) ((a)>(b))?(a):(b)
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

/* 
   Function: eno_init
   ------------------
   Initialize interpolation object
   order - interpolation order
   n     - data length
*/
eno eno_init (int order, int n)
{
    eno ent;
    int i;
    
    ent = (eno) sf_alloc(1,sizeof(*ent));
    ent->order = order;
    ent->n = n;
    ent->diff = (float**) sf_alloc(order,sizeof(float*));
    for (i = 0; i < order; i++) {
	ent->diff[i] = sf_floatalloc(n-i);
    }
  
    return ent;
}

/* 
   Function: eno_close
   -------------------
   Free internal storage
*/
void eno_close (eno ent)
{
    int i;

    for (i = 0; i < ent->order; i++) {
	free(ent->diff[i]);
    }
    free (ent->diff);
    free (ent);
}

/* 
   Function: eno_set
   -----------------
   Set the interpolation table
   ent  - ENO object
   c[n] - data
   Note: c can be changed or freed afterwords
*/
void eno_set (eno ent, float* c)
{
    int i, j;
    
    for (i=0; i < ent->n; i++) {
	/* copy the initial data */
	ent->diff[0][i] = c[i];
    }
    
    for (j=1; j < ent->order; j++) {
	for (i=0; i < ent->n-j; i++) {
	    /* compute difference tables */
	    ent->diff[j][i] = ent->diff[j-1][i+1] - ent->diff[j-1][i];
	}
    }
}

/* 
   Function: eno_apply
   -------------------
   Apply interpolation
   ent  - ENO object
   i    - grid location
   x    - offset from grid
   f    - data value (output)
   f1   - derivative value (output)
   what - flag of what to compute: 
   FUNC - function value
   DER  - derivative value
   BOTH - both
*/
void eno_apply (eno ent, int i, float x, float *f, float *f1, der what) 
{
    int j, k, i1, i2, n;
    float s, s1, y, w, g, g1;
    
    i2 = MAX(0,MIN(i,ent->n-ent->order));
    i1 = MIN(i2,MAX(0,i-ent->order+2));
    
    w = fabsf(ent->diff[ent->order-1][i1]);
    for (j=i1+1; j <=i2; j++) {
	g = fabsf(ent->diff[ent->order-1][j]);
	if (w > g) w = g;
    }
    
    /* loop over starting points */
    for (g = 0., g1 = 0., n = 0, j = i1; j <= i2; j++) {
	if (fabsf(ent->diff[ent->order-1][j]) > w) continue;
	n++;
        
	y = x + i - j;
	
	/* loop to compute the polynomial */
	for (s = 1., s1 = 0., k=0; k < ent->order; k++) {
	    if (what != FUNC) {
		g1 += s1*ent->diff[k][j];
		s1 = (s + s1*(y-k))/(k+1.);
	    }
	    if (what != DER) g += s*ent->diff[k][j];
	    s *= (y-k)/(k+1.);
	}
    }
    
    if (what != DER) *f = g/n;
    if (what != FUNC) *f1 = g1/n;
}

/* 	$Id: eno.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */

