#include <rsf.h>

#include "eno2.h"
#include "eno.h"

/* concrete data type */
struct Eno2 {
    int order, ng, n1, n2;
    eno jnt, *ent;
    float *f, *f1;
};

/* 
   Function: eno2_init
   -------------------
   Initialize interpolation object
   order  - interpolation order
   n1, n2 - data dimensions
*/
eno2 eno2_init (int order, int n1, int n2)
{
    eno2 pnt;
    int i2;
    
    pnt = (eno2) sf_alloc(1,sizeof(*pnt));
    pnt->order = order; 
    pnt->n1 = n1; 
    pnt->n2 = n2;
    pnt->ng = 2*order-2;
    if (pnt->ng > pnt->n2) sf_error("%s: ng=%d is too big",__FILE__,pnt->ng);
    pnt->jnt = eno_init (order, pnt->ng);
    pnt->f  = sf_floatalloc(pnt->ng);
    pnt->f1 = sf_floatalloc(pnt->ng);
    pnt->ent = (eno*) sf_alloc(n2,sizeof(eno));
    for (i2 = 0; i2 < n2; i2++) {
	pnt->ent[i2] = eno_init (order, n1);
    }

    return pnt;
}

/* 
   Function: eno2_set
   ------------------
   Set the interpolation table
   pnt       - ENO object
   c[n2][n1] - data
   Note: c can be changed or freed afterwords
*/
void eno2_set (eno2 pnt, float** c)
{
    int i2;
    
    for (i2 = 0; i2 < pnt->n2; i2++) {
	eno_set (pnt->ent[i2], c[i2]);
    }
}

/* 
   Function: eno2_set1
   -------------------
   Set the interpolation table
   pnt      - ENO object
   c[n2*n1] - data
   Note: c can be changed or freed afterwords
*/
void eno2_set1 (eno2 pnt, float* c)
{
    int i2;
    
    for (i2 = 0; i2 < pnt->n2; i2++) {
	eno_set (pnt->ent[i2], c+i2*(pnt->n1));
    }
}

/* 
   Function: eno2_close
   --------------------
   Free internal storage
*/
void eno2_close (eno2 pnt)
{
    int i2;
    
    eno_close (pnt->jnt);
    for (i2 = 0; i2 < pnt->n2; i2++) {
	eno_close (pnt->ent[i2]);
    }
    free (pnt->f);
    free (pnt->f1);
    free (pnt->ent);
    free (pnt);
}

/* 
   Function: eno2_apply
   --------------------
   Apply interpolation
   pnt   - ENO object
   i,j   - grid location
   x,y   - offsets from grid
   f     - data value (output)
   f1[2] - derivative value (output)
   what  - flag of what to compute: 
            FUNC - function value
	    DER  - derivative value
	    BOTH - both
*/
void eno2_apply (eno2 pnt, int i, int j, float x, float y, 
		 float* f, float* f1, der what)
{
    int k, b2;
    float g;
    
    if (j-pnt->order+2 < 0) {
	b2 = 0;
    } else if (j+pnt->order-1 > pnt->n2-1) {
	b2 = pnt->n2 - pnt->ng;
    } else {
	b2 = j-pnt->order+2;
    }
    
    j -= b2;
    
    for (k = 0; k < pnt->ng; k++) {
	if (what != FUNC) {
	    eno_apply (pnt->ent[b2+k],i,x,pnt->f+k, pnt->f1+k,BOTH);
	} else {
	    eno_apply (pnt->ent[b2+k],i,x,pnt->f+k, pnt->f1+k,FUNC);
	}
    }
    
    eno_set (pnt->jnt,pnt->f);
    eno_apply (pnt->jnt,j,y,f,f1+1,what);
    
    if (what != FUNC) {
	eno_set (pnt->jnt,pnt->f1);
	eno_apply(pnt->jnt,j,y,f1,&g,FUNC);
    }
}

/* 	$Id: eno2.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
