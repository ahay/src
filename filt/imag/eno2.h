#ifndef _eno2_h
#define _eno2_h

/*
  File: eno2.h
  ------------
  ENO 2-D interpolation
*/


#include "eno.h"

/* abstract data type */
typedef struct Eno2 *eno2;

/* 
   Function: eno2_init
   -------------------
   Initialize interpolation object
   order  - interpolation order
   n1, n2 - data dimensions
*/
eno2 eno2_init (int order, int n1, int n2);

/* 
   Function: eno2_set
   ------------------
   Set the interpolation table
   pnt       - ENO object
   c[n2][n1] - data
   Note: c can be changed or freed afterwords
*/
void eno2_set (eno2 pnt, float** c);

/* 
   Function: eno2_set1
   -------------------
   Set the interpolation table
   pnt      - ENO object
   c[n2*n1] - data
   Note: c can be changed or freed afterwords
*/
void eno2_set1 (eno2 pnt, float* c);

/* 
   Function: eno2_close
   --------------------
   Free internal storage
*/
void eno2_close (eno2 pnt);

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
		 float* f, float* f1, der what);

#endif
