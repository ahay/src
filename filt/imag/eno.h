#ifndef _eno_h
#define _eno_h

/*
  File: eno.h
  -----------
  ENO 1-D interpolation
*/

/* abstract data type */
typedef struct Eno *eno;

/* flag values */
typedef enum {FUNC, DER, BOTH} der;

/* 
   Function: eno_init
   ------------------
   Initialize interpolation object
   order - interpolation order
   n     - data length
*/
eno eno_init (int order, int n);

/* 
   Function: eno_close
   -------------------
   Free internal storage
*/
void eno_close (eno ent);

/* 
   Function: eno_set
   -----------------
   Set the interpolation table
   ent  - ENO object
   c[n] - data
   Note: c can be changed or freed afterwords
*/
void eno_set (eno ent, float* c);

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
void eno_apply (eno ent, int i, float x, float *f, float *f1, der what);

#endif
