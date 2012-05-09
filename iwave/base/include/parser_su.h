/* 
parser_su.h
Igor Terentyev.
********************************************************************************
Wrapper to provide SU GETPARS functions.

Notes:
1) ps_destroyargs() function is not present in the SU package.
   I beleive SU has a memory leak, because allocated memory is not freed.
   My ps_initargs(..) function, if called subsequently, 
   calls ps_destroyargs() function and then allocates new memory.
   So ps_destroyargs() is not needed between subsequent ps_initargs(..) calls.
   Call it only in the end of the program or to free the memory 
   (if ps_getpar.. functions will not be used anymore).

2) ps_countparval(..) counts commas (like SU). But there may be empty values
   between commas, etc. Counting functions from parser package behave
   differently - they count number of readable values (which depends on type).
   
3) Extra functions to work with ireal data type: ps_getparreal, ps_getnparreal.

4) Unlike SU, ps_getparstring allocates memory and returns a copy of the string.
   SU is inconsistent here, because its getparstringarray allocates and makes
   a copy and getparstring does not.
*/
/*============================================================================*/

#ifndef __PARSER_SU_H_
#define __PARSER_SU_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
/*----------------------------------------------------------------------------*/

void ps_destroyargs();                            /* Not in the SU package */

void ps_initargs(int argc, char **argv);

int ps_getparstring(char *name, char **p);        /* ALLOCATES UNLIKE SU */
int ps_getparint(char *name, int *p);
int ps_getparuint(char *name, unsigned int *p);
int ps_getparshort(char *name, short *p);
int ps_getparushort(char *name, unsigned short *p);
int ps_getparlong(char *name, long *p);
int ps_getparulong(char *name, unsigned long *p);
int ps_getparfloat(char *name, float *p);
int ps_getpardouble(char *name, double *p);
int ps_getparreal(char *name, ireal *p);           /* Not in the SU package */

int ps_getnparstring(int n, char *name, char **p);/* ALLOCATES UNLIKE SU */
int ps_getnparint(int n, char *name, int *p);
int ps_getnparuint(int n, char *name, unsigned int *p);
int ps_getnparshort(int n, char *name, short *p);
int ps_getnparushort(int n, char *name, unsigned short *p);
int ps_getnparlong(int n, char *name, long *p);
int ps_getnparulong(int n, char *name, unsigned long *p);
int ps_getnparfloat(int n, char *name, float *p);
int ps_getnpardouble(int n, char *name, double *p);
int ps_getnparreal(int n, char *name, ireal *p);   /* Not in the SU package */

int ps_getnpar(int n, char *name, char *type, void *p);

int ps_countparname(char *name);
int ps_countparval(char *name);
int ps_countnparval(int n, char *name);

/* For ProMAX */
void ps_getPar(char *name, char *type, void *p);
/*----------------------------------------------------------------------------*/
/* 
NOT IMPLEMENTED.
*/
/*int ps_getparstringarray (char *name, char **p);
int ps_getnparstringarray (int n, char *name, char **p);*/
/*----------------------------------------------------------------------------*/
/* 
FOR TESTS.
*/
void* ps_getpararray();                           /* Not in the SU package */
/*----------------------------------------------------------------------------*/

#endif /*__PARSER_SU_H_*/
