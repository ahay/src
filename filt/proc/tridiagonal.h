#ifndef _tridiagonal_h
#define _tridiagonal_h

typedef struct Tris *tris;

tris tridiagonal_init (int n);
void tridiagonal_define (tris slv, float* diag, float* offd);
void tridiagonal_const_define (tris slv, float diag, float offd);
void tridiagonal_solve (tris slv, float* b);
void tridiagonal_close (tris slv);

#endif

/* 	$Id: tridiagonal.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
