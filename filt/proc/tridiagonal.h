#ifndef _tridiagonal_h
#define _tridiagonal_h

typedef struct Tris *tris;

tris tridiagonal_init (int n);
void tridiagonal_define (tris slv, float* diag, float* offd);
void tridiagonal_const_define (tris slv, float diag, float offd);
void tridiagonal_solve (tris slv, float* b);
void tridiagonal_close (tris slv);

#endif
