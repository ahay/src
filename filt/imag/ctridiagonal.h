#ifndef _ctridiagonal_h
#define _ctridiagonal_h

#include <rsf.h>

typedef struct CTris *ctris;

ctris ctridiagonal_init (int n);
void ctridiagonal_define (ctris slv, 
			  float complex* diag, float complex* offd);
void ctridiagonal_const_define (ctris slv, 
				float complex diag, float complex offd);
void ctridiagonal_solve (ctris slv, float complex* b);
void ctridiagonal_close (ctris slv);

#endif
