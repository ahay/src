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

/* 	$Id: ctridiagonal.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
