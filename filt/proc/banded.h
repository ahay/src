#ifndef _banded_h
#define _banded_h

typedef struct Bands *bands;

bands banded_init (int n, int band);
void banded_define (bands slv, float* diag, float** offd);
void banded_const_define (bands slv, float diag, const float* offd);
void banded_solve (bands slv, float* b);
void banded_close (bands slv);

#endif

/* 	$Id: banded.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
