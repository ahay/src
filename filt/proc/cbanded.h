#ifndef _cbanded_h
#define _cbanded_h

#include <rsf.h>

/* Inverting Hermitian positive definite matrix.

In define: offd are lower subdiagonals.
*/

void cbanded_init (int n_in, int band_in);

void cbanded_const_define (float diag, const float complex *offd);
void cbanded_define (const float *diag, float complex **offd);
void cbanded_solve (float complex *b);
void cbanded_close (void);

#endif

/* 	$Id: cbanded.h,v 1.1 2004/03/13 06:11:12 fomels Exp $	 */
