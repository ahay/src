#ifndef _cbanded_h
#define _cbanded_h

void cbanded_init (int n_in, int band_in);
void cbanded_const_define (float complex diag, const float complex *offd);
void cbanded_solve (float complex *b);
void cbanded_close (void);

#endif

/* 	$Id: cbanded.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
