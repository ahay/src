#ifndef _sf_decart_h
#define _sf_decart_h

/* index transform (vector to matrix) and its inverse */
void sf_line2cart( int dim, const int* nn, int i, int* ii);

int sf_cart2line( int dim, const int* nn, const int* ii);

int sf_first_index (int i, int j, int dim, const int *n, const int *s);

#endif

/* 	$Id: decart.h,v 1.3 2003/09/29 14:34:55 fomels Exp $	 */
