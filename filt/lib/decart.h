#ifndef _sf_decart_h
#define _sf_decart_h

/* index transform (vector to matrix) and its inverse */
void sf_line2cart( int dim, const int* nn, int i, int* ii);

int sf_cart2line( int dim, const int* nn, const int* ii);

#endif
