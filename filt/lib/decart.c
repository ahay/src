#include "decart.h"

/* index transform (vector to matrix) and its inverse */
void sf_line2cart( int dim, const int* nn, int i, int* ii) {
    int axis;
 
    for (axis = 0; axis < dim; axis++) {
	ii[axis] = i%nn[axis];
	i /= nn[axis];
    }
}

int sf_cart2line( int dim, const int* nn, const int* ii) {
    int i, axis;

    if (dim < 1) return 0;

    i = ii[dim-1];
    for (axis = dim-2; axis >= 0; axis--) {
	i = i*nn[axis] + ii[axis];
    }
    return i;
}

int sf_first_index (int i, int j, int dim, const int *n, const int *s)
{
    int i0, n123, k, ii;

    n123 = 1;
    i0 = 0;
    for (k=0; k < dim; k++) {
	if (k == i) continue;
	ii = (j/n123)%n[k]; /* to cartesian */
	n123 *= n[k];	
	i0 += ii*s[k];      /* back to line */
    }

    return i0;
}
