#ifndef _sf_alloc_h
#define _sf_alloc_h

#include <stdlib.h>

#include "c99.h"

/*@out@*/ void *sf_alloc (size_t n, size_t size);
/*@out@*/ void *sf_realloc (void* ptr, size_t n, size_t size);
/*@out@*/ char *sf_charalloc (size_t n);
/*@out@*/ unsigned char *sf_ucharalloc (size_t n);
/*@out@*/ int *sf_intalloc (size_t n);
/*@out@*/ float *sf_floatalloc (size_t n);
/*@out@*/ bool *sf_boolalloc (size_t n);
/*@out@*/ bool **sf_boolalloc2 (size_t n1, size_t n2);

#ifndef __cplusplus
/*@out@*/ float complex *sf_complexalloc (size_t n);
/*@out@*/ float complex **sf_complexalloc2 (size_t n1, size_t n2);
#endif

/*@out@*/ float **sf_floatalloc2 (size_t n1, size_t n2);
/*@out@*/ float ***sf_floatalloc3 (size_t n1, size_t n2, size_t n3);
/*@out@*/ int **sf_intalloc2 (size_t n1, size_t n2);
/*@out@*/ char **sf_charalloc2 (size_t n1, size_t n2);
/*@out@*/ unsigned char **sf_ucharalloc2 (size_t n1, size_t n2);

#endif
