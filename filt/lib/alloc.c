#include <stdlib.h>

#include "alloc.h"
#include "error.h"
#include "c99.h"

/*@out@*/ void *sf_alloc (size_t n, size_t size)
{
    void *ptr; 
  
    ptr = malloc (n*size);
    if (NULL == ptr)
	sf_error ("%s: cannot allocate %d bytes:", __FILE__, n*size);
    return ptr;
}

void *sf_realloc (void* ptr, size_t n, size_t size)
{
    void *ptr2; 
  
    ptr2 = realloc (ptr,n*size);
    if (NULL == ptr2)
	sf_error ("%s: cannot reallocate %d bytes:", __FILE__, n*size);
    return ptr2;
}

/*@out@*/ char *sf_charalloc (size_t n) 
{
    char *ptr;
    ptr = (char*) sf_alloc (n,sizeof(char));
    return ptr;
}

/*@out@*/ unsigned char *sf_ucharalloc (size_t n) 
{
    unsigned char *ptr;
    ptr = (unsigned char*) sf_alloc (n,sizeof(unsigned char));
    return ptr;
}

/*@out@*/ int *sf_intalloc (size_t n) 
{
    int *ptr;
    ptr = (int*) sf_alloc (n,sizeof(int));
    return ptr;
}

/*@out@*/ float *sf_floatalloc (size_t n) 
{
    float *ptr;
    ptr = (float*) sf_alloc (n,sizeof(float));
    return ptr;
}

/*@out@*/ float complex *sf_complexalloc (size_t n) 
{
    float complex *ptr;
    ptr = (float complex*) sf_alloc (n,sizeof(float complex));
    return ptr;
}

/*@out@*/ float complex **sf_complexalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    float complex **ptr;
    
    ptr = (float complex**) sf_alloc (n2,sizeof(float complex*));
    ptr[0] = sf_complexalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*@out@*/ bool *sf_boolalloc (size_t n) 
{
    bool *ptr;
    ptr = (bool*) sf_alloc (n,sizeof(bool));
    return ptr;
}

/*@out@*/ bool **sf_boolalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    bool **ptr;
    
    ptr = (bool**) sf_alloc (n2,sizeof(bool*));
    ptr[0] = sf_boolalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*@out@*/ float **sf_floatalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    float **ptr;
    
    ptr = (float**) sf_alloc (n2,sizeof(float*));
    ptr[0] = sf_floatalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*@out@*/ float ***sf_floatalloc3 (size_t n1, size_t n2, size_t n3) 
{
    size_t i3;
    float ***ptr;
    
    ptr = (float***) sf_alloc (n3,sizeof(float**));
    ptr[0] = sf_floatalloc2 (n1,n2*n3);
    for (i3=1; i3 < n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}

/*@out@*/ int **sf_intalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    int **ptr;
    
    ptr = (int**) sf_alloc (n2,sizeof(int*));
    ptr[0] = sf_intalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*@out@*/ char **sf_charalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    char **ptr;
    
    ptr = (char**) sf_alloc (n2,sizeof(char*));
    ptr[0] = sf_charalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*@out@*/ unsigned char **sf_ucharalloc2 (size_t n1, size_t n2) 
{
    size_t i2;
    unsigned char **ptr;
    
    ptr = (unsigned char**) sf_alloc (n2,sizeof(unsigned char*));
    ptr[0] = sf_ucharalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

