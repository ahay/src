/* Convenience allocation programs. */
/*
Copyright (C) 2004 University of Texas at Austin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdlib.h>
/*^*/

#include "alloc.h"
#include "error.h"

#include "_bool.h"
#include "c99.h"
/*^*/

    /*@out@*/ void *sf_alloc (size_t n    /* number of elements */, 
			      size_t size /* size of one element */)
/*< output-checking allocation >*/
{
    void *ptr; 
  
    ptr = malloc (n*size);
    if (NULL == ptr)
	sf_error ("%s: cannot allocate %d bytes:", __FILE__, n*size);
    return ptr;
}

void *sf_realloc (void* ptr   /* previous data */, 
		  size_t n    /* number of elements */, 
		  size_t size /* size of one element */)
/*< output-checing reallocation >*/
{
    void *ptr2; 
  
    ptr2 = realloc (ptr,n*size);
    if (NULL == ptr2)
	sf_error ("%s: cannot reallocate %d bytes:", __FILE__, n*size);
    return ptr2;
}

/*@out@*/ char *sf_charalloc (size_t n /* number of elements */)
/*< char allocation >*/ 
{
    char *ptr;
    ptr = (char*) sf_alloc (n,sizeof(char));
    return ptr;
}

/*@out@*/ unsigned char *sf_ucharalloc (size_t n /* number of elements */)
/*< unsigned char allocation >*/ 
{
    unsigned char *ptr;
    ptr = (unsigned char*) sf_alloc (n,sizeof(unsigned char));
    return ptr;
}

/*@out@*/ int *sf_intalloc (size_t n /* number of elements */)
/*< int allocation >*/  
{
    int *ptr;
    ptr = (int*) sf_alloc (n,sizeof(int));
    return ptr;
}

/*@out@*/ float *sf_floatalloc (size_t n /* number of elements */)
/*< float allocation >*/ 
{
    float *ptr;
    ptr = (float*) sf_alloc (n,sizeof(float));
    return ptr;
}

#ifndef __cplusplus
/*^*/

/*@out@*/ float complex *sf_complexalloc (size_t n /* number of elements */) 
/*< complex allocation >*/
{
    float complex *ptr;
    ptr = (float complex*) sf_alloc (n,sizeof(float complex));
    return ptr;
}

/*@out@*/ float complex **sf_complexalloc2 (size_t n1 /* fast dimension */, 
					    size_t n2 /* slow dimension */)
/*< complex 2-D allocation, out[0] points to a contiguous array >*/ 
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

#endif
/*^*/

/*@out@*/ bool *sf_boolalloc (size_t n /* number of elements */)
/*< bool allocation >*/
{
    bool *ptr;
    ptr = (bool*) sf_alloc (n,sizeof(bool));
    return ptr;
}

/*@out@*/ bool **sf_boolalloc2 (size_t n1 /* fast dimension */, 
				size_t n2 /* slow dimension */)
/*< bool 2-D allocation, out[0] points to a contiguous array >*/
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

/*@out@*/ bool ***sf_boolalloc3 (size_t n1 /* fast dimension */, 
				  size_t n2 /* slower dimension */, 
				  size_t n3 /* slowest dimension */)
/*< bool 3-D allocation, out[0][0] points to a contiguous array >*/ 
{
    size_t i3;
    bool ***ptr;
    
    ptr = (bool***) sf_alloc (n3,sizeof(bool**));
    ptr[0] = sf_boolalloc2 (n1,n2*n3);
    for (i3=1; i3 < n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}


/*@out@*/ float **sf_floatalloc2 (size_t n1 /* fast dimension */, 
				  size_t n2 /* slow dimension */)
/*< float 2-D allocation, out[0] points to a contiguous array >*/ 
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

/*@out@*/ float ***sf_floatalloc3 (size_t n1 /* fast dimension */, 
				   size_t n2 /* slower dimension */, 
				   size_t n3 /* slowest dimension */)
/*< float 3-D allocation, out[0][0] points to a contiguous array >*/ 
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

/*@out@*/ int **sf_intalloc2 (size_t n1 /* fast dimension */, 
			      size_t n2 /* slow dimension */)
/*< float 2-D allocation, out[0] points to a contiguous array >*/  
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

/*@out@*/ char **sf_charalloc2 (size_t n1 /* fast dimension */, 
				size_t n2 /* slow dimension */) 
/*< char 2-D allocation, out[0] points to a contiguous array >*/
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

/*@out@*/ unsigned char **sf_ucharalloc2 (size_t n1 /* fast dimension */, 
					  size_t n2 /* slow dimension */)
/*< unsigned char 2-D allocation, out[0] points to a contiguous array >*/
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

/* 	$Id$	 */
