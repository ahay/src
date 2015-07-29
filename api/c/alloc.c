 /* Convenience allocation programs. */
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2007 Colorado School of Mines  

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

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include <stdlib.h>
/*^*/

#include "alloc.h"
#include "error.h"

#include "_bool.h"
#include "c99.h"
/*^*/

/*------------------------------------------------------------*/
/*@out@*/ void *sf_alloc (size_t n    /* number of elements */, 
			  size_t size /* size of one element */)
	  /*< output-checking allocation >*/
{
    void *ptr; 
    
    size *= n;
    
    if (0>=size) sf_error("%s: illegal allocation (%d bytes)",__FILE__,size);
    
    ptr = malloc (size);

    if (NULL == ptr)
	sf_error ("%s: cannot allocate %lu bytes:", __FILE__,size);

    return ptr;
}

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
/*@out@*/ char *sf_charalloc (size_t n /* number of elements */)
	  /*< char allocation >*/ 
{
    char *ptr;
    ptr = (char*) sf_alloc (n,sizeof(char));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ unsigned char *sf_ucharalloc (size_t n /* number of elements */)
	  /*< unsigned char allocation >*/ 
{
    unsigned char *ptr;
    ptr = (unsigned char*) sf_alloc (n,sizeof(unsigned char));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ short *sf_shortalloc (size_t n /* number of elements */)
	  /*< short allocation >*/  
{
    short *ptr;
    ptr = (short*) sf_alloc (n,sizeof(short));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ int *sf_intalloc (size_t n /* number of elements */)
	  /*< int allocation >*/  
{
    int *ptr;
    ptr = (int*) sf_alloc (n,sizeof(int));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ off_t *sf_largeintalloc (size_t n /* number of elements */)
	  /*< sf_largeint allocation >*/  
{
    off_t *ptr;
    ptr = (off_t*) sf_alloc (n,sizeof(off_t));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ float *sf_floatalloc (size_t n /* number of elements */)
	  /*< float allocation >*/ 
{
    float *ptr;
    ptr = (float*) sf_alloc (n,sizeof(float));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ sf_complex *sf_complexalloc (size_t n /* number of elements */) 
	  /*< complex allocation >*/
{
    sf_complex *ptr;
    ptr = (sf_complex*) sf_alloc (n,sizeof(sf_complex));
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ sf_complex **sf_complexalloc2 (size_t n1 /* fast dimension */, 
					 size_t n2 /* slow dimension */)
	  /*< complex 2-D allocation >*/ 
{
    size_t i2;
    sf_complex **ptr;
    
    ptr = (sf_complex**) sf_alloc (n2,sizeof(sf_complex*));
    ptr[0] = sf_complexalloc (n1*n2);
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ sf_complex ***sf_complexalloc3 (size_t n1 /* fast dimension */, 
					  size_t n2 /* slower dimension */, 
					  size_t n3 /* slowest dimension */)
	  /*< complex 3-D allocation >*/ 
{
    size_t i3;
    sf_complex ***ptr;
    
    ptr = (sf_complex***) sf_alloc (n3,sizeof(sf_complex**));
    ptr[0] = sf_complexalloc2 (n1,n2*n3);
    for (i3=1; i3 < n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ sf_complex ****sf_complexalloc4 (size_t n1 /* fast dimension */, 
					   size_t n2 /* slower dimension */, 
					   size_t n3 /* slower dimension */, 
					   size_t n4 /* slowest dimension */)
	  /*< complex 4-D allocation >*/ 
{
    size_t i4;
    sf_complex ****ptr;
    
    ptr = (sf_complex****) sf_alloc (n4,sizeof(sf_complex***));
    ptr[0] = sf_complexalloc3 (n1,n2,n3*n4);
    for (i4=1; i4 < n4; i4++) {
	ptr[i4] = ptr[0]+i4*n3;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ sf_complex *****sf_complexalloc5 (size_t n1 /* fast dimension */, 
					    size_t n2 /* slower dimension */, 
					    size_t n3 /* slower dimension */, 
					    size_t n4 /* slower dimension */,
					    size_t n5 /* slowest dimension */)
	  /*< complex 5-D allocation >*/ 
{
    size_t i5;
    sf_complex *****ptr;
    
    ptr = (sf_complex*****) sf_alloc (n5,sizeof(sf_complex****));
    ptr[0] = sf_complexalloc4 (n1,n2,n3,n4*n5);
    for (i5=1; i5 < n5; i5++) {
	ptr[i5] = ptr[0]+i5*n4;
    }
    return ptr;
}

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
/*@out@*/ float ****sf_floatalloc4 (size_t n1 /* fast dimension */, 
				    size_t n2 /* slower dimension */, 
				    size_t n3 /* slower dimension */, 
				    size_t n4 /* slowest dimension */)
/*< float 4-D allocation, out[0][0][0] points to a contiguous array >*/ 
{
    size_t i4;
    float ****ptr;
    
    ptr = (float****) sf_alloc (n4,sizeof(float***));
    ptr[0] = sf_floatalloc3 (n1,n2,n3*n4);
    for (i4=1; i4 < n4; i4++) {
	ptr[i4] = ptr[0]+i4*n3;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ float *****sf_floatalloc5 (size_t n1 /* fast dimension */, 
				     size_t n2 /* slower dimension */, 
				     size_t n3 /* slower dimension */, 
				     size_t n4 /* slower dimension */,
				     size_t n5 /* slowest dimension */)
/*< float 5-D allocation, out[0][0][0][0] points to a contiguous array >*/ 
{
    size_t i5;
    float *****ptr;
    
    ptr = (float*****) sf_alloc (n5,sizeof(float****));
    ptr[0] = sf_floatalloc4 (n1,n2,n3,n4*n5);
    for (i5=1; i5 < n5; i5++) {
	ptr[i5] = ptr[0]+i5*n4;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ float ******sf_floatalloc6 (size_t n1 /* fast dimension */, 
				      size_t n2 /* slower dimension */, 
				      size_t n3 /* slower dimension */, 
				      size_t n4 /* slower dimension */,
				      size_t n5 /* slower dimension */,
				      size_t n6 /* slowest dimension */)
/*< float 6-D allocation, out[0][0][0][0][0] points to a contiguous array >*/ 
{
    size_t i6;
    float ******ptr;
    
    ptr = (float******) sf_alloc (n6,sizeof(float*****));
    ptr[0] = sf_floatalloc5 (n1,n2,n3,n4,n5*n6);
    for (i6=1; i6 < n6; i6++) {
	ptr[i6] = ptr[0]+i6*n5;
    }
    return ptr;
}

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
/*@out@*/ int ***sf_intalloc3 (size_t n1 /* fast dimension */, 
			       size_t n2 /* slower dimension */, 
			       size_t n3 /* slowest dimension */)
/*< int 3-D allocation, out[0][0] points to a contiguous array >*/ 
{
    size_t i3;
    int ***ptr;
    
    ptr = (int***) sf_alloc (n3,sizeof(int**));
    ptr[0] = sf_intalloc2 (n1,n2*n3);
    for (i3=1; i3 < n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}

/*------------------------------------------------------------*/
/*@out@*/ int ****sf_intalloc4 (size_t n1 /* fast dimension */, 
			       size_t n2 /* slower dimension */, 
			       size_t n3 /* slower dimension */,
			       size_t n4 /* slowest dimension */ )
/*< int 4-D allocation, out[0][0][0] points to a contiguous array >*/ 
{
    size_t i4;
    int ****ptr;
    
    ptr = (int****) sf_alloc (n4,sizeof(int***));
    ptr[0] = sf_intalloc3 (n1,n2,n3*n4);
    for (i4=1; i4 < n4; i4++) {
	ptr[i4] = ptr[0]+i4*n3;
    }
    return ptr;
}

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
/*@out@*/ unsigned char ***sf_ucharalloc3 (size_t n1 /* fast dimension */, 
					   size_t n2 /* slower dimension */, 
					   size_t n3 /* slowest dimension */)
/*< unsigned char 3-D allocation, out[0][0] points to a contiguous array >*/ 
{
    size_t i3;
    unsigned char ***ptr;
    
    ptr = (unsigned char***) sf_alloc (n3,sizeof(unsigned char**));
    ptr[0] = sf_ucharalloc2 (n1,n2*n3);
    for (i3=1; i3 < n3; i3++) {
	ptr[i3] = ptr[0]+i3*n2;
    }
    return ptr;
}

