 /* More Convenience allocation programs. */
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

#include <rsf.h>
#include "extraalloc.h"

/*^*/

float **sf_floatrealloc2 (float ** ptr,
			  size_t n1 /* fast dimension */, 
			  size_t n2 /* slow dimension */)
/*< reallocate float** array.  only n2, slowest dimension, can change 
    ptr[0][0] points to a contiguous array >*/ 
{
    size_t i2;
    float **ptr2;
    
    if (ptr[1] != ptr[0] + n1){
      sf_error("realloc2 cannot change n1. new n1=%d.",n1);
    }
    ptr2 = (float**) realloc (ptr,n2*sizeof(float*));
    ptr2[0] = realloc (ptr2[0],n1*n2*sizeof(float));
    for (i2=1; i2 < n2; i2++) {
	ptr2[i2] = ptr2[0]+i2*n1;
    }
    return ptr2;
}

char  **sf_alloc2(size_t n1bytes, /* sizeof a vector.  for a 2D float array this
				    will be n1*sizeof(float).  for an array of
				    seqy traces with will be 240 byte header +
				    nt*sizeof(float) */
                  size_t n2 /* slow dimension */)
/*< allocate a two d array.  to allocate a 2d array of floats:
    my2Dfloatarray=(float**) sf_alloc2(n1*sizeof(float),n2);
    to allocate an array of segy traces:
    myarraysegytraces=(segytrace*)sfalloc(sizeof(segytrace)+
                                             (n1-1)*sizeof(float),
					     n2); >*/

{
    size_t i2;
    char **ptr;
    
    ptr = (char**) sf_alloc (n2,sizeof(char*));
    ptr[0] = (char*)sf_alloc (n1bytes*n2,sizeof(char));
    for (i2=1; i2 < n2; i2++) {
	ptr[i2] = ptr[0]+i2*n1bytes;
    }
    return ptr;
}

char  **sf_realloc2(void** ptrin,
		    size_t n1bytes, /* sizeof a vector.  for a 2D float array 
				       this will be n1*sizeof(float).  for an 
				       array of segy traces with will be 240 
				       byte header + nt*sizeof(float) */
                    size_t n2 /* slow dimension */)
/*< reallocate char** array.  only n2, slowest dimension, can change 
    ptr[0][0] points to a contiguous array >*/ 
{
    size_t i2;
    char **ptr=(char**)ptrin;
    CHAR **ptr2;
    
    if (ptr[1] != ptr[0] + n1bytes){
      sf_error("realloc2 cannot change n1bytes. new n1bytes=%d.",n1bytes);
    }
    ptr2 = (char**) realloc (ptr,n2*sizeof(char*));
    ptr2[0] = (char*)realloc (ptr2[0],n1bytes*n2*sizeof(char));
    for (i2=1; i2 < n2; i2++) {
	ptr2[i2] = ptr2[0]+i2*n1bytes;
    }
    return ptr2;
}
		    
