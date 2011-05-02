/* Implementation of aligned malloc and free.*/
/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#include "mm_malloc.h"

#include <stdlib.h>
/*^*/

#include <errno.h>
/*----------------------------------------------------------------------------*/

void* _mm_malloc(size_t size, size_t align)
/*< malloc >*/
{
    void *malloc_ptr;
    void *aligned_ptr;

    /* Error if align is not a power of two.  */
    if ( align & (align - 1) )
    {
        errno = EINVAL;
        return (void*)0;
    }

    if ( size == 0 ) return (void*)0;

    /* Assume malloc'd pointer is aligned at least to sizeof (void*).
    If necessary, add another sizeof (void*) to store the value
    returned by malloc. Effectively this enforces a minimum alignment
    of sizeof double. */     
    if ( align < 2 * sizeof(void*) ) align = 2 * sizeof(void*);

    malloc_ptr = malloc(size + align);
    if ( !malloc_ptr ) return (void*)0;

    /* Align  We have at least sizeof (void *) space below malloc'd ptr. */
    aligned_ptr = (void*)( ((size_t)malloc_ptr + align) & ~((size_t)align - 1) );

    /* Store the original pointer just before p. */    
    ((void**)aligned_ptr)[-1] = malloc_ptr;

    return aligned_ptr;
}
/*----------------------------------------------------------------------------*/

void _mm_free (void *aligned_ptr)
/*< free >*/
{
    if ( aligned_ptr ) free( ((void**)aligned_ptr)[-1] );
}
/*----------------------------------------------------------------------------*/

