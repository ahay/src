/* Finite difference stencil. */

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

/* 
stencil.c
Igor Terentyev.
*/
/*============================================================================*/

#include <trip/base.h>

#include "stencil.h"


/**
 \section stencil
 Stencil: finite difference stencil description object.

 Stencil is described by masks.
 Each mask [ip,ir] describes which grid points of array #ip are used
 for computation of array #ir.
 It is not required to store empty masks.

 Requirements:
<ul>
 <li>Along each coordinate axis all arrays have same grid step. Grid step can
   be different for different axes.</li>

 <li> Arrays can be shifted relative to each other (e.g. staggered). Shifts are 
   not reflected in the stencil info and array indices are always integer.</li>
  
 <li> It is user's responsibility to make assumptions about shifts.</li>
</ul>

 Example: 
 Acoustic wave equations:
 \f[
   \rho(\mathbf{x})\frac{\partial\mathbf{v}}{\partial t} + \nabla p = 0,
 \f]
 \f[ 
   \frac{1}{\kappa(\mathbf{x})}\frac{\partial p}{\partial t} + 
   \nabla \cdot \mathbf{v} =  f(\mathbf{x},t;\mathbf{x}_s).
 \f]

  2-4 staggered grid, 2D.\n
 - Assumptions:
  - Pressure \f$ p \f$ with index (0,0) corresponds to ( 0 , 0 ) point.
  - Velocity \f$v_x\f$ with index (0,0) corresponds to (1/2, 0 ) point.
  - Velocity \f$v_y\f$ with index (0,0) corresponds to ( 0 ,1/2) point.
 <pre>
  E.g. (see first mask below):\n
     variable:   P   P V P   P 
        index:  -1   0 0 1   2 
      stencil:  -o---o-X-o---o-
 </pre>
  
 - Masks: \n
 <pre>
    [ip]        [ir]                                   
     \f$p\f$ mask for \f$v_x\f$: { (-1,0)  ( 0,0)  (1,0)  (2,0) }.
     \f$p\f$ mask for \f$v_y\f$: { (0,-1)  ( 0,0)  (0,1)  (0,2) }.
     \f$v_x\f$ mask for \f$p\f$: { (-2,0)  (-1,0)  (0,0)  (1,0) }.
     \f$v_y\f$ mask for \f$p\f$: { (0,-2)  (0,-1)  (0,0)  (0,1) }.
     \f$v_x\f$ mask for \f$v_y\f$: {}.
     \f$v_y\f$ mask for \f$v_x\f$: {}.
     \f$p\f$ mask for \f$p\f$: {}.
     \f$v_x\f$ mask for \f$v_x\f$: {}.
     \f$v_y\f$ mask for \f$v_y\f$: {}.
 </pre>
 */

#ifndef _sf_stencil_h

typedef struct
{
  /** index of participating array */
  int ip;
  /** index of recomputed array */
  int ir;
  /** index set size */
  int n;
  /** pointer to the index set of length n */
  IPNT *s;
} STENCIL_MASK;
/*^*/

typedef struct
{
  /** number of masks in the mask array */
  int nmask;
  /** mask array */
  STENCIL_MASK *masks;
} STENCIL;
/*^*/

#endif

/*----------------------------------------------------------------------------*/

int mask_setnull(STENCIL_MASK *mask)
/*< Set all fields in a mask to zeros (NO DEALLOCATION). >*/
{
    memset((void*)mask, 0, sizeof(STENCIL_MASK));

    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_create(STENCIL_MASK *mask, int ip, int ir, int n)
/*<*
 * Create mask.
 *
 * @param [out] mask - (STENCIL_MASK *) mask pointer
 * @param [in]  ip - (int) index of participating array
 * @param [in]  ir - (int) index of recomputed array
 * @param [in]  n - index set size 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 >*/
{
    mask_setnull(mask);                 /* empty mask */
    if ( n < 0 ) return E_BADINPUT;

    if ( n > 0 )                        /* allocate memory */
    {
        mask->s = (IPNT*)malloc(n * sizeof(IPNT));
        if ( mask->s == NULL ) return E_ALLOC;
    }
    
    mask->n = n;
    mask->ip = ip;
    mask->ir = ir;

    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_destroy(STENCIL_MASK *mask)
/*<*
 * Destroy mask (STORAGE DEALLOCATION).
 * Free allocated memory pointed by s.
 >*/
{
    free(mask->s);
    mask_setnull(mask);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_set(STENCIL_MASK *mask, int i, const IPNT ind)
/*<*
 * Get the i'th entry in the index set and store it in ind.
 * ind = mask->s[i]
 >*/
{
    if ( (unsigned int)i >= (unsigned int)(mask->n) ) return E_BADINDEX;
    IASN(mask->s[i], ind);
        
    return 0;
}
/*----------------------------------------------------------------------------*/

int mask_get(STENCIL_MASK *mask, int i, IPNT ind)
/*<*
 * Get the i'th entry in the index set and store it in ind.
 * ind = mask->s[i]
 >*/
{
    if ( (unsigned int)i >= (unsigned int)(mask->n) ) return E_BADINDEX;
    IASN(ind, mask->s[i]);
        
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_setnull(STENCIL *sten)
/*< Set all fields in a stencil zeros >*/
{
    memset((void*)sten, 0, sizeof(STENCIL));
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_create(STENCIL *sten, int nmask)
/*<*
 * Create stencil (STORAGE ALLOCATION).
 * 
 * @param [out] sten - (STENCIL *) stencil pointer
 * @param [in] nmask - (int) number of masks
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 >*/
{
    int m;
    
    sten_setnull(sten);                 /* empty stencil */
    if ( nmask < 0 ) return E_BADINPUT;

    if ( nmask > 0 )                    /* allocate memory */
    {
        sten->masks = (STENCIL_MASK*)malloc(nmask* sizeof(STENCIL_MASK));
        if ( sten->masks == NULL ) return E_ALLOC;
        for ( m = 0; m < nmask; ++m ) mask_setnull(sten->masks + m);
    }
    
    sten->nmask = nmask;
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_destroy(STENCIL *sten)
/*< Destroy stencil (STORAGE DEALLOCATION). >*/
{
    int m;
    
    for ( m = 0; m < sten->nmask; ++m ) mask_destroy(sten->masks + m);
    free(sten->masks);
    sten_setnull(sten);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_set(STENCIL *sten, int imask, STENCIL_MASK *mask)
/*<*
 * Set the imask'th entry in the mask set to be (*mask).
 * sten->masks[imask] = *mask
 >*/
{
    if ( (unsigned int)imask >= (unsigned int)(sten->nmask) ) return E_BADINDEX;
        
    sten->masks[imask] = *mask;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_get(STENCIL *sten, int imask, STENCIL_MASK *mask)
/*< Get the imask'th entry in the mask set and store it in (*mask). >*/
{
    if ( (unsigned int)imask >= (unsigned int)(sten->nmask) ) return E_BADINDEX;
        
    *mask = sten->masks[imask];
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sten_out(STENCIL *sten, FILE* stream, const char* ind2str_fun(int))
/*<* 
 * Output stencil (for tests).
 * 
 * @param [in] sten - (STENCIL *) stencil pointer
 * @param [in] stream - (FILE *) file pointer
 * @param [in] ind2str_fun (const char *(int)) a function pointer transferring 
 *             an integer (int) to the corresponding field variable name (char *)
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 >*/
{
    int err, m, i, d;
    STENCIL_MASK mask;
    IPNT ind;
    
    for ( m = 0; m < sten->nmask; ++m )
    {
        err = sten_get(sten, m, &mask);
        if ( err ) return err;
		if ( ind2str_fun )
            printf("mask %d:  ip = %s  ir = %s  n = %d\n", m,
			       ind2str_fun(mask.ip), ind2str_fun(mask.ir), mask.n);
		else
            printf("mask %d:  ip = %d  ir = %d  n = %d\n", m,
				   mask.ip, mask.ir, mask.n);
    
        for ( i = 0; i < mask.n; ++i )
        {
            err = mask_get(&mask, i, ind);
            if ( err ) return err;
            printf(" [");
            for ( d = 0; d < RARR_MAX_NDIM; ++d ) printf(" %d", ind[d]);
            printf(" ]");
        }
		printf("\n");
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
