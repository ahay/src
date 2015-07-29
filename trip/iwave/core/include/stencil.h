/* 
stencil.h
Igor Terentyev.
********************************************************************************
*/
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
/*============================================================================*/

#ifndef __STENCIL_H_
#define __STENCIL_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
/*----------------------------------------------------------------------------*/
/*
Mask structure.
Indices in the index set are of ip array.

int ip  :  index of participating array.
int ir  :  index of recomputed array.
int n   :  index set size.
IPNT *s :  index set.
*/
/**
 * Mask structure.
 * Indices in the index set are of ip array.
 */
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
/*----------------------------------------------------------------------------*/
/*
Stencil structure.

int nmask           :  number of masks in the mask array.
STENCIL_MASK *masks :  mask array.
*/
/**
 * Stencil structure.
 */
typedef struct
{
  /** number of masks in the mask array */
  int nmask;
  /** mask array */
  STENCIL_MASK *masks;
} STENCIL;
/*----------------------------------------------------------------------------*/
/*
STENCIL_MASK *mask :  mask pointer.
int i              :  set element index./
IPNT ind           :  set element (copied as data, not as pointer).
*/
/**
 * Set all fields in a mask to zeros (NO DEALLOCATION). 
 */
int mask_setnull(STENCIL_MASK *mask);
/**
 * Create mask.
 *
 * @param [out] mask - (STENCIL_MASK *) mask pointer
 * @param [in]  ip - (int) index of participating array
 * @param [in]  ir - (int) index of recomputed array
 * @param [in]  n - index set size 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int mask_create(STENCIL_MASK *mask, int ip, int ir, int n); /* ALLOCATION */
/**
 * Destroy mask (STORAGE DEALLOCATION).
 * Free allocated memory pointed by s.
 */
int mask_destroy(STENCIL_MASK *mask);                       /* DEALLOCATION */
/**
 * Set the i'th entry in the index set to be ind.
 * mask->s[i] = ind
 */
int mask_set(STENCIL_MASK *mask, int i, const IPNT ind);
/**
 * Get the i'th entry in the index set and store it in ind.
 * ind = mask->s[i]
 */
int mask_get(STENCIL_MASK *mask, int i, IPNT ind);
/*----------------------------------------------------------------------------*/
/*
Do not destroy masks that were set/get in the stencil. 
Use sten_destroy which will destroy them.

STENCIL *sten      :  stencil pointer.
STENCIL_MASK *mask :  mask pointer.
int imask          :  mask index.
*/
/**
 * Set all fields in a stencil zeros
 */
int sten_setnull(STENCIL *sten);
/**
 * Create stencil (STORAGE ALLOCATION).
 * 
 * @param [out] sten - (STENCIL *) stencil pointer
 * @param [in] nmask - (int) number of masks
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int sten_create(STENCIL *sten, int nmask); /* ALLOCATION */
/**
 * Destroy stencil (STORAGE DEALLOCATION).
 */
int sten_destroy(STENCIL *sten);           /* DEALLOCATION, INCLUDING MASKS */
/**
 * Set the imask'th entry in the mask set to be (*mask).
 * sten->masks[imask] = *mask
 */
int sten_set(STENCIL *sten, int imask, STENCIL_MASK *mask);
/**
 * Get the imask'th entry in the mask set and store it in (*mask). 
 */
int sten_get(STENCIL *sten, int imask, STENCIL_MASK *mask);
/*----------------------------------------------------------------------------*/
/* 
Output stencil (for tests).
*/
/** 
 * Output stencil (for tests).
 * 
 * @param [in] sten - (STENCIL *) stencil pointer
 * @param [in] stream - (FILE *) file pointer
 * @param [in] ind2str_fun (const char *(int)) a function pointer transferring 
 *             an integer (int) to the corresponding field variable name (char *)
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 */
int sten_out(STENCIL *sten, FILE* stream, const char* ind2str_fun(int));
/*----------------------------------------------------------------------------*/

#endif /*__STENCIL_H_*/
