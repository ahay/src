/* 
rarray.h
Igor Terentyev.
********************************************************************************
*/
/** \file rarray.h
 * Dimension information, array type and its operations.
 * Functions:
 * - Allocation/deallocation.
 * - Creation of virtual arrays.
 * - Text and binary output of array or array slice.
 * - Element access (get/set).
  
 * Notes:
 * - Binary output depends on the definition of \ref ireal.
 * - Output to file erases preexisting file content.
 * - IMPORTANT. All array functions work for allocated array only.
 *   Declared arrays are used only to store array boundaries; allocate
 *   before performing operations on declared arrays.
 */
/*============================================================================*/

#ifndef __RARRAY_H_
#define __RARRAY_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "exchangeinfo.h"

/*----------------------------------------------------------------------------*/
/**
 * Information about one dimension of the array.
 *
 * Note: 
 * - n is non-negative and n = ge - gs + 1.
 * - n0 is non-negative and n0 = ge0 - gs0 + 1.
 */
typedef struct
{
  /** the working (computational virtual) array size */
  int n;    
  /** global start index of the working (computational virtual) array */
  int gs;   
  /** global end index of the working (computational virtual) array */
  int ge;   
  /** size of the allocated array, [used for stride] */
  int n0;   
  /** global start index of the allocated array [used for bounds checks] */
  int gs0;  
  /** global end index of the allocated array [used for bounds checks] */
  int ge0;  
} INFODIM;
/*----------------------------------------------------------------------------*/
/**
 * Array type.
 * Layout: dimensions of smaller index are faster.
 *
 */
typedef struct
{
  /** number of dimensions, e.g.\ 1,2,3 */
  int ndim; 
  /** beginning pointer of the working (computational virtual) array */
  ireal *_s;
  /** beginning pointer of the allocated array */            
  ireal *_s0;      
  /** 
   * dimension information. RARR_MAX_NDIM: 
   * maximum dimension number 
   */
  INFODIM _dims[RARR_MAX_NDIM];  
  /* added 03.11.12 WWS: multidim array access */
  /* revised 03.13 WWS: offset index ranges [s0,e0]*/
  ireal * restrict _s1;
#if RARR_MAX_NDIM > 1
  ireal ** restrict _s2;
  ireal ** _s02;
#endif
#if RARR_MAX_NDIM > 2
  ireal *** restrict _s3;
  ireal *** _s03;
#endif
#if RARR_MAX_NDIM > 3
  ireal **** restrict _s4;
  ireal **** _s04;
#endif
#if RARR_MAX_NDIM > 4
  ireal ***** restrict _s5;
  ireal ***** _s05;
#endif
} RARR;
/*----------------------------------------------------------------------------*/
/*
Functions' parameters.

RARR *arr        :  array pointer.
int ndim         :  number of dimensions.
int idim         :  dimension number.
IPNT gs, ge      :  global indices of the array start, end.
IPNT ls, le      :  local indices of the array start, end.
IPNT n           :  array size.
IPNT os, oe      :  array start, end offsets (>0 - inside, <0 - outside).
IPNT gi          :  global indices.
IPNT li          :  local indices.
ireal r           :  array value.
FILE* stream     :  output stream.
const char* path :  file name.

int return       :  error code.
*/
/*----------------------------------------------------------------------------*/
/** 
 * Set all fields in an array to zeros (NO DEALLOCATION). 
 * @param[out] arr - (RARR *)
 * @return 0
 */
int ra_setnull(RARR *arr);
/*----------------------------------------------------------------------------*/
/* 
Create array (STORAGE ALLOCATION).
Calls ra_declare and then ra_allocate.
*/
typedef int (RA_CREATE_FUN)(RARR *arr, int ndim, IPNT v1, IPNT v2);
/*----------------------------------------------------------------------------*/
/** 
 * Create array (given gs and n). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_create_s(RARR *arr, int ndim, IPNT gs, IPNT n);
/*----------------------------------------------------------------------------*/
/** 
 * Create array (given ge and n). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_create_e(RARR *arr, int ndim, IPNT ge, IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Create array (given gs and ge). 
 *
 * Declaration + allocation. Set n=n0, gs=gs0, ge=ge0 and _s=_s0.
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_create(RARR *arr, int ndim, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/*
Declare array. Works like create, but does not allocate memory.
Use ra_allocate to allocate memory for declared array.
*/
/** 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_declare_s(RARR *arr, int ndim, IPNT gs, IPNT n);
/*----------------------------------------------------------------------------*/
/** 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @param [in] n  - (IPNT) sizes of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_declare_e(RARR *arr, int ndim, IPNT ge, IPNT n);
/*----------------------------------------------------------------------------*/
/** 
 * Declare array (given gs and n). 
 * Works like create, but does not allocate memory.
 * Also set n=n0, gs=gs0, ge=ge0.
 * Use ra_allocate to allocate memory.
 *
 * @param [out] arr - (RARR *)
 * @param [in] ndim - (int) number of dimensions
 * @param [in] gs - (IPNT) global start indices of the array in all dimensions
 * @param [in] ge - (IPNT) global end indices of the array in all dimensions
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_declare(RARR *arr, int ndim, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/**
 * Allocate array. 
 * Allocate memory. Set _s=_s0.
 * 
 * @param [out] arr - (RARR *)
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_allocate(RARR *arr);
/*----------------------------------------------------------------------------*/
/** 
 * Destroy array (STORAGE DEALLOCATION). 
 * Free the allocated memory pointed by _s0
 * 
 * @param[out] arr - (RARR *)
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_destroy(RARR *arr);
/*----------------------------------------------------------------------------*/
/* 
Set array (NO STORAGE ALLOCATION). 
*/
typedef int (RA_SET_FUN)(RARR *arr, const IPNT v1, const IPNT v2);
/**
 * Reset the working (computational virtual) array (given gs and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  gs  - (IPNT) global start indices of the working (computational virtual) array in all dimensions
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_greset_s(RARR *arr, const IPNT gs, const IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Reset the working (computational virtual) array (given ge and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  ge  - (IPNT) global end indices of the working (computational virtual) array in all dimensions
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_greset_e(RARR *arr, const IPNT ge, const IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Reset the working (computational virtual) array (given gs and ge) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  gs  - (IPNT) global start indices of the working (computational virtual) array in all dimensions.
 * @param[in]  ge  - (IPNT) global end indices of the working (computational virtual) array in all dimensions.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_greset(RARR *arr, const IPNT gs, const IPNT ge);
/*----------------------------------------------------------------------------*/
/**
 * Reset the working (computational virtual) array (given os and n) (NO STORAGE ALLOCATION).
 *
 * new_gs = gs + os, new_ge = ge - oe.
 * @param[out] arr - (RARR *)
 * @param[in]  os  - (IPNT) start index offsets (forward) of the working (computational virtual) array in all dimensions.
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_offset_s(RARR *arr, const IPNT os, const IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Reset the working (computational virtual) array (given oe and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  oe  - (IPNT) end index offsets (backward) of the working (computational virtual) array in all dimensions.
 * @param[in]  n   - (IPNT) sizes of the working (computational virtual) array in all dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_offset_e(RARR *arr, const IPNT oe, const IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Reset the working (computational virtual) array (given oe and n) (NO STORAGE ALLOCATION).
 *
 * @param[out] arr - (RARR *)
 * @param[in]  os  - (IPNT) start index offsets (forward) of the working (computational virtual) array in all dimensions.
 * @param[in]  oe  - (IPNT) end index offsets (backward) of the working (computational virtual) array in all dimensions.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_offset(RARR *arr, const IPNT os, const IPNT oe);
/*----------------------------------------------------------------------------*/
/** 
 * Dump array information.
 *
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 if successful, else error code as in base/include/utils.h.
 */
int ra_dump(const RARR *arr, FILE* stream);
/*----------------------------------------------------------------------------*/
/**
 * Output the working (computational virtual) array to stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_print(RARR *arr, FILE* stream);
/*----------------------------------------------------------------------------*/
/**
 * Output the working (computational virtual) array to a file
 *
 * Format: formatted ASCII
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_fprint(RARR *arr, const char *path);
/*----------------------------------------------------------------------------*/
/**
 * Output the working (computational virtual) array to stream
 *
 * Format: binary
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_write(RARR *arr, FILE* stream);
/*----------------------------------------------------------------------------*/
/**
 * Output the working (computational virtual) array to a file
 *
 * Format: binary
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_fwrite(RARR *arr, const char *path);
/*----------------------------------------------------------------------------*/
/**
 * read the working (computational virtual) array from a binary stream
 *
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_read(RARR *arr, FILE* stream);
/*----------------------------------------------------------------------------*/
/**
 * read the working (computational virtual) array from a binary file
 *
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_fread(RARR *arr, const char *path);
/*----------------------------------------------------------------------------*/
/**
 * Output a slice of the working (computational virtual) array to a stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_printslice(RARR *arr, FILE* stream, int idim, int islice);
/*----------------------------------------------------------------------------*/
/**
 * Output a slice of the working (computational virtual) array to a stream.
 *
 * Format: formatted ASCII.
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_fprintslice(RARR *arr, const char *path, int idim, int islice);
/*----------------------------------------------------------------------------*/
/**
 * Output a slice of the working (computational virtual) array to a file.
 *
 * Format: binary.
 * @param [in] arr - (RARR *)
 * @param [in] stream - (FILE *)
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_writeslice(RARR *arr, FILE* stream, int idim, int islice);
/*----------------------------------------------------------------------------*/
/**
 * Output a slice of the working (computational virtual) array to a file.
 *
 * Format: binary.
 * @param [in] arr - (RARR *)
 * @param [in] path - (const char *) file name
 * @param [in] idim - (int) dimension number. the idim'th index is fixed
 * @param [in] islice - (int) the fixed index 
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_fwriteslice(RARR *arr, const char *path, int idim, int islice);
/*----------------------------------------------------------------------------*/
/* 
Get/set value. 
If index is out of bounds, these functions:
  1) Do not perform get/set operation.
  2) Write error info into stderr.
*/
/*----------------------------------------------------------------------------*/
/**
 * Get value at a local index relative to gs.
 *
 * [No difference with ra_gget, since gs is alway 0]
 * @param [in] arr - (RARR *)
 * @param [in] li - (IPNT) local index relative to gs
 * @return the value at the specified entry, else REAL_NAN.
 */
ireal ra_get(const RARR *arr, IPNT li);
/*----------------------------------------------------------------------------*/
/**
 * Get value at a global index.
 *
 * @param [in] arr - (RARR *)
 * @param [in] gi - (IPNT) global index
 * @return the value at the specified entry, else REAL_NAN.
 */
ireal ra_gget(const RARR *arr, IPNT gi);
/*----------------------------------------------------------------------------*/
/**
 * Set value at a local index relative to gs.
 *
 * [No difference with ra_gset, since gs is alway 0]
 * @param [out] arr - (RARR *)
 * @param [in] li - (IPNT) local index relative to gs
 * @param [in] r - (ireal) the value to be set
 */
void ra_set(RARR *arr, IPNT li, ireal r);
/*----------------------------------------------------------------------------*/
/**
 * Set value at a global index.
 *
 * @param [out] arr - (RARR *)
 * @param [in] gi - (IPNT) global index
 * @param [in] r - (ireal) the value to be set
 */
void ra_gset(RARR *arr, IPNT gi, ireal r);
/*----------------------------------------------------------------------------*/
/**
 * Get axis lengths of the working (computational virtual) array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] n - (IPNT)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_size(const RARR *arr, IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Get size of the working (computational virtual) array.
 *
 * @param [in] arr - arg
 * @param [out] n -  size of array = prod of axis lengths
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_datasize(const RARR *arr, size_t * n);
/*----------------------------------------------------------------------------*/
/**
 * Get axis lengths of the allocated array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] n - (IPNT)
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_size(const RARR *arr, IPNT n);
/*----------------------------------------------------------------------------*/
/**
 * Get size of the allocated array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] n -  size of array = prod of axis lengths
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_datasize(const RARR *arr, size_t * n);
/*----------------------------------------------------------------------------*/
/**
 * Get the global start and end indices of the working (computational virtual) array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] gs - (IPNT) start indices
 * @param [out] ge - (IPNT) end indices
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_gse(const RARR *arr, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/**
 * Get the start and end indices of the working (computational virtual) array 
 * RELATIVE to the allocated array. That is, return s = gs-gs0, e = ge-gs0.
 * These are the correct index limits to loop over the computational array by
 * indexing into the multidimensinal array representation of the allocated array. 
 * For example, if a is 2D, then you can assign to its computational part by the 
 * nested loop
 * <p>
 *   ra_gse(&a,s,e);
 *   for (i[1]=s[1];i[1]<=e[1];i[1]++) 
 *     for (i[0]=s[0];i[0]<=e[0];i[0]++) 
 *       a._s02[i[1]][i[0]] = ...
 * <p>
 * @param [in] arr - (RARR *)
 * @param [out] s - (IPNT) start indices
 * @param [out] e - (IPNT) end indices
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_se(const RARR *arr, IPNT s, IPNT e);
/*----------------------------------------------------------------------------*/
/**
 * Get the start and end indices of the allocated array.
 *
 * @param [in] arr - (RARR *)
 * @param [out] gs - (IPNT) start indices
 * @param [out] ge - (IPNT) end indices
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_gse(const RARR *arr, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/** 
 * Relative index shift from source rarray to target array. Used to align
 * multidimensional index references: for example, if ds is retrieved by 
 * <p>
 * ra_ds(&tgt,&src,ds)
 * <p>
 * and tgt, src are 3D rarrays, then   
 * <p>
 * tgt._s3[i[2]+ds[2]][i[1]+ds[1]][i[0]+ds[0]]
 * <p>
 * and
 * <p>
 * src._s3[i[2]][i[1]][i[0]]
 * <p>
 * refer to values at the same global grid location
 *
 * @param [in] tgt - (RARR *) target rarray (shifted indices)
 * @param [in] src - (RARR *) source rarray (unshifted indices)
 * @param [out] ds - (IPNT) index shift from source to target
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_ds(const RARR *tgt, const RARR * src, IPNT ds);
/*----------------------------------------------------------------------------*/
/**
 * Check if a local index relative to gs (and idim) is within bounds of the allocated array.
 * 
 * @param [in] arr - (RARR *)
 * @param [in] idim - (int) dimension number.
 * @param [in] li - local index relative to gs
 * @return 0 if within bounds, else error code as in base/include/utils.h.
 */
int ra_checkbound(const RARR *arr, int idim, int li);
/*----------------------------------------------------------------------------*/
/**
 * Check if a gloabal index (and idim) is within bounds of the allocated array.
 *
 * @param [in] arr - (RARR *)
 * @param [in] idim - (int) dimension number. 
 * @param [in] gi - global index
 * @return 0 if within bounds, else error code as in base/include/utils.h.
 */
int ra_gcheckbound(const RARR *arr, int idim, int gi);
/*----------------------------------------------------------------------------*/
/**
 * Get number of dimensions.
 *
 * @param [in] arr - (RARR *)
 * @param [out] ndim - (int *) number of dimensions
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_ndim(const RARR *arr, int *ndim);
/*----------------------------------------------------------------------------*/
/**
 * Array empty query.
 *
 * @param [in] arr - (RARR *)
 * @param [out] empty - (int *) 0: nonempty, 1: empty
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_empty(RARR *arr, int *empty);
/*----------------------------------------------------------------------------*/
/**
 * Set the working (computational virtual) array empty
 * @param [in] arr - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_setempty(RARR *arr);
/*----------------------------------------------------------------------------*/
/**
 * Populates exchange info. 
 *
 * Creates MPI_Datatype inside - do not forget to destoy.
 * @param[in] arr - (RARR *)
 * @param[out] einfo - (EXCHANGEINFO *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_setexchangeinfo(RARR *arr, EXCHANGEINFO *einfo);
/*----------------------------------------------------------------------------*/
/**
 *  Checks if the two working (computational virtual) arrays of arr1
 *  and arr2 overlap.
 * 
 * @param [in] arr1, arr2 - (RARR *)
 * @param [out] overlap - (int *) 0: not overlap, 1: overlap
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_overlap(RARR *arr1, RARR *arr2, int *overlap);
/*----------------------------------------------------------------------------*/
/**
 * Set the working (computational virtual) array's dimension info of
 * arr1 to be that of the overlap part of arr1 and arr2.
 *
 * @param [in,out] arr1 - (RARR *)
 * @param [in] arr2 - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_setoverlap(RARR *arr1, RARR *arr2);
/*----------------------------------------------------------------------------*/
/**
 * Set the entries of the computational array all zero
 * @param[in] arr - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_zero(RARR *arr);
/*----------------------------------------------------------------------------*/
/**
 * Set the entries of the allocated array all zero
 * @param[in] arr - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_zero(RARR *arr);
/*----------------------------------------------------------------------------*/
/**
 * Copy all of the members of the source struct to the target struct 
 *
 * No allocation occurs. tgt works as a reference of src.
 * 
 * @param [in] src - (RARR const *) src cannot be pointed to other array, but one 
 *                    can change the array pointed by it
 * @param [out] tgt - (RARR *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_deepcopy(RARR const * src, RARR *tgt);
/*----------------------------------------------------------------------------*/
/**
 * Copy the contents of the source VIRTUAL array to the destination
 * VIRTUAL array. Both arrays must have the same structure and be
 * allocated.
 *
 * @param [in] src - (RARR const *) data copied from 
 * @param [out] tgt - (RARR *) data copied to
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_copy(RARR *arr_des, RARR *arr_src);
/*----------------------------------------------------------------------------*/
/**
 * Copy the contents of the source ALLOCATED arry to the destination
 * ALLOCATED array.  Both arrays must have the same structure and be
 * allocated.
 *
 * @param [in] src - (RARR const *) data copied from 
 * @param [out] tgt - (RARR *) data copied to
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_copy(RARR *arr_des, RARR *arr_src);
/*----------------------------------------------------------------------------*/
/**
 * scale ALLOCATED array elements
 * @param [in] tgt - (RARR *) target data
 * @param [in] fac - (ireal) scale factor
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_scale(RARR * tgt, ireal fac);
/*----------------------------------------------------------------------------*/
/**
 * part of infrastructure to give RDOM minimal DataContainer characteristics
 *
 * note that inner product is over entire ALLOCATED array, not just 
 * virtual array
 *
 * @param [in] arr1 - (RARR const *) input array 1
 * @param [in] arr2 - (RARR const *) input array 2
 * @param [out] ip  - (ireal) inner product - unscaled l2
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int ra_a_inner(RARR const * arr1, RARR const * arr2, ireal * ip);
/*----------------------------------------------------------------------------*/
/**
 * part of infrastructure to give RDOM minimal DataContainer characteristics
 *
 * saxpy for virtual arrays
 *
 * @param [in] arrx - (RARR const *) input array
 * @param [in] arry - (RARR *) input/output array = y <- ax+y
 * @param [in] a - (ireal) multiplier
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */  
int ra_axpy(RARR * arry, RARR const * arrx, ireal a);
/*-----------------------------------------------------------------------------*/
/** 
 * comparison function - required to make RARR behave like a proper metadata 
 * object
 * @return 0 if arguments have same ndim, gs and ge, 1 otherwise
 */
int ra_compare_meta(const RARR * a, const RARR * b); 
/*-----------------------------------------------------------------------------*/

#endif /*__RARRAY_H_*/
