/* 
rdomain.h
Igor Terentyev.
********************************************************************************
*/
/**
 * \file rdomain.h
 * Domain type and its operations.
 * 
 * Replicates array functions and has additional group operations.
 */
/*============================================================================*/

#ifndef __RDOMAIN_H_
#define __RDOMAIN_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "rarray.h"

/*----------------------------------------------------------------------------*/
/*
Domain.

int narr  :  number of arrays in the domain.
RARR _s[] :  array storage pointer.
*/
/**
 * Domain type.
 */
typedef struct
{
  /** number of arrays in the domain */
  int narr;
  /** array storage pointer. RDOM_MAX_NARR: 
   * maximum number of arrays in one domain.*/
  RARR _s[RDOM_MAX_NARR];
} RDOM;
/*----------------------------------------------------------------------------*/
/*
Functions' parameters.

RDOM *dom         :  domain pointer.
int narr          :  number of arrays.
int iarr          :  array index
int ndim          :  number of dimensions.
int idim          :  dimension number.
IPNT dgs[], dge[] :  array start/end global indices pointer.
IPNT dn[]         :  array sizes pointer.
IPNT gs, ge       :  array start/end global indices pointer. 
IPNT n            :  array sizes pointer.
IPNT os, oe       :  array start, end offsets (>0 - inside, <0 - outside).
IPNT li           :  local indices.
IPNT gi           :  global indices.
int islice        :  slice index.
ireal r            :  array value.
FILE* stream      :  output stream.
const char* path  :  file name.

int return        :  error code.
*/
/*----------------------------------------------------------------------------*/
/**
 * Set domain (all fields) to zeros. 
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @return 0
 */
int rd_a_setnull(RDOM *dom);
/**
 * Set all fields in the specified array of the domain to zeros
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in]  iarr - (int) array index
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h. One of the errors arises when iarr>=narr.
 */
int rd_setnull(RDOM *dom, int iarr);
/*----------------------------------------------------------------------------*/
/* 
Create domain arrays / next array (STORAGE ALLOCATION). 
*/
/**
 * Create all arrays in a domain (given start indices (dgs) and sizes (dn)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_create_s(RDOM *dom, int narr, IPNT dgs[], IPNT dn[]);
/**
 * Create all arrays in a domain (given end indices (dge) and sizes (dn)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_create_e(RDOM *dom, int narr, IPNT dge[], IPNT dn[]);
/**
 * Create all arrays in a domain (given start indices (dgs) and end indices (dge)).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_create(RDOM *dom, int narr, IPNT dgs[], IPNT dge[]);
/**
 * Create next array in a domain (given start indices (gs) and sizes (n) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_create_s(RDOM *dom, IPNT gs, IPNT n);
/**
 * Create next array in a domain (given end indices (ge) and sizes (n) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_create_e(RDOM *dom, IPNT ge, IPNT n);
/**
 * Create next array in a domain (given start indices (gs) and end indices (ge) of the next array).  
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_create(RDOM *dom, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/* 
Declare domain arrays / next array. 
*/
/**
 * Declare all arrays in a domain (given dgs and dn).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_declare_s(RDOM *dom, int narr, IPNT dgs[], IPNT dn[]);
/**
 * Declare all arrays in a domain (given dge and dn).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @param [in] dn - (IPNT []) an IPNT vector storing the sizes for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_declare_e(RDOM *dom, int narr, IPNT dge[], IPNT dn[]);
/**
 * Declare all arrays in a domain (given dgs and dge).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] narr - (int) number of arrays
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_declare(RDOM *dom, int narr, IPNT dgs[], IPNT dge[]);
/**
 * Declare the next array in a domain (given gs and n).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_declare_s(RDOM *dom, IPNT gs, IPNT n);
/**
 * Declare the next array in a domain (given ge and n).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_declare_e(RDOM *dom, IPNT ge, IPNT n);
/**
 * Declare the next array in a domain (given gs and ge).
 *
 * Works like create, but does not allocate memory.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_declare(RDOM *dom, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/*
Allocate all arrays/array.
NOTE: rd_a_allocate supresses E_ALREADYALLOC errors!
*/
/**
 * Allocate memory for all arrays in a domain.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h. NOTE: rd_a_allocate supresses E_ALREADYALLOC errors!
 */
int rd_a_allocate(RDOM *dom);
/** 
 * Allocate memory for a specified array in a domain.
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h. One of the errors arises when iarr>=narr.
 */
int rd_allocate(RDOM *dom, int iarr);
/*----------------------------------------------------------------------------*/
/** 
 * Destroy domain (STORAGE DEALLOCATION). 
 */
int rd_a_destroy(RDOM *dom);
/*----------------------------------------------------------------------------*/
/* 
Set array (NO STORAGE ALLOCATION). 
*/
/**
 * Reset all the working (computational virtual) arrays in a domain (given dgs and dge) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] dgs - (IPNT []) an IPNT vector storing the global start indices for every working (computational virtual) array
 * @param [in] dge - (IPNT []) an IPNT vector storing the global end indices for every working (computational virtual) array
 * * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_a_greset(RDOM *dom, IPNT dgs[], IPNT dge[]);
/**
 * Reset a specified working (computational virtual) arrays in a domain (given gs and n) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_greset_s(RDOM *dom, int iarr, IPNT gs, IPNT n);
/**
 * Reset a specified working (computational virtual) arrays in a domain (given ge and n) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @param [in] (n - (IPNT) sizes of the next array in this domain
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_greset_e(RDOM *dom, int iarr, IPNT ge, IPNT n);
/**
 * Reset a specified working (computational virtual) arrays in a domain (given gs and ge) 
 * (NO STORAGE ALLOCATION).
 *
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gs - (IPNT) the global start indices for the next array
 * @param [in] ge - (IPNT) the global end indices for the next array
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int rd_greset(RDOM *dom, int iarr, IPNT gs, IPNT ge);
/**
 * Reset a specified working (computational virtual) array in a domain (given os and n) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset_s.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] os  - (IPNT) start index offsets (forward) of the working (computational virtual) array
 * @param [in] 	n  - (IPNT) sizes of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 */
int rd_offset_s(RDOM *dom, int iarr, IPNT os, IPNT n);
/**
 * Reset a specified working (computational virtual) array in a domain (given oe and n) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset_e.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] oe - (IPNT) end index offsets (backward) of the working (computational virtual) array
 * @param [in] 	n  - (IPNT) sizes of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 */
int rd_offset_e(RDOM *dom, int iarr, IPNT oe, IPNT n);
/**
 * Reset a specified working (computational virtual) array in a domain (given os and oe) 
 * (NO STORAGE ALLOCATION). Refer to \ref ra_offset.
 * 
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] os  - (IPNT) start index offsets (forward) of the working (computational virtual) array
 * @param [in] oe - (IPNT) end index offsets (backward) of the working (computational virtual) array
 * @return 0 on successful completion, else error code as in base/include/utils.h. 
 */
int rd_offset(RDOM *dom, int iarr, IPNT os, IPNT oe);
/*----------------------------------------------------------------------------*/
/* Resize all arrays (NO STORAGE ALLOCATION). */
/*----------------------------------------------------------------------------*/
/* 
Dump domain information. 
*/
/**
 * Dump information of all arrays in a domain  
 */
int rd_a_dump (const RDOM *dom, FILE* stream);
/**
 * Dump information of a specified arrays in a domain
 *
 * @param [in] dom - (const RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0
 */
int rd_dump(const RDOM *dom, int iarr, FILE* stream);
/*----------------------------------------------------------------------------*/
/* 
Formatted and binary output of array / all arrays in domain.
*/
/**
 * Output all the working (computational virtual) arrays in a domain to a stream
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_a_print(RDOM *dom, FILE* stream);
/**
 * Output all the working (computational virtual) arrays in a domain to a file
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_a_fprint(RDOM *dom, const char *path);
/**
 * Output each the working (computational virtual) array in a domain to a corresponding file
 *
 * the i'th array is stored in the file named 'str(path)+str(i)'
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] path - (const char *) file name (general)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_a_fsprint(RDOM *dom, const char *path);
/**
 * Output a specified working (computational virtual) array in a domain to a stream
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_print(RDOM *dom, int iarr, FILE* stream);
/**
 * Output a specified working (computational virtual) arrays in a domain to a file
 *
 * Format: formatted ASCII
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_fprint(RDOM *dom, int iarr, const char *path);
/**
 * Output a specified working (computational virtual) array in a domain to a stream
 *
 * Format: binary
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] stream - (FILE *) file pointer
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_write(RDOM *dom, int iarr, FILE* stream);
/**
 * Output a specified working (computational virtual) array in a domain to a file
 *
 * Format: binary
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] path - (const char *) file name
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_fwrite(RDOM *dom, int iarr, const char *path);
/*----------------------------------------------------------------------------*/
/* 
Formatted and binary output of array slice to stream/file. 
*/
/**
 * Output a slice of a specified working (computational virtual) array in a domain to a stream.
 *
 * check if iarr < narr. Then call \ref ra_printslice.
 */
int rd_printslice(RDOM *dom, int iarr, FILE* stream, int idim, int islice);
/**
 * Output a slice of a specified working (computational virtual) array in a domain to a stream.
 *
 * check if iarr < narr. Then call \ref ra_fprintslice.
 */
int rd_fprintslice(RDOM *dom, int iarr, const char *path, int idim, int islice);
/**
 * Output a slice of a specified working (computational virtual) array in a domain to a binary stream.
 *
 * check if iarr < narr. Then call \ref ra_writeslice.
 */
int rd_writeslice(RDOM *dom, int iarr, FILE* stream, int idim, int islice);
/**
 * Output a slice of a specified working (computational virtual) array in a domain to a binary file.
 *
 * check if iarr < narr. Then call \ref ra_fwriteslice.
 */
int rd_fwriteslice(RDOM *dom, int iarr, const char *path, int idim, int islice);
/*----------------------------------------------------------------------------*/
/* 
Get/set value. 
If index or array number is out of bounds, these functions:
  1) Do not perform get/set operation.
  2) Write error info into stderr.
*/
/**
 * Get value at a local index relative to gs in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_get.
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] li - (IPNT) local index relative to gs
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 */

ireal rd_get(const RDOM *dom, int iarr, IPNT li);
/**
 * Get value at a global index in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_gget.
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gi - (IPNT) global index
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 */
ireal rd_gget(const RDOM *dom, int iarr, IPNT gi);
/**
 * Set value at a local index relative to gs in a specified working (computational virtual) array.
 *
 * Refer to \ref ra_set.
 * [No difference with ra_gset, since gs is alway 0]
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] li - (IPNT) local index relative to gs
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 */
void rd_set(RDOM *dom, int iarr, IPNT li, ireal r);
/**
 * Set value at a global index in a specified working (computational virtual) array
 *
 * Refer to \ref ra_gset.
 * @param [out] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [in] gi - (IPNT) global index
 * @param [in] r - (ireal) the value to be set
 * @return the value at the specified entry, else error code as in base/include/utils.h.
 */
void rd_gset(RDOM *dom, int iarr, IPNT gi, ireal r);
/*----------------------------------------------------------------------------*/
/* 
Get size, gloabal start/end indices. 
gs, ge can be NULL.
*/
/**
 * Get size of a specified working (computational virtual) array in a domain.
 *
 * Refer to \ref ra_size. 
 */
int rd_size(RDOM *dom, int iarr, IPNT n);
/**
 * Get size of a specified allocated array in a domain.
 *
 * Refer to \ref ra_a_size. 
 */
int rd_a_size(RDOM *dom, int iarr, IPNT n);
/**
 * Get the start and end indices of a specified working (computational virtual) array in a domain
 *
 * Refer to \ref ra_gse.
 */
int rd_gse(const RDOM *dom, int iarr, IPNT gs, IPNT ge);
/**
 * Get the start and end indices of a specified allocated array in a domain
 *
 * Refer to \ref ra_a_gse.
 */
int rd_a_gse(const RDOM *dom, int iarr, IPNT gs, IPNT ge);
/*----------------------------------------------------------------------------*/
/*
Get number of dimensions/arrays
*/
/**
 * Get number of dimensions of a specified array.
 * 
 * @param [in] dom - (const RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [out] ndim - (int *) number of dimensions
 * 
 * @return 0 on successful completion, 
 * if iarr >= narr, return \ref E_BADARRINDEX as in base/include/utils.h.
 *
 */
int rd_ndim(const RDOM *dom, int iarr, int *ndim);
/**
 * Get number of arrays in a domain.
 */
int rd_a_narr(const RDOM *dom, int *narr);
/*----------------------------------------------------------------------------*/
/*
Array empty query.
*/
/**
 * Set a specified working (computational virtual) array in a domain empty.
 *
 * Error arises if iarr >= narr.
 */
int rd_setempty(RDOM *dom, int iarr);
/**
 * empty query for a specified array in a domain.
 *
 * @param [in] dom - (RDOM *) domain pointer
 * @param [in] iarr - (int) array index
 * @param [out] empty - (int *) 0: nonempty, 1: empty
 * @return 0
 */
int rd_empty(RDOM *dom, int iarr, int *empty);
/*----------------------------------------------------------------------------*/
/*
Populates exchange info. Creates MPI_Datatype inside - do not forget to destroy.
*/
/**
 * Populates exchange info for a specified array in a domain. 
 *
 * Creates MPI_Datatype inside - do not forget to destroy.
 * Refers to \ref ra_setexchangeinfo and \ref IMODEL::ld_r and IMODEL::ld_s.
 * @param[in] dom - (RDOM *) domain pointer
 * @param[in] iarr - (int) array index
 * @param[out] einfo - (EXCHANGEINFO *)
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_setexchangeinfo(RDOM *dom, int iarr, EXCHANGEINFO *einfo);
/*----------------------------------------------------------------------------*/
/*
Checks if arrays overlap.
*/
/**
 * Checks if two specified working (computational virtual) arrays in two domains overlap.
 * 
 * Refers to \ref ra_overlap.
 * @param [in] dom1, dom2 - (RARR *) domain pointers
 * @param [in] iarr1, iarr2 - (int) array indices 
 * @param [out] overlap - (int *) 0: not overlap, 1: overlap
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 */
int rd_overlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2, int *overlap);
/*----------------------------------------------------------------------------*/
/*
Set arr1 := arr1 intesect arr2
*/
/**
 * Set the first working (computational virtual) array's dimension info in dom1 to be that of the overlap part 
 * of the two working (computational virtual) arrays in dom1 and dom2.
 *
 * Refer to \ref ra_setoverlap.
 * @param [in,out] dom1 - (RDOM *) domain pointer
 * @param [in] dom2 - (RDOM *) domain pointer
 * @param [in] iarr1, iarr2 - (int) array indices
 * @return 0 on successful completion, else error code as in base/include/utils.h.
 * One of the errors arises when iarr1 >= narr1 or iarr2 >= narr2.
 */
int rd_setoverlap(RDOM *dom1, int iarr1, RDOM *dom2, int iarr2);
/*----------------------------------------------------------------------------*/
/*
 * part of infrastructure to give RDOM minimal DataContainer characteristics
 *
 * note that inner product is over ALLOCATED arrays, not just 
 * virtual arrays
 *
 * @param [in] d1 - (RDOM const *) input array 1
 * @param [in] d2 - (RDOM const *) input array 2
 * @param [out] ip  - (ireal) inner product - unscaled l2
 * @return 0 - normal return
 */
int rd_a_inner(RDOM const * dom1, RDOM const * dom2, ireal * ip);

/*----------------------------------------------------------------------------*/
/*
 * part of infrastructure to give RDOM minimal DataContainer characteristics
 *
 * note: zeros ALL ALLOCATED arrays, not just virtual arrays
 *
 * @param [in] dom - (RDOM const *) input array 
 * @return 0 - normal return
 */
int rd_a_zero(RDOM * dom);

/*----------------------------------------------------------------------------*/
/*
 * part of infrastructure to give RDOM minimal DataContainer characteristics
 *
 * note: scales ALLOCATED array by index
 *
 * @param [in] dom - (RDOM const *) input array 
 * @return 0 - normal return
 */
int rd_a_scale(RDOM * dom, int iarr, ireal fac);

#endif /*__RDOMAIN_H_*/
