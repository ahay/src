#ifndef __XWAVE_CRST_22K__
#define __XWAVE_CRST_22K__

#include "fd.h"

/**
 * Create stencils of a 2-2k FD scheme for a 1st order wave equations provided the 
 * dependent matrix and grid types of variables (static and dynamic) are given.
 *
 * @param [in] stream - (FILE *) file stream for dumping run-time messages and warnings
 * @param [in] k      - (int) spatial order 2k
 * @param [in] ndim   - (int) grid dimension
 * @param [in] m_size - (int) number of arrays in rdomain
 * @param [in] gtype  - (IPNT *) grid type array, length = m_size
 * @param [in] stendm - (int **) stencil dependency matrix, m_size x m_size
 * @param [in] isdyn  - (int (*)(int)) - fcn ptr: returns true if index pts to dyn arr
 * @param [out] sten  - (STENCIL *) STENCIL struct pointer
 */
int create_sten2_2k(FILE * stream, 
		    //		    int k, int ndim, int m_size, 
		    int k, int ndim,
		    IPNT gtype[RDOM_MAX_NARR], int stendm[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    int (*isdyn)(int), 
		    //		    int (*getindices)(int),
		    STENCIL * sten);

#endif
