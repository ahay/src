#ifndef __XWAVE_CRST_22K__
#define __XWAVE_CRST_22K__

#include "fd.h"

/** value dependency */
#define DEP_F      4
/** z-derivative dependency */
#define DEP_DFDZ   1
/** x-derivative dependency */
#define DEP_DFDX   2
/** y-derivative dependency */
#define DEP_DFDY   3
/** dependency on all derivs (Laplacian) */
#define DEP_LAPF   5

/**
 * Create stencils of a 2-2k FD scheme for a 1st order wave equations provided the 
 * dependent matrix and grid types of variables (static and dynamic) are given.
 *
 * @param [in] stream - (FILE *) file stream for dumping run-time messages and warnings
 * @param [in] ic     - (IWaveInfo const &) info struct including access to dynamic flags
 * @param [in] k      - (int) spatial order 2k
 * @param [in] ndim   - (int) grid dimension
 * @param [in] gtype  - (IPNT *) grid type array, length = m_size
 * @param [in] stendm - (int **) stencil dependency matrix, m_size x m_size
 * @param [out] sten  - (STENCIL *) STENCIL struct pointer
 */
int create_sten2_2k(FILE * stream, 
		    IWaveInfo const & ic,
		    int k, int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    int stendm[RDOM_MAX_NARR][RDOM_MAX_NARR], 
		    STENCIL * sten);

#endif
