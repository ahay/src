#ifndef __XWAVE_FD_FUN__
#define __XWAVE_FD_FUN__

#include "model.h"
#include "mtuple.hh"

/*----------------------------------------------------------------------------*/
/**
 * General FD model creator (IMODEL struct), suitable for 1st order
 * wave equations FD modeling.  Mimics sgn_modelcrea but more
 * general, does the following operations:
 * 
 *    - read grid info
 *    - make action list
 *    - create stencil
 *    - compute size of the domain on its host processor according
 *    to the parallel info Cartesian grid dimensions (cdims) and
 *    Cartesian rank (crank)
 *    - allocated memory, NOTE: STORAGE ALLOCATION occurs only once,
 *    that is for \ref IMODEL::ld_a and set the dimension
 *    information for IMODEL::ld_r[i], IMODEL::ld_s[i], but let \n
 *    IMODEL::ld_r[i]._s[j]._s = IMODEL::ld_s[i]._s[j]._s =
 *    IMODEL::ld_a[i]._s[j]._s,\n where the first '_s' is RARR type
 *    and the second '_s' is (ireal *) type - etc
 *
 * @param [in] cdims  - (IPNT) cartesian grid dimensions in MPI communicator.
 * @param [in] crank  - (IPNT) cartesian rank of this processor in MPI communicator.
 * @param [in] par    - (PARARRAY *) parameter arrary pointer.
 * @param [in] stream - (FILE *) stream to output comments (created by driver).
 * @paoram [out] model - (IMODEL *)  IMODEL pointer.
 * 
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model, IWaveConnectorBase const & i);

int fd_getndim(PARARRAY * pars, FILE * stream, IMODEL * model, int * ndim);

int fd_modeldest(IMODEL * model);

//int fd_isdyn(IMODEL * model, int i);

/**
 * Computes computational domain size for dynamic arrays
 *
 * @param [in] cdims - (IPNT) cartesian grid dimension in MPI communicator.
 * @param [in] crank - (IPNT) cartesian rank of this processor in MPI communicator.
 * @param [out] dgs,dge  - (IPNT []) start and end indices for arrays.
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h
 */

  
#endif
