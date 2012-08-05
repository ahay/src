#ifndef __XWAVE_FD_FUN__
#define __XWAVE_FD_FUN__

#include "model.h"

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
int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model);

int fd_getndim(PARARRAY * pars, FILE * stream, IMODEL * model, int * ndim);

int fd_modeldest(IMODEL * model);

int fd_isdyn(IMODEL * model, int i);

/**
 * Computes computational domain size for dynamic arrays
 *
 * @param [in] cdims - (IPNT) cartesian grid dimension in MPI communicator.
 * @param [in] crank - (IPNT) cartesian rank of this processor in MPI communicator.
 * @param [out] dgs,dge  - (IPNT []) start and end indices for arrays.
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h
 */

/**
 * Not implemented yet
 */
/*< interface: IMODEL.minfo */
int fd_modelinfo(FILE * stream, IMODEL * model);

/**
 * Load media parameters and time step
 *
 * @param [in] stream - (FILE *) file stream for dumping run-time messages and warnings
 * @param [in, out] model - (IMODEL *) IMODEL struct pointer
 * @param [in] par - (PARARRAY *) parameter array pointer
 * @param[in]   panelindex (int) - panel index of extended model (always be 0 for non-extended model) 
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h
 */
/*< interface: IMODEL.mread */
int fd_mread    (FILE * stream, IMODEL * model, PARARRAY * par, int panelindex);
  
typedef ireal (*MED_PARA_FUN)(int n, ireal *a);

/*
 * function to read media parameters individually or in group:
 * 'medkey' contains the keywords in the parameter-table with which the 
 * model data files are associated and then used to compute the media 
 * parameters whoses indices are indicated by 'ind_med'
 * para(ind_med) = fun_array(data(medkey)) where fun_array defines the 
 * relation between the data from model files and the mat data in fd model
 * implicitly requires nfun=nmed
 * read material parameters
 */
/** Not yet implemented */
int fd_set_media(PARARRAY * par, FILE * stream, IMODEL * model, 
		 int *ind_med, int nmed, const char *medkey[], int nstr, 
		 MED_PARA_FUN *fun_array, int nfun);

/*
 * if mult turns off, perform averaging of media parameters according to 
 * the grid type 
 */
/** Not yet implemented */
int fd_avg_media(PARARRAY * par, FILE * stream, IMODEL * model, 
		 int *ind_med, int nmed, const char *medkey);
  
#endif
