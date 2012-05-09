#ifndef __IWAVE_GRIDIO
#define __IWAVE_GRIDIO

/* uncomment to force use of file manager */
#define IWAVE_USE_FMGR

#ifdef SUXDR
  #include <rpc/rpc.h>
#endif
#include "grid.h"
#include "offsets.h"
#include "parser.h"
#ifdef IWAVE_USE_FMGR
  #include "iwave_fopen.h"
#endif

/** read grid from SEP77/RSF header file 
@param[out] g (grid *) - grid to be initialized
@param[in]  fname (char *) - name of RSF header file
@param[in]  fp (FILE *) - verbose output strea
@return 0 on success, else nonzero error code (see \ref utils.h)
*/
int read_grid(grid * g, char * fname, FILE * fp);

/** initialize grid from PARARRAY object 
@param[out] g (grid *) - grid to be initialized
@param[in] par (PARARRAY) - param array from which to read grid data
@param[in] fp (FILE *) - verbose output stream
@return 0 on success, else nonzero error code (see \ref utils.h)
*/
int par_grid(grid * g, PARARRAY par, FILE * fp);

/** extend by constant along axis ax an array of dimension dim defined
    by (rags,ran), assuming that its subarray defined by
    defined by (gs,n) has been initialized to correct values. If the 
    subarray is void (i.e. if (rags,ran) does not overlap (gs,n)) or if
    subarray is the entire array defined by (rags,ran), than this function
    is a no-op.

    Written dimensionally, to avoid offset computations. All axis lengths
    assumed to be correctly representable as ints. No particular attention
    paid to efficiency - it is assumed that this routine will be account 
    for an infinitesimal part of the flops of an application.
@param[out] a (ireal *) - array to be extended
@param[in] rags (IPNT) - global start indices
@param[in] ran  (IPNT) - global axis lengths
@param[in] gs  (IPNT)  - start indices of already initialized subarray
@param[in] n (IPNT)    - axis lengths of already initialized subarray
@param[in] dim (int)   - dimension of grid
@param[in]  ax (int)   - axis along which to extend by const
@return 0 on success, else nonzero error code (see \ref utils.h)
 */
int extend_array(ireal * a, IPNT rags, IPNT ran, IPNT gs, IPNT n, int dim, int ax);

/** adjoint-extend by constant along axis ax an array of dimension dim
    defined by (rags,ran) - replaces boundary elements of array by sum
    over extended axis, elements outside array by zero. If the
    subarray is void (i.e. if (rags,ran) does not overlap (gs,n)) or
    if subarray is the entire array defined by (rags,ran), than this
    function is a no-op.

    Written dimensionally, to avoid offset computations. All axis lengths
    assumed to be correctly representable as ints. No particular attention
    paid to efficiency - it is assumed that this routine will be account 
    for an infinitesimal part of the flops of an application.
@param[out] a (ireal *) - array to be extended
@param[in] rags (IPNT) - global start indices
@param[in] ran  (IPNT) - global axis lengths
@param[in] gs  (IPNT)  - start indices of already initialized subarray
@param[in] n (IPNT)    - axis lengths of already initialized subarray
@param[in] dim (int)   - dimension of grid
@param[in]  ax (int)   - axis along which to extend by const
@return 0 on success, else nonzero error code (see \ref utils.h)
 */
int adj_extend_array(ireal * a, IPNT rags, IPNT ran, IPNT gs, IPNT n, int dim, int ax);

/** read array from SEP77/RSF file structure. 
 *
 *  Preconditions:
 *  <ol>
 *  <li>file describes SEP77/RSF data structure, including (n,d,o) grid
 *  parameters, "in" key with data file name as value, and
 *  "data_format" key with either "native_float" or "xdr_float" as
 *  parameters.</li>
 *  <li>file pointed to by "in=" contains float data, either native or
 *  xdr. No other data types are currently admitted.</li>
 * </ol>
 *  Postconditions:
 *  <ol>
 *  <li>intersection of array and file grids computed; data from
 *  intersection read from file into sub-rarray
 *  corresponding to intersection</li>
 *  </ol>
 *  @param[out]  a (ireal *)     - array to be read
 *  @param[in]     gs (IPNT)     - global indices of axis starts
 *  @param[in]    n  (IPNT)      - global axis lengths
 *  @param[in]   fname  (char *) - file from which to read data
 *  @param[in]     extend (int)  - extension flag - extend along all axes in decreasing axis order if set
 *  @param[in]   fp (FILE *)     - verbose output parameter
 *  @param[in]   panelindex (int) - panel index of extended model (always be 0 for non-extended model) 
 *  @return 0 on success, else nonzero error code (see \ref utils.h)
 */
int rsfread(ireal * a, 
	    IPNT gs, 
	    IPNT n,
	    char * fname, 
	    int extend, 
	    FILE * fp,
	    int panelindex       /* D.S. 01.01.11: extended-model related*/
	    );

/** write array to RSF file structure
 *
 *  Preconditions:
 *  <ol>
 *  <li>file describes SEP77/RSF data structure, including (n,d,o) grid
 *  parameters, "in" key with data file name as value, and
 *  "data_format" key with either "native_float" or "xdr_float" as
 *  parameters.</li>  
 *  <li>file pointed to by "in=" contains float data, either native or
 *  xdr. No other data types are currently admitted. Correct string
 *  signifying type (native_float or xdr_float) submitted as arg type.</li>
 *  <li>files must both exist; data file will be overwritten, but header
 *  file is left intact. [Major change 19.01.11 WWS]</li>
 *  </ol>
 *  Postconditions:
 *  <ol>
 *  <li>intersection of array and file grids computed; data from
 *  intersection written from sub-rarray into file sector
 *  corresponding to intersection, </li>
 *  </ol>
 *  @param[in]   a (ireal *)     - array to be written
 *  @param[in]   gs (IPNT)       - global indices of axis starts
 *  @param[in]   n (IPNT)        - global axis lengths
 *  @param[in]   fname (char *)  - file to which to write data
 *  @param[in]     extend (int)  - extension flag - adjoint-extend along all axes in increasing axis order if set
 *  @param[in]   fp (FILE *)     - verbose output parameter
 *  @param[in]   panelindex (int) - panel index of extended model (always be 0 for non-extended model) 
 *  @return 0 on success, else nonzero error code (see \ref utils.h)
 */
int rsfwrite(ireal * a, 
	     IPNT gs, 
	     IPNT n, 
	     char * fname, 
	     int extend,
	     FILE * fp,
	     int panelindex        /* D.S. 01.01.11: extended-model related*/
	     );

#endif /* __SEAM_GRIDIO */
