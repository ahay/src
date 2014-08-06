/** @file
 */
#ifndef __IWAVE_MPI_GRIDIO
#define __IWAVE_MPI_GRIDIO

#ifdef SUXDR
#include <rpc/rpc.h>
#endif
#include "grid.h"
#include "gridio.h"
#include "offsets.h"
#include "parser.h"
#include "usempi.h"
#define N_SLOWEST_AXIS_SLICE 1  /**< number of data slice read in the slowest axis each time */

#ifdef IWAVE_USE_MPI

/** read rarray from SEP77/RSF file structure by  
    reading a chunk of data, then broadcasting to every process 
 *  
 *  @param[in]   panelindex (int) - panel index of extended model (always be 0 for non-extended model) 
*/
int mpirsfread(ireal * a, const IPNT _gs, const IPNT _n, 
	       const char * fname, float scale, FILE * fp,
	       int panelindex      /* D.S. 01.01.11: extended-model related */
	       );

int mpirsfwrite(ireal * a, const IPNT _gs, const IPNT _n, 
		const char * fname, float scale, FILE * fp,
	       int panelindex      /* D.S. 01.01.11: extended-model related */
	       );
#endif /* IWAVE_USE_MPI */

#endif /* __IWAVE_MPI_GRIDIO */
