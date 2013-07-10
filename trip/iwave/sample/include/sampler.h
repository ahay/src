#ifndef __IWAVE_SAMPLER__
#define __IWAVE_SAMPLER__

#include "trace_term.h"

/** \section sampler Sampler

Sampler Interface. Defines lightweight wrapper around trace handler. Principal purpose: enable transfer of trace data to multiple grids representing the same field. Confusingly, the adjoint is simply transfer from one copy of the grid to the trace, since all copies of the grid are the same.

The struct records the sample index array and (optional) multiplier index array (used in run method), offset coordinate vector (used in intialization), and load flag (used in both).
*/

typedef struct {
  TRACE_TERM t;
  IPNT sindex;
  RPNT mindex;
  RPNT scoord;
  int load;
} SAMPLER;

/** constructor is mostly pass-through to \ref traceterm_construct,
    also initialises struct data.

    @param[out] t - (SAMPLER *) sampler object to be constructed.
    @param[in] par - (PARARRAY *) parameter array containing filenames etc.
    @param[in] sindex - (IPNT) up to dim indices of fields to be sampled
    @param[in] mindex - (RPNT) coefficienqt array for linear comb of fields spec'd in sindex
    @param[in] scoord - (RPNT) coordinate offset vector, applied to all s/r positions
    @param[in] load - (int) if set, run method loads data from traces into grid, else saves data from grid into traces. [mod 02/13] Can be set two ways: if > 0, loads data into internal \ref tracegeom buffer using cubic spline interpolation onto internal time grid; if < 0, uses adjoint cubic spline interpolation.
    @param[in] hdrkey - (const char *) key identifying header file name in parameter table 
    @param[in] datakey - (const char *) key identifying data file name in parameter table 
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
 */
int sampler_construct(SAMPLER * s,
		      PARARRAY * par,
		      const IPNT sindex,
		      const RPNT mindex,
		      const RPNT scoord,
		      int load,
		      const char * hdrkey,
		      const char * datakey,
		      FILE * stream);

/** post-construction initialization - mostly a pass-through to 
    traceterm_init.

    @param[out] t - (SAMPLER *) to be initialized
    @param[in] m - (IMODEL *) conveys spatial grid information
    @param[in] par - (PARARRAY *) parameter array containing filenames etc.
    @param[in] stream = (FILE *) verbose output stream


    @return 0 on successful completion, else error code
*/
int sampler_init(SAMPLER * s,
		 IMODEL * m,
		 PARARRAY * par, 
		 FILE * stream);

/** pass-through to traceterm_destroy 

    @param[out] t - (SAMPLER *) to be destroyed

    @return 0 on successful completion, else error code
*/

int sampler_destroy(SAMPLER * s);

/** in load case, sets field index to each element of index
    array in turn, checks to see that this index has not been
    visited before, and if not delegates to traceterm_run. For
    save case, straight delegation, and samples only the first 
    sindex component, ie. the array with index sindex[0]. 

    @param[out] t - (SAMPLER *) trace data to be sampled
    @param[in] m - (IMODEL *) grid data to be sampled

    @return 0 on successful completion, else error code
*/
int sampler_run(SAMPLER * s, IMODEL * m);

/** Print method - delegates to \ref traceterm_fprint

    @param[out] t - (SAMPLER *) to be printed
    @param[in] stream = (FILE *) verbose output stream

*/
void sampler_fprint(SAMPLER * s, FILE * stream);

#endif
