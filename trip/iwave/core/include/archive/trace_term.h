#ifndef __SEAM_TRACE_TERMINATOR__
#define __SEAM_TRACE_TERMINATOR__

#include "traceio.h"
#include "parser.h"
#include "grid.h"
#include "model.h"

/** \section trace_term Trace Handler

    Struct stores total number of time steps,
    start and stop times for recording (index), a tracegeom struct,
    interpolation order, domain index of the field to be sampled, and
    a flag to indicate whether the action is load (read data from
    file, insert into grid by adj interp), or save (extract from file
    by interp, write to file). Associated functions initialize and
    manipulate this data.
*/
typedef struct {
  int istart;          /* start index for recording */
  int istop;           /* stop index for recording */
  tracegeom tg;        /* trace geometry */ 
  int order;           /* sampling order - default = 0*/
  int index;           /* sampling index - default = 0*/
  ireal mult;          /* sample multiplier - default = REAL_ONE */
  int load;            /* load flag: 0=save, 1=load, default=0 */
} TRACE_TERM;

/** Trace handler default constructor.  Calling this function default-
    constructs the \ref TRACE_TERM argument by delegation to the \ref
    traceterm_construct function. The load, sampling, and field index
    flags are reserved, and start and stop time indices given default
    values (=0). Most important role: implements options for separate
    header vs. file overwrite. If the hdrkey arg is NULL, or if it is
    given but the \ref PARARRAY does not contain a value for that key,
    then the trace headers are assumed to be present in the in- or
    out-put data file (pointed to by datakey). It is an error for the
    \ref PARARRAY not to define a value for the key datakey. If load is
    set then data file must contain trace data. If load is not set (save
    option, for recording traces) then data file may be empty on 
    on construction. 

    @param[out] t - (TRACE_TERM *) trace handler object to be constructed.
    @param[in] par - (PARARRAY *) parameter array containing filenames etc.
    @param[in] load - (int) if set, run method loads data from traces
    into grid, else saves data from grid into traces.
    @param[in] hdrkey - (const char *) key identifying header file
    name in parameter table
    @param[in] datakey - (const char *) key identifying data file name
    in parameter table
    @param[in] stream - (FILE *) verbose output stream

    @return 0 on successful completion, else error code
*/
/*
  removed 03.12
    @param[in] sindex - (int) domain index of the field to be sampled
    @param[in] mindex - (ireal) multiplier - applied to all samples
*/

int traceterm_construct(TRACE_TERM * t, 
			PARARRAY * par, 
			/* removed 03.12
			//			int sindex, 
			//			ireal mindex, */
			int load,
			const char * hdrkey,
			const char * datakey,
			FILE * stream);

/** Post-construction initialization. Extracts global and
    local grid params from \ref IMODEL argument, delegates to \ref
    init_tracegeom.
    
    @param[out] t - (TRACE_TERM *) to be initialized
    @param[in] m - (IMODEL *) conveys spatial grid information
    @param[in] par - (PARARRAY *) parameter array containing filenames etc.
    @param[in] stream = (FILE *) verbose output stream


    @return 0 on successful completion, else error code
*/
int traceterm_init(TRACE_TERM * t, IMODEL * m, PARARRAY * par, FILE * stream);

/** 
    Destructor - delegates to destroy_tracegeom 
    
    @param [out] t (TRACE_TERM *) to be destroyed

    @return 0 on successful completion, else error code
*/    
int traceterm_destroy(TRACE_TERM * t);

/** Extracts allocated and computational index limits from \ref IMODEL
    input for the \ref RARR (rarray, discrete field) corresponding to the
    stored index, passes this info to sampletraces. 

    @param[out] t - (TRACE_TERM *) trace data to be sampled
    @param[in] m - (IMODEL *) grid data to be sampled

    @return 0 on successful completion, else error code
*/
int traceterm_run(TRACE_TERM * t, IMODEL * m);

/** Print method - displays internal data of TRACE_TERM, then delegates 
    to fprint_tracegeom.

    @param[in] t - (TRACE_TERM *) to be printed
    @param[out] stream = (FILE *) verbose output stream
*/
void traceterm_fprint(TRACE_TERM const * t, FILE * stream);


#endif
