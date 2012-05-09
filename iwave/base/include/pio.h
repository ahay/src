#ifndef __SEAM_PIO__
#define __SEAM_PIO__

#define PIO_IOPROC 0

#include "utils.h"

/** page pio
    
Parallel versions of standard library i/o functions. All except
pio_fwrite have exactly the same signatures as the standard
counterparts, and have identical action and return values in serial
environment. In parallel environment (when IWAVE_USE_MPI is set) i/o
takes place on the i/o node (IOPROC) supplemented with appropriate
communication. The fwrite analog must specify the source process rank
in addition to the standard arguments.

In this version, all offset computations are reserved assumed to be
confined to IOPROC. Thus no offsets are communicated. Thus definition
of MPI_Offset is avoided, and only core MPI 1 features are used.
*/

/** Parallel fopen: opens file on i/o process, returns NULL on all
    other processes. */
FILE * pio_fopen(const char * restrict path, const char * restrict mode);

/** Parallel fclose: closes file on i/o process if one is open, no-op
    on other processes. */
int pio_fclose(FILE * pfp);

/** Tests whether file is attached to file pointer on i/o process. 
    Instead of testing for non-NULL, as one would with a FILE object 
    returned by fopen, use this test.*/ 
int pio_isopen(FILE * pfp);

/** Parallel fread: data read on i/o node and broadcast to other nodes
    in parallel case. */
int pio_fread(void * restrict ptr, 
	      size_t size, 
	      size_t nmemb, 
	      FILE * restrict pstream);

/** Parallel fwrite: data sent from process rank src to IOPROC, written 
    to file. */
int pio_fwrite(void * restrict ptr, 
	       size_t size, 
	       size_t nmemb, 
	       FILE * restrict pstream,
	       int src);

/** Parallel fseeko: seek happens on IOPROC only, so no need for MPI to 
    understand the off_t type */
int pio_fseeko(FILE *stream, off_t offset, int whence);

#endif
