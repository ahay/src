#ifndef __SEAM_TRACEIO
#define __SEAM_TRACEIO

/* uncomment to force use of file manager */
#define IWAVE_USE_FMGR

#include "utils.h"
#include "usempi.h"
#include "cubic.h"
/* added 05.03.10 */
#ifdef IWAVE_USE_FMGR
#include "iwave_fopen.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "su.h"
#include "header.h"
#include "segy.h"
#ifdef __cplusplus
}
#endif

/* max number of traces permitted in record */
#define MAX_TRACES 500000  
/* max number of records permitted in data */
#define MAX_RECS 1000
/* tolerance for shot position (fraction of space step) */
#define SRC_TOL 0.001
/* end-of-file flag - should be disinct from other error returns */
#define W_EOF -10

void mygethdval(segy * tr, const char * ch, Value * val);
void myputhdval(segy * tr, const char * ch, Value * val);

typedef struct s_offsegy {
  off_t m;
  segy tr;
} offsegy;

/** \section tracegeom 

Trace Geometry. Consists of struct defining sampling etc. and
storing essential trace header and i/o status information, together
with constructor, initializer, destructor, sampler, writer, and
printer functions. The struct stores information for only one shot
record at a time. Constructor connects the object to header (trace
input) and data (trace output) files, which may be the
same. Initializer populates the struct with data from the next shot
record, read from the header file; this data includes all essential
header words and may include trace samples. Sampler interacts with a
grid, extracting samples at trace locations and storing them in the
trace buffer data member of the struct. Writer dumps the current trace
buffer to the data file. These functions rely on an internal
traceserver object, which hides details of parallelism and data
distribution.
*/

typedef struct s_tracegeom {

  /** number of shot records in data */
  int nrec;
  /** index of current record (most recently accessed) */
  int irec;
  /** index of next record to be accessed */
  int xrec;
  /** index of first record (in group) */
  int first;
  /** index of last record (in group) */
  int last;
  /** number of traces per shot record */
  /*  int ntr[MAX_RECS]; */
  int * ntr;
  /** file offset of first trace in record */
  /*  off_t recoff[MAX_RECS]; */
  off_t * recoff;
  /** source coordinates - index convention is z=0, x=1, y=2 */
  /*  RPNT src[MAX_RECS]; */
  RPNT * src;
  /** number of traces in this record with receivers in grid cells */
  int ntraces;     
  /** receiver index array in grid */
  /*  IPNT ig[MAX_TRACES]; */
  IPNT * ig;
  /** receiver relative coordinates within grid cell (betw 0 and 1) */
  /*  RPNT rg[MAX_TRACES]; */
  RPNT * rg;
  /** source index array in grid */
  IPNT is;
  /** source relative x coordinate within grid cell (betw 0 and 1) */
  RPNT rs;
  /** axis order array */
  IPNT axord;
  /** number of time samples in simulation */
  int nt;
  /** number of time samples to be output (may differ from nt) */
  int ntout;
  /** time step for simulation */
  float dt;
  /** time step for output (may differ from dt) */
  float dtout;
  /** int repn of dt in microseconds */
  int dtmus;
  /** max time */
  float tmax;
  /** min time */
  float t0;
  /** min time for output (may differ from t0) */
  float t0out;
  /** tracl - trace number within line (SEGY) */
  /*  int tracl[MAX_TRACES]; */
  int * tracl;
  /** tracr - trace number within record (SEGY) */
  /*  int tracr[MAX_TRACES]; */
  int * tracr;
  /** fldr - field record numer (SEGY) */
  /*  int fldr[MAX_TRACES]; */
  int * fldr;
  /** trace number within field record (SEGY) */
  /*  int tracf[MAX_TRACES]; */
  int * tracf;
  /** scale factor for elevations (SEGY) */
  int scalel;
  /** scale factor for offsets (SEGY) */
  int scalco;
  /** buffer for trace data */
  float * buf;
  /** interpolation flag */
  int interp;
  /** grid dimension */
  int ndim;
  /** grid cell volume */
  ireal dvol;
  /** input file pointer */
  FILE * fpin;
  /** offset of trace in file */
  /*  off_t troff[MAX_TRACES]; */
  off_t * troff;
#ifdef IWAVE_USE_MPI
  /** datatype for trace augmented with file offset */
  MPI_Datatype p;
#endif
} tracegeom;

/** Trace geometry constructor. Calling this function <ul>
    <li>accesses the input (header) file and output (data) file. These
    files may be the same. If the output file has different length
    than the input file (in particular, if it did not exist on call)
    then this function copies the input file to the output file, so
    that traces can be written in the positions dictated by the input
    file structure.</li> <li>determines the number of distinct shot
    records, the number of traces in each, and the offsets of the
    first traces in these records, and returns these to the calling
    unit.</li> <li>uses for these services the initializer of an
    internal traceserver object, which hides low-level details of data
    distribution, interaction with MPI, and i/o.</li> </ul>

    @param[out] tg - (tracegeom *) trace geometry object to be constructed
    @param[in] fin - (char *) name of header (input) file.
    @param[in] dt  - (float) time step for internal time loop, used to reserve buffer space
    @param[in] tol - (float) tolerance for distinguishing source positions with l-infty norm; should be computed in calling unit as fraction of maximum space step of grid to be sampled, possibly using cpp macro SRC_TOL defined in the header file.
    @param[in] stream - (FILE *) output stream for diagnostic information.

    @return 0 on successful completion, else error code as in base/include/utils.h. 
*/

int construct_tracegeom(tracegeom * tg,
			const char * fin,
			float dt,
			float tol,
			FILE * stream);

/** Initialization for trace geometry object. Takes spatial geometry
    and start and end time from input file. Time sample rate specified
    externally (arg dt), determines buffer allocation. Time unit is ms. It
    is ASSUMED that any common spatial data (eg. scalel, scalco) is
    uniform, and may be read from any trace.

    The data to be sampled repeatedly to form traces are assumed to
    functions on a uniform rectangular grid. The construction of
    samples differentiates between a LOCAL grid, meaning the grid on
    which the data to be sampled are being updated by the process, and
    a GLOBAL grid, of which the LOCAL grid is a subgrid. Index tuples
    are defined by the GLOBAL grid. The sampling function (below) is
    supplied with the index tuple of the LOCAL grid origin, which
    enables conversion from global indices to local offset.

    This function does several tasks, using the tracegeom struct to 
    record its results:

    - timestepping setup: records time steps, determines number of 
      time steps necessary;

    - sampling setup: determines source and receiver positions
      relative to a GLOBAL grid, and sets up sampling coefficients
      ("rectangular barycentric coordinates") for traces which lie in a
      LOCAL subgrid. This task requires two sets of origin coordinates, one 
      for the GLOBAL grid, one for the LOCAL subgrid. Also allocates buffer 
      into which samples may be recorded;
    
    - trace output setup: records headers necessary to write
      output to traces, and computes output time sampling info 
      either from the input trace header or from input arguments;

    - data input: optionally initializes sampling buffer by reading 
      data from file.

    The routine may be called several times to initialize any number
    of time series on the same internal time grid. This internal
    (stepping) time grid is fixed by the requirement tha (1) the step
    be given by the input argument dt, and (2) t=0 is a sample point.
    
    arguments:
    
    @param[out] tg (\ref tracegeom *) trace geometry object to be initialized
    @param[in] og (\ref RPNT)
       coordinates of GLOBAL grid origin (grid point with all indices = 0) 
         relative to physical
         coordinate system. The physical coordinate system refers to the
	 physical domain, and is the coordinate system with respect to
	 which the source and receiver coordinates are defined. These
	 numbers (lengths, not indices!) give the coordinates of the global
	 computational grid origin in the physical system. 
    @param[in] n (\ref IPNT) number of samples on each axis for the LOCAL subgrid in which
      traces are to be sampled.
    @param[in] o (\ref RPNT) coordinates of LOCAL subgrid origin.
    @param[in] d (\ref RPNT) spatial steps (common to GLOBAL grid and LOCAL subgrid)
    @param[in] dt (float)
      time step in ms. input argument, NOT read from file - 
        computed by the simulator, actual time step between time
        levels of simulation. Recorded in tracegeom object as tg->dt.
    @param[in] ndim (int) dimension of GLOBAL and LOCAL grids.
    @param[in] usernt,usert0 (int, float)
      optional specification of number output samples and time of first output
      sample. 
      If set (i.e. usernt > 0), then
      <ol>
      <li> output sample rate (tg->dtout) set to simulator-supplied rate (dt);</li>
      <li> number of output samples (tg->ntout) set to usernt;</li>
      <li> number of simulation time steps (tg->nt) set to usernt;</li>
      <li> time of first output sample (tg->t0out) is set to usert0, possibly
      adjusted so that t=0 is a sample;</li>
      <li> time of first simulation sample (tg->t0) also set to usert0.</li>
      </ol>
      [this option is mostly useful in convergence tests and the like].
      If NOT set (i.e. usernt <= 0) then 
      <ol>
      <li> output sample rate (tg->dtout) read from file, and used together
      with cubic interpolation (see \ref cubic) to write traces;</li>
      <li> number of output samples (tg->ntout) read from file (usernt ignored);</li>
      <li> time of first output sample (tg->t0out) read from file (usert0 ignored);</li>
      <li> number of simulation time steps (tg->nt) computed from dt, tg->dtout, tg->t0out, and tg->ntout;</li>
      <li> time of first simulation sample (tg->t0) computed from tg->t0out and dt so that 
          simulation grid includes t=0 as a sample point.</li>
</ol>
      Thus in the second case (usernt NOT set) the simulator output is interpolated
      onto the time sampling of the file used to initialized the tracegeom.
    @param[in] initbuf (int)
      flag for data initiation. If set, data traces are read
        and transferred to the internal simulation (tg->dt,tg->nt,tg->t0) grid by cubic spline interpolation (initbuf > 0) or adjoint cubic
	interpolation (initbut < 0) (see \ref cubic). To avoid temp storage, this option implemented by
	another pass through file.

	ADDED 22.10.08: AXIS ORDERING

	@param[in] axord (int) The axord argument codes axis order as follows:
	<ol>
	<li>for 1D, who cares - z index is zero in any case, so check that axord[0]=0;</li>
	<li>for 2D, axord[0]=index of z, axord[1]=index of x - check that values are legit;</li>
	<li>for 3D, axord[0]=index of z, axord[1]=index of x, axord[2]=index of y - check that values are legit.
</ol>

@param[in] order (int) interpolation order for non-grid source and receiver positions - current legal values are 0 and 1
      
@param[in] stream (FILE *) verbose output stream

@return 0 on successful completion, else error code as in base/include/utils.h. 
 */
// added irec - external control of record 11.01.14
int init_tracegeom(tracegeom * tg, 
		   int irec,
		   RPNT og,
		   IPNT n, RPNT d, RPNT o, 
		   IPNT axord,
		   int order,
		   int ndim,
		   int initbuf,
		   FILE * stream);

/** destructor: frees all dynamically allocated memory, closes
    files. Should be called as part of "default constructor"*/
void destroy_tracegeom(tracegeom * tg);

/** default initialization (as distinct from construction)- sets
    tracelength to zero, all other internals to zero or other sensible
    default values. Leaves files open, does not reset whole-data-set
    dependent quantities. Should be called to initialized structure
    for new record. */
void setnull_tracegeom(tracegeom * tg);

/** diagnostic output to stream (fp). */
void fprint_tracegeom(tracegeom const * tg, FILE * fp);

/** diagnostic output to stdout. */
void print_tracegeom(tracegeom const * tg);

/** Sampling function. To be called at every time step. Index arrays
    defined relative to global grid. Requires both local (to each 
    process) strides and origins of allocated grid, to compute offset
    into array to be sampled, and strides and origins of computational
    grid (array of points to be updated by each processor. Only those
    summands in interpolation formulas computed corresponding to these
    computational grid points - since computational grids do not overlap,
    each interpolation term is computed by at most one processor. If no
    processor computes a term, then point contributing belongs to no 
    computational grid, i.e. lies on Dirichlet for the corresponding field.
    TODO: this idea certainly works for the pressure, but needs to be 
    checked for sampling velocity in the near-boundary layer.

    Version 1.5 addendum: the final argument is a scale factor,
    applied to every sample. For load (adjoint to save), an additional
    factor of cell volume ratio is applied - that is, time step /
    volume of spatial cell.

    @param[out] tg (tracegeom) trace geometry - output for save, input for
    load. defines sampling of grid, samples saved to buffer tg.buf in save
    mode, loaded from tg.buf in load mode.
    @param[in] order (int) order of interpolation rule.- either
    trilinear (order=1) or nearest neighbor (order=0 are supported;
    @param[in] load (int) if set, transfer data from traces to grid by adj 
    interp; else, transfer data from grid to traces by interp
    @param[in] it (int) time step index for sampling - function is
    no-op if it < 0 or it > tg.nt-1;
    @param[in] allocstrides (\ref IPNT) number of samples on each
    local ALLOCATED grid axis;
    @param[in] allocorigs (\ref IPNT) index tuple of local ALLOCATED grid
    origin relative to global grid;
    @param[in] strides (\ref IPNT) number of samples on each
    local COMPUTATIONAL grid axis (updated points);
    @param[in] origs (\ref IPNT) index tuple of local COMPUTATIONAL grid 
    origin relative to global grid;
    @param[in] field (float *) (input for save, output for load) 
    buffer containing grid data to be sampled.
    @param[in] mult (ireal) scale factor for samples - for load
    (adjoint of save), scaled by ratio of cell volumes
*/
/*
  outdated material:

    Version 1.1 addendum: this function now both samples the field to the
    output trace buffer (save mode), and adjoint-samples the trace buffer
    onto the field (load) mode. The load flag determines which option is
    taken. In load mode, an additional option is provided to scale the 
    updated field values by another multiplier field and multiplication by
    the time step. This option is switched by the multiplier field pointer; 
    a null value indicates no scaling. Global grid geometry data supplied,
    since this array may not have same size and shape as sampled field, so
    offsets must be computed independently.

    @param[in] allocstrides_mult (\ref IPNT) number of samples on each
    local ALLOCATED grid axis for multiplier array;
    @param[in] allocorigs_mult (\ref IPNT) index tuple of local ALLOCATED grid
    origin relative to global grid for multiplier array;
    @param[in] strides_mult (\ref IPNT) number of samples on each
    local COMPUTATIONAL grid axis (updated points) for multiplier array;
    @param[in] origs_mult (\ref IPNT) index tuple of local COMPUTATIONAL grid 
    origin relative to global grid for multiplier array;
    @param[in] mult (float *) buffer for multiplier field. ignored on save, 
    optionally used to scale input data on load. default value is null, which
    indicates no scaling.
*/
void sampletraces(tracegeom * tg, 
		  int order,
		  int load,
		  int it,
		  IPNT allocstrides,
		  IPNT allocorigs,
		  IPNT strides,
		  IPNT origs,
		  ireal * field,
		  /*		  IPNT allocstrides_mult,
		  //		  IPNT allocorigs_mult,
		  //		  IPNT strides_mult,
		  //		  IPNT origs_mult,
		  //		  int * mult); */
		  ireal mult);

/** Trace output. Sets up internal segy for writing, by transferring
    ntout, dtout, scalel, scalco, etc. to header fields. writes out
    traces. Uses axis steps and origin coordinates combined with index
    and relative coordinate data for source and receiver locations to
    compute trace headers. On return, fname will be SU file in working
    directory.
    
   @param[in] tg (tracegeom) header and trace data to be sampled;
   @param[in] d (\ref RPNT) axis steps in grid;
   @param[in] og (\ref RPNT) coordinate tuple for origin of GLOBAL grid;
   @param[in] stream (FILE *) verbose output stream

   @return 0 on successful write, else error code as in base/include/utils.h.
*/
int writetraces(tracegeom const * tg, 
		RPNT d,
		RPNT og,
		FILE * stream);

/* 13.03.11: compute first and last record indices for current group
*/
void calc_group(int * first, int * last, int nrec);

#endif
