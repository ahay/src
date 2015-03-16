/* 
   iwave.h
   Igor Terentyev.

   thoroughly mangled by WWS 01/09
   even further mangled by WWS & XW 11/10
   mashed by WWS 10/13-01/14
*/
/*============================================================================*/
/**
 * \file iwave.h
 * Wave simulation struct. 
 * 
 * Contains the \ref IMODEL struct and the parallel information.
 */
#ifndef __IWAVE_H_
#define __IWAVE_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "usempi.h"
#include "exchangeinfo.h"
#include "parser.h"
#include "fd.h"

/*----------------------------------------------------------------------------*/
/*
  Parallel information (MPI, MPI Cartesian grid, OpenMP).
  Cartesian grid for 1D/2D/3D.
*/
/**
 * Parallel information (MPI, MPI Cartesian grid, OpenMP).
 * 
 * Assumption: cartesian grid for 1D/2D/3D.
 */
typedef struct s_PARALLELINFO {
  int ndim;                /**< number of dimensions */
  int nnei;                /**< number of neighbor processors, it equals to IWAVE_NNEI */
  int wsize;               /**< Size of the MPI_COMM_WORLD */
  int wrank;               /**< Rank in the MPI_COMM_WORLD */            
  int threadsupp;          /**< Thread support level       */
  MPI_Comm ccomm;          /**< Cartesian communicator     */
  IPNT cdims;              /**< Cartesian grid dimensions  */
  IPNT crank;              /**< Cartesian rank             */
  int lrank;               /**< Linear rank in ccomm       */
  int sranks[IWAVE_NNEI];  /**< Ranks of send neighbors - subtracted from crank */
  int rranks[IWAVE_NNEI];  /**< Ranks of recv neighbors - added to crank        */
  EXCHANGEINFO seinfo[RDOM_MAX_NARR][IWAVE_NNEI]; /**< send buffers, corresponding to send neighbors */
  EXCHANGEINFO reinfo[RDOM_MAX_NARR][IWAVE_NNEI]; /**< recv buffers, corresponding to recv neighbors */
  /** 
   * Assign ranks flag.
   * - 0: make a new communicator to which topology information has been attached, and determin process coords
   *      in cartesian topology by calling MPI_Cart_coords,
   * - 1: read from a file the process coord in cartesian topology. In this case, ccomm = MPI_COMM_WORLD.
   * 
   * Refer to \ref createparallel.
   */
} PARALLELINFO;


/*----------------------------------------------------------------------------*/
/**
 * Stores MPI and OMP data for global access.
 *
 * @param [in] ts - (int) thread support level
 */
int initparallel(int ts);

/*----------------------------------------------------------------------------*/
/**
 * Stores global MPI and OMP data for global access.
 *
 * @param [in] ts - (int) thread support level
 */
void initparallel_global(int);

/*----------------------------------------------------------------------------*/
/**
 * Stores local (group, remote) MPI and OMP data for global access.
 *
 * reads cartesian grid info
 *
 * @param [in] ts - (int) thread support level
 */
int initparallel_local(PARARRAY, FILE *);

/**
 * Initializes those parts of the parallel info struct which are indep 
 * of phys grid dimn - depend only on cart comm
 * 
 * Called in \ref iwave_construct. 
 *
 * @param [out] pinfo - (PARALLELINFO *) parallel info struct
 * @param [in] stream - (FILE *) file pointer where to dump run-time message
 */
int initpinfo(PARALLELINFO *pinfo, FILE *stream);
/*----------------------------------------------------------------------------*/
/*
  Creates MPI Cartesian grid and computes neighbor processor ranks.
  Exits quitely if process is out of Cartesian grid.
  Sets OMP threads.
*/
/*
 * Determines the number of processes in each dimension, i.e.\ cdims and asgnrk 
 * in \ref PARALLELINFO. 
 */
/* 10.04.11: deprecated WWS
   int readparallel(PARALLELINFO *pinfo, PARARRAY *pars, FILE *stream);
*/
/*----------------------------------------------------------------------------*/
/*
 * Computes and stores neighbor processor ranks for correct physical grid dimn
 * Exits quietly if process is out of Cartesian grid and sets OMP threads.
 *
 * @param [out] pinfo - (PARALLELINFO *) parallel info pointer
 * @param [in] ndim - (int) physical grid dimension (as opposed to MPI cart comm dim)
 * @param [in] stream - (FILE *) file pointer where to dump run-time message
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int initexch(PARALLELINFO *pinfo, int ndim, FILE *stream);

/*----------------------------------------------------------------------------*/
/*
 * Destroy Cartesian grid and MPI-related.
 *
 * int fabort :  flag to abort.
 *
 * Destroy the parallel info.
 *
 * Frees the communicator and destroys the send and recv info. No STORAGE 
 * DEALLOCATION occurs at this point.
 */
/* renamed 10.04.11 WWS
   int destroyparallel(PARALLELINFO *pinfo, int fabort);
*/
/** resets data exchange info
 * 
 * @param [out] pinfo - (PARALLELINFO *) parallel info pointer
 * @param [in] fabort - (int) physical grid dimension (as opposed to MPI cart comm dim)
 */
int destroypinfo(PARALLELINFO *pinfo, int fabort);
/*----------------------------------------------------------------------------*/

/** 
 * Wave simulation struct. It includes:
 * - an IMODEL struct storing the field variables,
 * - an PARALLELINFO struct storing parallel info for MPI.
 */
typedef struct s_IWAVE {
  /** Parallel info struct */
  PARALLELINFO pinfo; 
  /** Model struct */
  IMODEL model;
  /** 
   * Calculate statistics. 
   * stats > 0 output the timing data after the end of time loop. 
   * Its parameter name is "stats" defined in ../src/state/iwavefun.c
   */
  int stats;
  /** 
   * No post-time-step function, e.g.\ applying boundary conditions. 
   * - 0: apply boundary conditions,
   * - 1: don't apply boundary conditions.
   *
   * Its parameter name is "nopts" defined in ../src/state/iwavefun.c
   * Refer to \ref iwave_run.
   */
  int nopts;
  /**
   * Level of outputing actions.
   * - 0: no output for computing and exchange actions,
   * - 1: print time index for computing actions: 
   * - 2,3,4,5: output for both computing and exchange actions:
   *    - computing: print time index and and the internal step index,
   *    - exchange: print which variable array is exchanged in an action.
   * - >5: also dump all arrays after every time step to a file.
   *
   * Its parameter name is "printact" defined in src/state/iwavefun.c
   */
  int printact;
#ifdef IWAVE_USE_MPI
  /* ALL OF THIS STUFF REMOVED TEMP 13.03.12 WWS */
  /** 
   * Starting time of a simulation.
   *
   * Recorded in \ref prepexch which is called in \ref iwave_construct.
   */
  /*    double wt0; */
  /**
   * The minimum, maximum and average cpu time consumed on a 
   * certain action of a given variable.
   *
   * Note: actually the average cpu time = stat_times[ir][ia][2] / stat_calls[ir][ia][2].
   */
  /*    double stat_times[RDOM_MAX_NARR][MODEL_MAX_NACTION][3];  min, max, ave */
  /**
   * Indicates the place (step) where the minimum and maximum cpu time consumed 
   * on a certain action of a given variable happpens, and the number of such action 
   * applied in the entire simulation.
   */
  /*    int stat_calls[RDOM_MAX_NARR][MODEL_MAX_NACTION][3]; */
#endif
} IWAVE;

/* MAJOR IWAVE INTERFACES */
/** constructor */
/** 
 * IWAVE struct constructor.  Does the following operations: -
 * parallel info (\ref PARALLELINFO) related operation: -
 * initializes the parallel info struct, references \ref initpinfo -
 * determines the number of processes in each dimension and asgnrk
 * in \ref PARALLELINFO, references \ref readparallel - creates MPI
 * Cartesian grid and computes neighbor processor ranks - model
 * struct (IMODEL) related operation: - sets model struct defaults,
 * references \ref im_construct - initializes model struct,
 * references \ref minit - read other misc flags of an IWAVE struct,
 * references \ref readparpars - create model struct, references
 * IMODEL::mcrea: - read grid info - make action list - create
 * stencil - compute size of the domain on its host processor
 * according to the parallel info Cartesian grid dimensions (cdims)
 * and Cartesian rank (crank) - allocated memory, NOTE: STORAGE
 * ALLOCATION occurs only once, that is for \ref IMODEL::ld_a and
 * set the dimension information for IMODEL::ld_r[i],
 * IMODEL::ld_s[i], but let \n IMODEL::ld_r[i]._s[j]._s =
 * IMODEL::ld_s[i]._s[j]._s = IMODEL::ld_a[i]._s[j]._s,\n where the
 * first '_s' is RARR type and the second '_s' is (ireal *) type -
 * etc - set up paralle information for communication - set the
 * MPI_DATATYPE for sending and receiving: - let\n
 * PARALLELINFO::reinfo[ia][i].buf = IMODEL::ld_r[i]._s[ia]._s (two
 * pointers)\n PARALLELINFO::seinfo[ia][i].buf =
 * IMODEL::ld_s[i]._s[ia]._s (two pointers) - constructor the
 * EXCHANGEINFO::type, which the vector datatype striding the
 * continuous memory to the pieces for communication (send and
 * receive)
 *
 * NOTE: as mentioned before, STORAGE ALLOCATION occurs only once,
 * that is for \ref IMODEL::ld_a. All other domains in IMODEL struct
 * is of different dimension information, but refer to the same
 * chunk of physical memories of (ireal *) type with IMODEL::ld_a.
 *
 * @param [out] state - (IWAVE *) IWAVE struct pointer
 * @param [in] pars - (PARARRAY *) parameter array pointer 
 * @param [in] stream - (FILE *) file pointer for dumping run-time message and warning
 * @param [in] minit - a function pointer that initializes an IMODEL struct. References 
 *                     \ref IMODEL::minit
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int iwave_construct(IWAVE * state, PARARRAY * pars, FILE * stream, IWaveInfo const & ic);

/** write */
/** 
 * Dump information of all arrays including the allocated, computational virtual, send 
 * and recv ones as well as parallel info. 
 */
void iwave_printf(IWAVE * state, PARARRAY * pars, FILE * stream);

/** assigns all dynamical fields to zero initial values
    can use same function for fwd, lin, adj
*/
/**
 * Resets all dynamical fields zeros, time index = it and internal 
 * step index = 0. 
 */
void iwave_dynamic_init(IWAVE * state,
			int it, IWaveInfo const & ic);

/**
 * STORAGE DEALLOCATION for IMODEL struct and frees the communicator and destroys the send 
 * and recv info.
 * 
 * @param state - (IWAVE *) target IWAVE struct
 */
void iwave_destroy(IWAVE * state, FD_MODELDEST d);

/* MORE IWAVE/PARALLELINFO HELPER FUNCTIONS */
/** 
 * Normal exit. 
 * should only be called after everything is done, so 
 * assume stream exists.
 */
void quietexit(PARARRAY * pars,FILE ** stream);

/**
 * Cleans up environment and exits.
 */
void abortexit(int err,PARARRAY * pars,FILE ** stream);

/**
 * Initializes parallel output stream.
 * Called in driver.
 * @param [out] stream - (FILE **) a pointer of FILE pointer
 * @param [in] rk - (int) rank in MPI_COMM_WORLD
 * @param [in] sz - (int) number of processes in the group of MPI_COMM_WORLD
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */

int initoutstream(FILE ** stream, int rk, int sz);
/**
 * Initializes input parameter structures and populates them.
 * Called in driver.
 */

int readinput(PARARRAY ** pars, FILE * stream, int argc, char **argv);
/** 
 * Reads miscellaneous flags from param array, such as stats, nopts and printact.
 * Called in \ref iwave_construct.
 */

void readparpars(PARARRAY * pars,
		 int * stats,
		 int * nopts,
		 int * printact,
		 FILE * out);

/** 
 * Static storage of MPI comm params.
 */
void storeparallel(PARALLELINFO *pinfo);

/**
 * Sets up parallel information for communication - sets the MPI_DATATYPE for sending and receiving.
 * 
 * Calls \ref rd_setexchangeinfo.
 */
int setrecvexchange(IMODEL * model, PARALLELINFO * pinfo, FILE * stream, IWaveInfo const & ic);

/**
 * Dumps parallel info.
 */
int dump_pi(PARARRAY * pars, PARALLELINFO *pinfo, FILE *stream);  

/**
 * Dumps information of all arrays including the allocated, and computational virtual ones.
 */

void dump_ac(PARARRAY * pars, IMODEL * model, FILE * stream);

/**
 * Dumps information of all arrays including the send and recv virtual ones.
 */
void dump_rs(PARARRAY * pars, IMODEL * model, FILE * stream, int sends);

#ifdef IWAVE_USE_MPI
/**
 * Prepares data exchanges: commit the MPI_DATATYPE.
 */
void prepexch(PARALLELINFO * pinfo, IMODEL * model, IWaveInfo const & ic);
/* REMOVED TEMP 13.03.12 WWS
//		int stats, 
//		double * wt0); */

/**
 * Prepare timing data at start of time loop. 
 */
/* REMOVED TEMP 13.03.12 WWS
//  void preptime(double stat_times[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
//		int stat_calls[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
//                IWAVE * state); */

/**
 * Outputs timing data after the end of time loop.
 * Refer to \ref IWAVE. Called in driver.
 */
/* REMOVED TEMP 13.03.12 WWS
//  void printtime(double stat_times[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
//                 int stat_calls[RDOM_MAX_NARR][MODEL_MAX_NACTION][3],
//                 IWAVE * state, FILE * stream); */

#endif

#endif /*__IWAVE_H_*/

