/* 
model.h
Igor Terentyev.
********************************************************************************
Model function types.
Model structure.
*/
/*============================================================================*/

/**
 * \file model.h
 *
 * Model structure and its related functions
 */
#ifndef __IWAVE_MODEL_H_
#define __IWAVE_MODEL_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "parser.h"
#include "stencil.h"
/* WWS 22.02.08 */
#include "grid.h"
#include "rdomain.h"

/** limits maximum number of actsions definable for any scheme */
#define MAX_ACTIONS 100

/*----------------------------------------------------------------------------*/
/**
 * Time step action pair. Each action pair couples a dynamic array,
 * identified by its index (arr) in the model RDOMAIN, to an action
 * (action). The possible actions are identified by cpp definitions in
 * \ref model_action.h .  These pairs correspond to updates of each
 * array in the timestep loop.  For example, typical actions are
 * COMPUTE and EXCHANGE, signifying serial update code and data
 * exchange via MPI respectively. So the TIMESTEPINFOPAIR (10,COMPUTE)
 * would indicated that the array 10 in the RDOMAIN should be updated
 * by calling whatever timestep function does that job.
 */
typedef struct
{
  /** array index */
  int arr;
  /** action code */
  int action;
} TIMESTEPINFOPAIR;

/**
 * Time step information
 * 
 * Includes variable indicator set and all action pairs within a time step.
 */
typedef struct {
  /** 
   * number of dynamic arrays in a simulation (representing dynamic
   * arrays like pressures, velocities, stresses). These are the
   * arrays that are to be updated during the timestepping loop.
   * Must be less than RDOM_MAX_NARR, since these arrays are amongst
   * those represented in the RDOMAIN for each model instance.
   */
  int narr;
  
  /** 
   * RDOMAIN indices of the dynamic arrays. Only the part
   * arrs[0,narr-1] is used. arrs[i] for i in [0,narr-1] is the
   * index in the model RDOMAIN of the ith dynamic array. This
   * device decouples the ordering of domains in the RDOMAIN from
   * the list of dynamic arrays for which actions must be specified.
   */
  int arrs[RDOM_MAX_NARR];
  
  /**
   * number of action pairs, corresponding to distinct actions 
   * during timestep loop.
   */
  int npair;

  /** array of action pairs, defining simulation timestep */
  TIMESTEPINFOPAIR pairs[MAX_ACTIONS];

} TIMESTEPINFO;

/*----------------------------------------------------------------------------*/
/**
 * Time structure for use with iwave simulators.
 * 
 * Includes both time index and internal step flag. Associated functions 
 * advance and compare these and also return true time. \ref init_step, 
 * \ref next_step, \ref less_than and \ref get_time.
 */
typedef struct {
  int it;    /**< time index */
  int iv;    /**< internal step index */
  ireal dt;  /**< real time step */
} TIMESTEPINDEX;

/**
 * Initialize time structure.
 *
 * @param [out] ts - (TIMESTEPINDEX *) 
 * @param [in] it - (int) time index
 */
void init_step(TIMESTEPINDEX * ts, int it);

/**
 * Advance simulation time.
 * 
 * If finishing all actions in a time step, increase time index
 * <em>it</em> and reset internal step index <em>iv</em> zero,
 * otherwise increase <em>iv</em>.  @param [out] ts - (TIMESTEPINDEX
 * *) @param [in] tsinfo - (TIMESTEPINFO)
 */
void next_step(TIMESTEPINDEX * ts, TIMESTEPINFO tsinfo);

/**
 * Compare two time structures.
 *
 * @return 1, if time index of t1 is less than that of t2, or the two
 * time indices equal and the internal step index of t1 is less than
 * that of t2. Otherwise, return 0.
 */
int less_than(TIMESTEPINDEX t1, TIMESTEPINDEX t2);

/**
 * Get current simulation time.
 *
 * @return current simulation time.
 */
ireal get_time(TIMESTEPINDEX ts);
  
/*----------------------------------------------------------------------------*/
/*
int ndim          :  number of dimensions.
STENCIL sten      :  stencil;
void *specs       :  model specific data.
TIMESTEP_FUN ts   :  timestep function.
void *tspars      :  timestep function parameters.
RDOM ld_a, ld_c   :  allocated domain, computational virtual domain.
RDOM *ld_s, *ld_r :  recv, send virtual domains.
int nt            :  number of timesteps.

BIG TODO: remove ld_s, ld_r to particular models and get rid of domains;
          they should be special objects that are prepared before send:
          copy overlaped data before sending, etc.
*/
/* WWS 22.08.02: added g, dt, it, removed ndim, nt */
/* Igor, Feb 26: added physical domain ld_p
RDOM ld_p
*/
/**
 * Model - data members.
 *
 * Contains almost all data structures for a simulation except the 
 * parallel information. Also contains pointer to \ref FD_MODEL
 */
typedef struct IMODEL {
  /**
   * pointer to struct containing pointers to functions used to
   * implement model behaviour. Must be initialized by
   * FD_MODEL::fd_model_init to point to \ref FD_MODEL object.
   * Declared here as void* because FD depends on IMODEL.  This also
   * opens up the possibility that future non-FD models might use
   * different collections of functions.
   */
  void * specs;

  /** 
   * Grid struct - global primal simulation grid. All other structures
   * follow from this one. Usually read from file, using the grid
   * package i/o functions. If several grids are involved in a
   * simulation (as is the case for staggered grid schemes, for
   * example), this grid provides all information necessary to create
   * all other grids.
   */
  grid g;
  /** 
   * Local primal grid = the grid implicit in the local computational 
   * domain. For parallel simulation via domain decomp, a subgrid
   * of g.
   */
  grid gl;
  /** artificial boundary layer widths */
  IPNT nls,nrs;
  /** Number of neighbor processors */
  int nnei;
  /** Time step information */
  TIMESTEPINDEX tsind;
  /** Time structure */
  TIMESTEPINFO tsinfo;

  /** Allocated domain */
  RDOM ld_a;
  /** Computational virtual domain */
  RDOM ld_c;
  /** Physical virtual domain - may be same as computational virtual domain */
  RDOM ld_p;
  /** Array of send virtual domains */
  RDOM *ld_s;
  /** Array of recv virtual domains */
  RDOM *ld_r;
  
} IMODEL;

/*----------------------------------------------------------------------------*/

typedef int (*IMODELINIT_FUN)(PARARRAY *pars, FILE *stream, IMODEL *model);
/*----------------------------------------------------------------------------*/
/**
 *  Sets model struct defaults.
 */
int im_construct(IMODEL *model);
/*----------------------------------------------------------------------------*/
/**
 * Destroys model (STORAGE DEALLOCATION).
 */
int im_destroy(IMODEL *model);
/*----------------------------------------------------------------------------*/
/**
 * Allocates the domain array ld_s, ld_r and set number of neighbor processors (nnei).
 */
int im_setndim(IMODEL *model);
/*----------------------------------------------------------------------------*/

#endif /*__MODEL_H_*/
