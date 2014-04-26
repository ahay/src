/* 
model.h
Igor Terentyev.
********************************************************************************
Model function types.
Model structure.
*/
/*============================================================================*/

/*
 * Modified 03.12 by WWS - case switch in action list changes from
 * array to be updated to internal time step, to accommodate the many
 * cases where more than one array is updated, hence needs to be
 * exchanged, and perhaps different numbers in fwd, lin, and adj
 * (inversion apps)
 */
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
  int niv;   /**< number of internal steps */
  ireal dt;  /**< real time step */
} TIMESTEPINDEX;

/**
 * Advance simulation time.
 * 
 * If finishing all actions in a time step, increase time index
 * <em>it</em> and reset internal step index <em>iv</em> zero,
 * otherwise increase <em>iv</em>.  
 * @param [out] ts - (TIMESTEPINDEX) 
 */
void next_step(TIMESTEPINDEX * ts);

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
   * pointer to struct containing FD scheme parameters. Must be initialized by
   * IWaveInfo::fd_modelinit
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
  /** artificial boundary layer widths - THESE BELONG ELSEWHERE!! */
  IPNT nls,nrs;
  /** Number of neighbor processors */
  int nnei;
  /** Time step information */
  TIMESTEPINDEX tsind;

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
int im_destroy(IMODEL *model, void (*destr)(void **));
/*----------------------------------------------------------------------------*/
/**
 * Allocates the domain array ld_s, ld_r and set number of neighbor processors (nnei).
 */
int im_setndim(IMODEL *model);
/*----------------------------------------------------------------------------*/

#endif /*__MODEL_H_*/
