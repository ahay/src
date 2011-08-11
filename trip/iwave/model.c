/* Model structure. */
/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* 
model.c
Igor Terentyev.
*/
/*============================================================================*/

#include <trip/base.h>
#include <trip/grid.h>
/*^*/

#include "model.h"

#include "rdomain.h"
/*^*/

#ifndef _sf_model_h

/** limits maximum number of actsions definable for any scheme */
#define MAX_ACTIONS 100
/*^*/

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
/*^*/

/*----------------------------------------------------------------------------*/

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
/*^*/

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
/*^*/



#endif

void init_step(TIMESTEPINDEX * ts, int it) 
/*<*
 * Initialize time structure.
 *
 * @param [out] ts - (TIMESTEPINDEX *) 
 * @param [in] it - (int) time index
>*/
{
  ts->it=it;
  ts->iv=0;
}

void next_step(TIMESTEPINDEX * ts, TIMESTEPINFO tsinfo) 
/*<*
 * Advance simulation time.
 * 
 * If finishing all actions in a time step, increase time index
 * <em>it</em> and reset internal step index <em>iv</em> zero,
 * otherwise increase <em>iv</em>.  @param [out] ts - (TIMESTEPINDEX
 * *) @param [in] tsinfo - (TIMESTEPINFO)
 >*/
{
  if (ts->iv == tsinfo.npair-1) {
    ts->iv=0;
    ts->it++;
  }
  else {
    ts->iv++;
  }
}

int less_than(TIMESTEPINDEX t1, TIMESTEPINDEX t2) 
/*< *
 * Compare two time structures.
 *
 * @return 1, if time index of t1 is less than that of t2, or the two
 * time indices equal and the internal step index of t1 is less than
 * that of t2. Otherwise, return 0.
 >*/
{
  if ((t1.it<t2.it) || ((t1.it == t2.it) && (t1.iv<t2.iv))) return 1;
  return 0;
}

ireal get_time(TIMESTEPINDEX ts) 
/*< Get current simulation time. >*/
{
  return ts.it*ts.dt;
}

#ifndef _sf_model_h

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
/*^*/

#endif

/*int im_setnull(IMODEL *model)*/
int im_construct(IMODEL * model) 
/*<  Sets model struct defaults. >*/
{

  model->specs=NULL;

  /* WWS 22.02.08 removed zero of ndim, added null constr for grid */  
  init_default_grid(&(model->g));

  model->nnei = 0;
  //  sten_setnull(&(model->sten));

  model->tsinfo.narr = 0;
  rd_a_setnull(&(model->ld_a));
  rd_a_setnull(&(model->ld_c));
  rd_a_setnull(&(model->ld_p));
  
  model->ld_s = model->ld_r = NULL;

  init_step(&(model->tsind),INT_MAX);

  return 0;
}

/*----------------------------------------------------------------------------*/

int im_setndim(IMODEL *model)
/*< Allocates the domain array ld_s, ld_r and set number of neighbor processors (nnei). >*/
{
    int nnei, i, err;
    RDOM *p;

    err = gen_3n1(model->g.dim, &nnei);
    if ( err ) return err;

    p = (RDOM*)malloc(nnei * 2 * sizeof(RDOM));
    if ( p == NULL ) return E_ALLOC;
    
    for ( i = 0; i < nnei * 2; ++ i ) rd_a_setnull(p + i);

    if ( model->ld_s != NULL ) free(model->ld_s);
    model->ld_s = p;
    model->ld_r = model->ld_s + nnei;
    model->nnei = nnei;
    
    return 0;
}
/*----------------------------------------------------------------------------*/
