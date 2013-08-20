/*
  Author: XW 01-19-2011
  wangxin.tom@gmail.com
  esgn.h 
  mathematical problems: isotropic elastic wave equations (1st order velocity stress)
  methods: staggered grid
*/
/*===========================================================================*/

#ifndef __ESGN_H_
#define __ESGN_H_

/**
 * \file esg.h
 * Staggered-grid finite difference solver for isotropic elastic wave equations 
 * with PML
 *
 * <p>Note: in order to run this solver, turn off the option IWAVE_USE_FMGR 
 * in grid/include/gridio.h <- "this is fixed" </p>
 *
 * Possible field variables indices:
 * - D_P0  : pressure or normal stress component along axis-0 direction
 * - D_V0  : velocity
 * - D_S0  : shear stress component in axis-0, axis-1 plane
 * - D_MPij: pressure or normal stress multiplier (e.g., bulk modulus or Lame parameters)
 * - D_MSij : shear stress component multiplier (e.g.,Lame second parameter)
 * - D_MVx : velocity multiplier  (buoyancy)
 * 
 * e.g., 3D stagger grid for elastic media
 * <pre>
 *
 *                   NAN_____________ syz  
 *                   / .             /|
 *                /               /   |              
 *             /       .       /      |
 *         sxz--------------vz        | 
 *          |          .    |         |
 *          |        sxy  . |   .   . |
 *          |        .      |        /vy
 *          |     .         |     /            
 *          |  .            |  /   
 *          vx--------------px (py, pz)
 *    (i-1/2,j,k)        (i,j,k)
 *
 * </pre>
 */
/*============================================================================*/

#include "fd.h"
#include "esgn_read.h"
#include "create_sten_2k.h"
#include "esgsteps.h"
#include "model_action.h"
#include "defaults.h"

/*----------------------------------------------------------------------------
 * Parameters for time step function (elastic wave + PML).
 *----------------------------------------------------------------------------
 *           _____________________________________
 *          |    |         III               |    |
 *          |____|___________________________|____|
 *          |    |                           |    |
 *          |    |                           |    |
 *          |    |        physical           |    |
 *          | I  |         domain            | II |
 *          |    |                           |    |
 *          |    |                           |    |
 *          |    |                           |    |
 *          |    |                           |    |
 *          |____|___________________________|____|
 *          |    |         IV                |    |
 *          |____|___________________________|____|
 *
 * e.g., region I and III are overlaped at the upper-left corner for no special 
 *       reasons
 */
typedef struct {
  ireal dt;         /**< time step size */
  RPNT lam;         /**< courant params: dt/dx [dt/dy dt/dz] */
  int k;            /**< scheme order: 2*k */
  int ndim;         /**< dimension, copied from IMODEL.grid */
  RDOM *ld_pml;     /**< auxillary variables for PML (z,x,y axis then PML region I, II, ...), 
                       memory allocated in esg_readschemeinfo -> fd_mread -> iwave_static_init 
                    */
} ESGN_TS_PARS;

void esgn_ts_parcopy(void * tgt, const void * src);

/*----------------------------------------------------------------------------
 * FUNCTIONS MANDATED BY FD_MODEL
 *----------------------------------------------------------------------------
 */
int esg_getndim();
int esg_getnarr();
//int esg_getindices(int);
//const char * getnames(int i);
int esg_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
int esg_set_grid_type(FILE *, int, IPNT[RDOM_MAX_NARR]);
int esg_build_sten_dep(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
int esg_create_sten(FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *, void*);
int esg_assign_action_array(FILE *, IMODEL *);
int esg_alter_dom(int, IPNT, IPNT);
const char * asg_ind2str(int);  
int esg_modelpostts(int, IMODEL *);
int esg_modelinit(PARARRAY *pars, FILE *stream, IMODEL *model);
int esg_modeldest(IMODEL * model);
int esg_dynamic_init(IMODEL * model);

int esgn_gts2d_24(RDOM*, RDOM*, RDOM*, int, void *, int);
int esgn_gts3d_24(RDOM*, RDOM*, RDOM*, int, void *, int);
int esgn_gts2d_210(RDOM*, RDOM*, RDOM*, int, void *, int);
int esgn_gts3d_210(RDOM*, RDOM*, RDOM*, int, void *, int);
#endif /*__ESGN_H_*/
