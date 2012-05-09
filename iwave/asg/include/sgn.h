/* 
sgn.h
Igor Terentyev.
********************************************************************************
Acoustics. NPML. Staggered grid. 1D, 2D, 3D. 
Implementation with only bar variables.
Arrays:
    velocities, multipliers (from sg.h),
    pressures, etas.
Stencils:
    2-2, 2-4 standard.
*/
/*============================================================================*/

#ifndef __SGN_H_
#define __SGN_H_
/*----------------------------------------------------------------------------*/

#include "fd.h"
#include "sgn_read.h"
#include "create_sten_2k.h"
#include "model_action.h"
#include "defaults.h"

/*----------------------------------------------------------------------------
 * Parameters for time step function.
 *----------------------------------------------------------------------------
 */
typedef struct {
  ireal dt;      // time step - copied from IMODEL.tsinfo
  RPNT lam;      // courant params
  int k;         // scheme order
  int ndim;      // dimension, copied from IMODEL.grid
  int psingle;   // use interior scheme - currently inactive
  IPNT eflags;   // disused PML vs. acoustic flags
  IPNT lbc;      // flag left boundary conditions
  IPNT rbc;      // flag right boundary conditions
} SGN_TS_PARS;  

void sgn_ts_parcopy(void * tgt, const void * src);

/*----------------------------------------------------------------------------
 * FUNCTIONS MANDATED BY FD_MODEL
 *----------------------------------------------------------------------------
 */
//int asg_readgrid(PARARRAY *, FILE *, IMODEL *);
//int asg_readtimegrid(PARARRAY *, FILE *, IMODEL *);
int asg_readschemeinfo(PARARRAY *, FILE *, IMODEL *);
int asg_set_grid_type(FILE *, int, IPNT[RDOM_MAX_NARR]);
int asg_build_sten_dep(FILE *, int, int[RDOM_MAX_NARR][RDOM_MAX_NARR]);
//int asg_create_sten(FILE *, int, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *, void*);
int asg_create_sten(FILE *, int, IPNT[RDOM_MAX_NARR], int[RDOM_MAX_NARR][RDOM_MAX_NARR], STENCIL *, void*);
int asg_assign_action_array(FILE *, IMODEL *);
int asg_alter_dom(int, IPNT, IPNT);
//int asg_readmedia(PARARRAY *, FILE *, IMODEL *);
const char * asg_ind2str(int);  
int asg_noop(int, IMODEL *);
int asg_modelpostts(int, IMODEL *);
int asg_modelinit(PARARRAY *pars, FILE *stream, IMODEL *model);
int asg_modeldest(IMODEL * model);
void asg_refsubstep(int*,int*,const IMODEL *);

/*----------------------------------------------------------------------------
 * IGOR'S TIMESTEP CODE
 *----------------------------------------------------------------------------
 */
int sgn_ts1d_22(RDOM*, int, void *);
int sgn_ts2d_22(RDOM*, int, void *);
int sgn_ts3d_22(RDOM*, int, void *);
int sgn_ts1d_24(RDOM*, int, void *);
int sgn_gts2d_24(RDOM*, int, void *);
int sgn_ts3d_24(RDOM*, int, void *);
int sgn_ts1d_210(RDOM*, int, void *);
int sgn_ts3d_210(RDOM*, int, void *);
int sgn_ts2d_2K(RDOM*, int, void *);
int sgn_ts3d_2K(RDOM*, int, void *);

// Bill's timestep code
int duhasg24_2d(RDOM*, int, void *);
int duhasgfm24_2d(RDOM*, RDOM*, int, void *);
int duhasgam24_2d(RDOM*, RDOM*, int, void *);

// Dong's timestep code
int asg_fts2d_24(RDOM *, int, void *);
int asg_ftsm2d_24(RDOM *, RDOM *, int, void *);
int asg_atsm2d_24(RDOM *, RDOM *, int, void *);
int asg_atsm2d_24p01(RDOM *, RDOM *, void *);
int asg_atsm2d_24v0(RDOM *, RDOM *, void *);
int asg_atsm2d_24v1(RDOM *, RDOM *, void *);


#endif /*__SGN_H_*/
