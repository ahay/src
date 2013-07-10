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
#include "sgnpars.h"
#include "sgn_read.h"
#include "create_sten_2k.h"
#include "defaults.h"

/* this declaration must be publicly available to external drivers:
   modelinit is passed as function pointer argument to iwave_construct. 
   On return, all other FD_MODEL function pointers are initialized and
   point to functions defined in this package - this is the only 
   package function which must be used before being assigned to an 
   FD_MODEL data member.
*/
int asg_modelinit(PARARRAY *pars, FILE *stream, IMODEL *model);

#endif /*__SGN_H_*/
