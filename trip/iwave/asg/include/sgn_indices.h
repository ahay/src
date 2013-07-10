/* 
combines
sg.h
Igor Terentyev.
and
 
sgn.h
Igor Terentyev.
*/
/*============================================================================*/

#ifndef __SGN_INDICES_H_
#define __SGN_INDICES_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "stencil.h"
#include "rdomain.h"

/*----------------------------------------------------------------------------*/
/*
Indices of arrays - pure acoustic. Permutation of {0,...,3/5/7}.
D_P0  :  pressure.
D_V0  :  velocity.
D_MPx :  pressure multiplier (bulk modulus).
D_MVx :  velocity multiplier (buoyancy).
*/
#define D_P0  0
#define D_MP0 1
#define D_V0  2
#define D_MV0 3
#define D_V1  4
#define D_MV1 5
#define D_V2  6
#define D_MV2 7

/*----------------------------------------------------------------------------*/
/*
Indices of additional arrays due to NPML. Permutation of {4/6/8,...,5/10/15}.
D_P  :  pressure (bar).
D_V  :  velocity (bar).
D_MP :  pressure multiplier (bulk modulus).
D_MV :  velocity multiplier (buoyancy).
D_EV :  eta velocity. Allocated as 1D arrays: Nx1x1.
D_EP :  eta pressure. Allocated as 1D arrays: Nx1x1.
 
All pressure arrays must have same global coordinates - P0 is used to get
global coordinates and sizes of all pressure arrays.
*/
#define D_EP0 (2 + 2 * RARR_MAX_NDIM)  /* 8  */
#define D_EV0 (D_EP0 + 1)              /* 9  */
#define D_P1  (D_EV0 + 1)              /* 10 */
#define D_EP1 (D_P1  + 1)              /* 11 */
#define D_EV1 (D_EP1 + 1)              /* 12 */
#define D_P2  (D_EV1 + 1)              /* 13 */
#define D_EP2 (D_P2  + 1)              /* 14 */
#define D_EV2 (D_EP2 + 1)              /* 15 */

/*----------------------------------------------------------------------------*/
#if   RARR_MAX_NDIM == 1
    static const IPNT D_V  = {D_V0 };
    static const IPNT D_MV = {D_MV0};
#elif RARR_MAX_NDIM == 2
    static const IPNT D_V  = {D_V0 , D_V1 };
    static const IPNT D_MV = {D_MV0, D_MV1};
#elif RARR_MAX_NDIM >= 3
    static const IPNT D_V  = {D_V0 , D_V1 , D_V2 };
    static const IPNT D_MV = {D_MV0, D_MV1, D_MV2};
#endif

#if   RARR_MAX_NDIM == 1
    static const IPNT D_P  = {D_P0 };
    static const IPNT D_EP = {D_EP0};
    static const IPNT D_EV = {D_EV0};
#elif RARR_MAX_NDIM == 2
    static const IPNT D_P  = {D_P0 , D_P1 };
    static const IPNT D_EP = {D_EP0, D_EP1};
    static const IPNT D_EV = {D_EV0, D_EV1};
#elif RARR_MAX_NDIM >= 3
    static const IPNT D_P  = {D_P0 , D_P1 , D_P2 };
    static const IPNT D_EP = {D_EP0, D_EP1, D_EP2};
    static const IPNT D_EV = {D_EV0, D_EV1, D_EV2};
#endif
/*----------------------------------------------------------------------------*/

#endif
