/*
  Author: XW 01-19-2011
  wangxin.tom@gmail.com
  esgn_indices.h 
  mathematical problems: isotropic elastic wave equations (1st order velocity stress)
  methods: staggered grid
*/
/*============================================================================*/

#ifndef __ESGN_INDICES_H_
#define __ESGN_INDICES_H_
/*----------------------------------------------------------------------------*/

#include "utils.h"
#include "stencil.h"
#include "rdomain.h"

/*----------------------------------------------------------------------------*/
/*
Indices of arrays - isotropic elastic.
D_P0  :  pressure or normal stress component along axis-0 direction.
D_V0  :  velocity.
D_S0  :  shear stress component in axis-0, axis-1 plane
D_MPij:  pressure or normal stress multiplier (bulk modulus or Lame parameter).
D_MSij: shear stress component multiplier (Lame second parameter)
D_MVx :  velocity multiplier (buoyancy).
*/
/*< pz */
#define D_P0   0
#define D_MP00 1                 /* lambda + 2 mu */
#define D_MP01 2                 /* lambda */
#define D_MP02 D_MP01

/*< vz */
#define D_V0   3
#define D_MV0  4

/*< szx */
#define D_S0   5
#define D_MS0  6 

/*< px */
#define D_P1   7
#define D_MP10 D_MP01
#define D_MP11 D_MP00
#define D_MP12 D_MP01

/*< vx */
#define D_V1  8
#define D_MV1 9

/*< sxy */
#define D_S1   10
#define D_MS1  11

/*< py */
#define D_P2   12
#define D_MP20 D_MP01
#define D_MP21 D_MP01
#define D_MP22 D_MP00

/*< vy */
#define D_V2   13
#define D_MV2  14

/*< szy */
#define D_S2   15
#define D_MS2  16

/*----------------------------------------------------------------------------*/
/*
coefficients (ETA) in PML
D_EP      :  PML coefficient for pressure or normal stresses
D_EV      :  PML coefficient for velocity
D_ES      :  PML coefficient for shear stresses

D_E_PRIME :  PML coefficient for prime grids (z-axis, x-axis, y-axis)
D_E_DUAL  :  PML coefficient for dual grids  (z-axis, x-axis, y-axis)
 
All pressure arrays must have same global coordinates - P0 is used to get
global coordinates and sizes of all pressure arrays.
*/
#define D_E_PRIME0 17
#define D_E_DUAL0  18
#define D_E_PRIME1 19
#define D_E_DUAL1  20
#define D_E_PRIME2 21
#define D_E_DUAL2  22

/*----------------------------------------------------------------------------*/
#if RARR_MAX_NDIM == 2
static const IPNT D_V  = {D_V0 , D_V1 };
static const IPNT D_MV = {D_MV0, D_MV1};
static const IPNT D_EV = {D_E_DUAL0, D_E_DUAL1};
#elif RARR_MAX_NDIM >= 3
static const IPNT D_V  = {D_V0 , D_V1 , D_V2 };
static const IPNT D_MV = {D_MV0, D_MV1, D_MV2};
static const IPNT D_EV = {D_E_DUAL0, D_E_DUAL1, D_E_DUAL2};
#endif

#if RARR_MAX_NDIM == 2
static const IPNT D_P  = {D_P0 , D_P1 };
static const IPNT D_MP = {D_MP00, D_MP01};
static const IPNT D_EP = {D_E_PRIME0, D_E_PRIME1};
#elif RARR_MAX_NDIM >= 3
static const IPNT D_P  = {D_P0 , D_P1 , D_P2 };
static const IPNT D_MP = {D_MP00, D_MP01, -1};
static const IPNT D_EP = {D_E_PRIME0, D_E_PRIME1, D_E_PRIME2};
#endif

#if RARR_MAX_NDIM == 2
static const IPNT D_S  = {D_S0 , D_S1 };
static const IPNT D_MS = {D_MS0, D_MS1};
static const IPNT D_ES = {D_E_DUAL0, D_E_DUAL1};
#elif RARR_MAX_NDIM >= 3
static const IPNT D_S  = {D_S0 , D_S1 , D_S2 };
static const IPNT D_MS = {D_MS0, D_MS1, D_MS2};
static const IPNT D_ES = {D_E_DUAL0, D_E_DUAL1, D_E_DUAL2};
#endif
/*----------------------------------------------------------------------------*/

#endif
