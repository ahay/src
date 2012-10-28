/**
   time step functions
   Author: Xin Wang 
   wangxin.tom@gmail.com 
   time: 09/24/2011  
   2-4 stagger finite difference in 3D for isotropic elastic wave
*/

#include "esgsteps.h"
#include "esgn.h"

/*============================================================================*/

#define NSTORE 20000
/* static ireal _delta[NSTORE]; */

#ifdef VEC
#undef VEC
#endif

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/
/*
1/24 constant.
*/
#define C24 4.166666666666666666666667e-2
/*----------------------------------------------------------------------------*/
/**
 * Normal stress updater of 3D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24ns (RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Shear stress szx updater of 3D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24ss0(RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Shear stress sxy updater of 3D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24ss1(RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Shear stress szy updater of 3D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24ss2(RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Velocity component vz updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24v0 (RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Velocity component vx updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24v1 (RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Velocity component vy updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts3d_24v2 (RDOM *, RDOM *, RDOM *, void *, int);


int esgn_gts3d_24(RDOM *dom, RDOM *rdom, RDOM *cdom, int iarr, void *pars,int _fwd) {

  ESGN_TS_PARS * esgnpars = (ESGN_TS_PARS *) pars;  
  /* sanity test - storage mode */
#ifdef STORE
  if ((dom->_s)[D_P0]._dims[0].n > NSTORE) return E_NOTIMESTEP;
#endif
  if ( iarr == D_P1 || iarr == D_P2)  return 0;
  if ( iarr == D_P0 ) return esgn_gts3d_24ns (dom, rdom, cdom, esgnpars, _fwd);
  if ( iarr == D_S0 ) return esgn_gts3d_24ss0(dom, rdom ,cdom, esgnpars, _fwd);
  if ( iarr == D_S1 ) return esgn_gts3d_24ss1(dom, rdom ,cdom, esgnpars, _fwd);
  if ( iarr == D_S2 ) return esgn_gts3d_24ss2(dom, rdom ,cdom, esgnpars, _fwd);
  if ( iarr == D_V0 ) return esgn_gts3d_24v0 (dom, rdom, cdom, esgnpars, _fwd);
  if ( iarr == D_V1 ) return esgn_gts3d_24v1( dom, rdom, cdom, esgnpars, _fwd);
  if ( iarr == D_V2 ) return esgn_gts3d_24v2( dom, rdom, cdom, esgnpars, _fwd);
    
  return E_NOTIMESTEP;
}

/*----------------------------------------------------------------------------*/
#ifdef VEC

int esgn_gts3d_24ns (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  
  fprintf(stderr, "Warning: VEC version 3D normal stress updater not implemented\n");
  return 0;
}

int esgn_gts3d_24ss0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {

  fprintf(stderr, "Warning: VEC version 3D shear stress s0 updater not implemented\n");
  return 0;
}

int esgn_gts3d_24ss1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {

  fprintf(stderr, "Warning: VEC version 3D shear stress s1 updater not implemented\n");
  return 0;
}

int esgn_gts3d_24ss2(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {

  fprintf(stderr, "Warning: VEC version 3D shear stress s2 updater not implemented\n");
  return 0;
}

int esgn_gts3d_24v0 (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  
  fprintf(stderr, "Warning: VEC version 3D v0 updater not implemented\n");
  return 0;
}

int esgn_gts3d_24v1 (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  fprintf(stderr, "Warning: VEC version 3D v1 updater not implemented\n");
  return 0;
}

int esgn_gts3d_24v2 (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  fprintf(stderr, "Warning: VEC version 3D v2 updater not implemented\n");
  return 0;
}

#else
/*----------------------------------------------------------------------------*/
/** normal stresses */
int esgn_gts3d_24ns(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, px_a, py_a, pz_a, mp00_a, mp01_a, vx_a, vy_a, vz_a;
  int px_aa, py_aa, pz_aa, mp00_aa, mp01_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
  int px_pml_I_a,   py_pml_I_a,   pz_pml_I_a,   px_pml_I_aa,   py_pml_I_aa,   pz_pml_I_aa;
  int px_pml_II_a,  py_pml_II_a,  pz_pml_II_a,  px_pml_II_aa,  py_pml_II_aa,  pz_pml_II_aa;
  int px_pml_III_a, py_pml_III_a, pz_pml_III_a, px_pml_III_aa, py_pml_III_aa, pz_pml_III_aa; 
  int px_pml_IV_a,  py_pml_IV_a,  pz_pml_IV_a,  px_pml_IV_aa,  py_pml_IV_aa,  pz_pml_IV_aa;
  int px_pml_V_a,   py_pml_V_a,   pz_pml_V_a,   px_pml_V_aa,   py_pml_V_aa,   pz_pml_V_aa;
  int px_pml_VI_a,  py_pml_VI_a,  pz_pml_VI_a,  px_pml_VI_aa,  py_pml_VI_aa,  pz_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend;
  register ireal * restrict _px_x_I,   * restrict _px_y_I,   * restrict _px_z_I; 
  register ireal * restrict _py_x_I,   * restrict _py_y_I,   * restrict _py_z_I;
  register ireal * restrict _pz_x_I,   * restrict _pz_y_I,   * restrict _pz_z_I;
  register ireal * restrict _px_x_II,  * restrict _px_y_II,  * restrict _px_z_II;
  register ireal * restrict _py_x_II,  * restrict _py_y_II,  * restrict _py_z_II;
  register ireal * restrict _pz_x_II,  * restrict _pz_y_II,  * restrict _pz_z_II;
  register ireal * restrict _px_x_III, * restrict _px_y_III, * restrict _px_z_III;
  register ireal * restrict _py_x_III, * restrict _py_y_III, * restrict _py_z_III;
  register ireal * restrict _pz_x_III, * restrict _pz_y_III, * restrict _pz_z_III;
  register ireal * restrict _px_x_IV,  * restrict _px_y_IV,  * restrict _px_z_IV;
  register ireal * restrict _py_x_IV,  * restrict _py_y_IV,  * restrict _py_z_IV;
  register ireal * restrict _pz_x_IV,  * restrict _pz_y_IV,  * restrict _pz_z_IV;
  register ireal * restrict _px_x_V,   * restrict _px_y_V,   * restrict _px_z_V;
  register ireal * restrict _py_x_V,   * restrict _py_y_V,   * restrict _py_z_V;
  register ireal * restrict _pz_x_V,   * restrict _pz_y_V,   * restrict _pz_z_V;
  register ireal * restrict _px_x_VI,  * restrict _px_y_VI,  * restrict _px_z_VI;
  register ireal * restrict _py_x_VI,  * restrict _py_y_VI,  * restrict _py_z_VI;
  register ireal * restrict _pz_x_VI,  * restrict _pz_y_VI,  * restrict _pz_z_VI;
  register ireal * restrict _mp00, * restrict _mp01;
  register ireal * restrict _vx3;
  register ireal * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
  register ireal * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  register ireal * restrict _epx, * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, dt2, vx3, vx2, vx1, vx0, dfdx, dfdy, dfdz, etaxdt, etaydt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  nz = s[D_P0]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_P0]._dims[0].gs;
  gxe = s[D_P0]._dims[0].ge;
  gys = s[D_P0]._dims[1].gs;
  gye = s[D_P0]._dims[1].ge;
  gzs = s[D_P0]._dims[2].gs;
  gze = s[D_P0]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_P0,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_P0]._dims[2].gs;
    gze_pml_I = s_pml[D_P0]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_P0,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_P0]._dims[2].gs;
    gze_pml_II = s_pml[D_P0]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_P0,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_P0]._dims[1].gs;
    gye_pml_III = s_pml[D_P0]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_P0,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_P0]._dims[1].gs;
    gye_pml_IV = s_pml[D_P0]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_P0,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_P0]._dims[0].gs;
    gxe_pml_V = s_pml[D_P0]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_P0,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_P0]._dims[0].gs;
    gxe_pml_VI = s_pml[D_P0]._dims[0].ge;
  }
  
  px_a = s[D_P0]._dims[0].n0 - nx;
  py_a = s[D_P1]._dims[0].n0 - nx;
  pz_a = s[D_P2]._dims[0].n0 - nx;

  mp00_a = s[D_MP00]._dims[0].n0 - nx;
  mp01_a = s[D_MP01]._dims[0].n0 - nx;
  
  vx_a = s[D_V0]._dims[0].n0 - nx;
  vy_a = s[D_V1]._dims[0].n0 - nx;
  vz_a = s[D_V2]._dims[0].n0 - nx;

  px_pml_I_a = ld_pml[0]._s[D_P0]._dims[0].n0 - nx;
  py_pml_I_a = ld_pml[0]._s[D_P1]._dims[0].n0 - nx;
  pz_pml_I_a = ld_pml[0]._s[D_P2]._dims[0].n0 - nx;
  
  px_pml_II_a = ld_pml[3]._s[D_P0]._dims[0].n0 - nx;
  py_pml_II_a = ld_pml[3]._s[D_P1]._dims[0].n0 - nx;
  pz_pml_II_a = ld_pml[3]._s[D_P2]._dims[0].n0 - nx;
  
  px_pml_III_a = ld_pml[6]._s[D_P0]._dims[0].n0 - nx;
  py_pml_III_a = ld_pml[6]._s[D_P1]._dims[0].n0 - nx;
  pz_pml_III_a = ld_pml[6]._s[D_P2]._dims[0].n0 - nx;
  
  px_pml_IV_a = ld_pml[9]._s[D_P0]._dims[0].n0 - nx;
  py_pml_IV_a = ld_pml[9]._s[D_P1]._dims[0].n0 - nx;
  pz_pml_IV_a = ld_pml[9]._s[D_P2]._dims[0].n0 - nx;
  
  px_pml_V_a = ld_pml[12]._s[D_P0]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  py_pml_V_a = ld_pml[12]._s[D_P1]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  pz_pml_V_a = ld_pml[12]._s[D_P2]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  
  px_pml_VI_a = ld_pml[15]._s[D_P0]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);
  py_pml_VI_a = ld_pml[15]._s[D_P1]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);
  pz_pml_VI_a = ld_pml[15]._s[D_P2]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);  

#pragma omp parallel private(\
  tsz,tid,iy,_px,_py,_pz,_pxend,_mp00,_mp01,_epx,_epy,_epz,_vx3,_vy3,_vy2,_vy1,_vy0,_vz3,_vz2,_vz1,_vz0,vx3,vx2,vx1,vx0, \
  _px_x_I,  _px_y_I,  _px_z_I,  _py_x_I,  _py_y_I,  _py_z_I,  _pz_x_I,  _pz_y_I,  _pz_z_I,   \
  _px_x_II, _px_y_II, _px_z_II, _py_x_II, _py_y_II, _py_z_II, _pz_x_II, _pz_y_II, _pz_z_II,  \
  _px_x_III,_px_y_III,_px_z_III,_py_x_III,_py_y_III,_py_z_III,_pz_x_III,_pz_y_III,_pz_z_III, \
  _px_x_IV, _px_y_IV, _px_z_IV, _py_x_IV, _py_y_IV, _py_z_IV, _pz_x_IV, _pz_y_IV, _pz_z_IV,  \
  _px_x_V,  _px_y_V,  _px_z_V,  _py_x_V,  _py_y_V,  _py_z_V,  _pz_x_V,  _pz_y_V,  _pz_z_V,   \
  _px_x_VI, _px_y_VI, _px_z_VI, _py_x_VI, _py_y_VI, _py_z_VI, _pz_x_VI, _pz_y_VI, _pz_z_VI,   \
  dfdx,dfdy,dfdz,etaxdt,etaydt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      px_aa = (tsz * s[D_P0]._dims[1].n0 - ny) * s[D_P0]._dims[0].n0;
      py_aa = (tsz * s[D_P1]._dims[1].n0 - ny) * s[D_P1]._dims[0].n0;
      pz_aa = (tsz * s[D_P2]._dims[1].n0 - ny) * s[D_P2]._dims[0].n0;
      
      mp00_aa = (tsz * s[D_MP00]._dims[1].n0 - ny) * s[D_MP00]._dims[0].n0;
      mp01_aa = (tsz * s[D_MP01]._dims[1].n0 - ny) * s[D_MP01]._dims[0].n0;
  
      vx_aa = (tsz * s[D_V0]._dims[1].n0 - ny) * s[D_V0]._dims[0].n0;
      vy_aa = (tsz * s[D_V1]._dims[1].n0 - ny) * s[D_V1]._dims[0].n0;
      vz_aa = (tsz * s[D_V2]._dims[1].n0 - ny) * s[D_V2]._dims[0].n0;
      
      px_pml_I_aa   = (tsz * ld_pml[0]._s[D_P0]._dims[1].n0 - ny) * ld_pml[0]._s[D_P0]._dims[0].n0;
      py_pml_I_aa   = (tsz * ld_pml[0]._s[D_P1]._dims[1].n0 - ny) * ld_pml[0]._s[D_P1]._dims[0].n0;
      pz_pml_I_aa   = (tsz * ld_pml[0]._s[D_P2]._dims[1].n0 - ny) * ld_pml[0]._s[D_P2]._dims[0].n0;
      
      px_pml_II_aa  = (tsz * ld_pml[3]._s[D_P0]._dims[1].n0 - ny) * ld_pml[3]._s[D_P0]._dims[0].n0;
      py_pml_II_aa  = (tsz * ld_pml[3]._s[D_P1]._dims[1].n0 - ny) * ld_pml[3]._s[D_P1]._dims[0].n0;
      pz_pml_II_aa  = (tsz * ld_pml[3]._s[D_P2]._dims[1].n0 - ny) * ld_pml[3]._s[D_P2]._dims[0].n0;
      
      px_pml_III_aa = (tsz * ld_pml[6]._s[D_P0]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_P0]._dims[0].n0;
      py_pml_III_aa = (tsz * ld_pml[6]._s[D_P1]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_P1]._dims[0].n0;
      pz_pml_III_aa = (tsz * ld_pml[6]._s[D_P2]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_P2]._dims[0].n0;
      
      px_pml_IV_aa  = (tsz * ld_pml[9]._s[D_P0]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_P0]._dims[0].n0;
      py_pml_IV_aa  = (tsz * ld_pml[9]._s[D_P1]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_P1]._dims[0].n0;
      pz_pml_IV_aa  = (tsz * ld_pml[9]._s[D_P2]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_P2]._dims[0].n0;
      
      px_pml_V_aa   = (tsz * ld_pml[12]._s[D_P0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_P0]._dims[0].n0;
      py_pml_V_aa   = (tsz * ld_pml[12]._s[D_P1]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_P1]._dims[0].n0;
      pz_pml_V_aa   = (tsz * ld_pml[12]._s[D_P2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_P2]._dims[0].n0;
      
      px_pml_VI_aa  = (tsz * ld_pml[15]._s[D_P0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_P0]._dims[0].n0;
      py_pml_VI_aa  = (tsz * ld_pml[15]._s[D_P1]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_P1]._dims[0].n0;
      pz_pml_VI_aa  = (tsz * ld_pml[15]._s[D_P2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_P2]._dims[0].n0;  
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_P0]._dims;
    _px     = s[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_P1]._dims;
    _py     = s[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_P2]._dims;
    _pz     = s[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    dims    = s_pml[D_P0]._dims;
    _px_x_I = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P1]._dims;
    _py_x_I = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P2]._dims;
    _pz_x_I = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[1]._s;
    dims    = s_pml[D_P0]._dims;
    _px_y_I = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P1]._dims;
    _py_y_I = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P2]._dims;
    _pz_y_I = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;

    s_pml   = ld_pml[2]._s;
    dims    = s_pml[D_P0]._dims;
    _px_z_I = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P1]._dims;
    _py_z_I = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims    = s_pml[D_P2]._dims;
    _pz_z_I = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
   
    dims    = s[D_MP00]._dims;
    _mp00   = cs[D_MP00]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_MP01]._dims;
    _mp01   = cs[D_MP01]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_V0]._dims;
    _vx3    = rs[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 1;

    dims    = s[D_V1]._dims;
    _vy2    = rs[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy3    = _vy2 + dims[0].n0;
    _vy1    = _vy2 - dims[0].n0;
    _vy0    = _vy1 - dims[0].n0;

    dims    = s[D_V2]._dims;
    _vz2    = rs[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz3    = _vz2 + dims[1].n0 * dims[0].n0;
    _vz1    = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0    = _vz1 - dims[1].n0 * dims[0].n0;


    _epx    = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EP[2]]._s + (gzs_pml_I + tid - s[D_EP[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _pxend = _px + nx; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_I) = ((*_px_x_I) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_I) = ((*_px_y_I) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_px_z_I) = ((*_px_z_I) * (1.0f - etazdt) + dfdz*(*_mp01))/(1.0f + etazdt);
            *_px       = *_px_x_I + *_px_y_I + * _px_z_I;

            (*_py_x_I) = ((*_py_x_I) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_I) = ((*_py_y_I) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);
            (*_py_z_I) = ((*_py_z_I) * (1.0f - etazdt) + dfdz*(*_mp01))/(1.0f + etazdt);
            *_py       = *_py_x_I + *_py_y_I + * _py_z_I;
            
            (*_pz_x_I) = ((*_pz_x_I) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_I) = ((*_pz_y_I) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_pz_z_I) = ((*_pz_z_I) * (1.0f - etazdt) + dfdz*(*_mp00))/(1.0f + etazdt);
            *_pz       = *_pz_x_I + *_pz_y_I + * _pz_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;
          _px_x_I++; _px_y_I++; _px_z_I++; 
          _py_x_I++; _py_y_I++; _py_z_I++;
          _pz_x_I++; _pz_y_I++; _pz_z_I++;
          
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _px += px_a; _py += py_a; _pz += pz_a; 
        _mp00 += mp00_a; _mp01 += mp01_a; 
        _vx3 += vx_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _epx -= nx;
        _px_x_I += px_pml_I_a; _py_x_I += py_pml_I_a; _pz_x_I += pz_pml_I_a;
        _px_y_I += px_pml_I_a; _py_y_I += py_pml_I_a; _pz_y_I += pz_pml_I_a;
        _px_z_I += px_pml_I_a; _py_z_I += py_pml_I_a; _pz_z_I += pz_pml_I_a; 
      }
      _px += px_aa; _py += py_aa; _pz += pz_aa; 
      _mp00 += mp00_aa; _mp01 += mp01_aa; 
      _vx3 += vx_aa; 
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; 
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
      _epy -= ny; _epz += tsz;
      _px_x_I += px_pml_I_aa; _py_x_I += py_pml_I_aa; _pz_x_I += pz_pml_I_aa;
      _px_y_I += px_pml_I_aa; _py_y_I += py_pml_I_aa; _pz_y_I += pz_pml_I_aa;
      _px_z_I += px_pml_I_aa; _py_z_I += py_pml_I_aa; _pz_z_I += pz_pml_I_aa; 
    }
    /*************************************************************************************/

    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _px   += (gze_pml_I + 1 + tid - iz) * s[D_P0]._dims[0].n0 * s[D_P0]._dims[1].n0;
    _py   += (gze_pml_I + 1 + tid - iz) * s[D_P1]._dims[0].n0 * s[D_P1]._dims[1].n0;
    _pz   += (gze_pml_I + 1 + tid - iz) * s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _mp00 += (gze_pml_I + 1 + tid - iz) * s[D_MP00]._dims[0].n0 * s[D_MP00]._dims[1].n0;
    _mp01 += (gze_pml_I + 1 + tid - iz) * s[D_MP01]._dims[0].n0 * s[D_MP01]._dims[1].n0;
    _vx3  += (gze_pml_I + 1 + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vy2  += (gze_pml_I + 1 + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vz2  += (gze_pml_I + 1 + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;

    _vz3 = _vz2 + s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vz1 = _vz2 - s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vz0 = _vz1 - s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;

    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gze_pml_I + 1 + tid - s[D_EP[2]]._dims[0].gs);        /* 1D */
    
    s_pml     = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_P0]._dims;
    _px_x_III = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_x_III = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_x_III = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[7]._s;
    dims      = s_pml[D_P0]._dims;
    _px_y_III = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_y_III = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_y_III = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[8]._s;
    dims      = s_pml[D_P0]._dims;
    _px_z_III = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_z_III = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_z_III = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_P0]._dims;
    _px_x_IV  = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_x_IV  = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_x_IV  = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[10]._s;
    dims      = s_pml[D_P0]._dims;
    _px_y_IV  = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_y_IV  = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_y_IV  = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_P0]._dims;
    _px_z_IV  = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_z_IV  = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_z_IV  = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_P0]._dims;
    _px_x_V   = s_pml[D_P0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_x_V   = s_pml[D_P1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_x_V   = s_pml[D_P2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[13]._s;
    dims      = s_pml[D_P0]._dims;
    _px_y_V   = s_pml[D_P0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_y_V   = s_pml[D_P1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_y_V   = s_pml[D_P2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[14]._s;
    dims      = s_pml[D_P0]._dims;
    _px_z_V   = s_pml[D_P0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_z_V   = s_pml[D_P1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_z_V   = s_pml[D_P2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_P0]._dims;
    _px_x_VI  = s_pml[D_P0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_x_VI  = s_pml[D_P1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_x_VI  = s_pml[D_P2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[16]._s;
    dims      = s_pml[D_P0]._dims;
    _px_y_VI  = s_pml[D_P0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_y_VI  = s_pml[D_P1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_y_VI  = s_pml[D_P2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[17]._s;
    dims      = s_pml[D_P0]._dims;
    _px_z_VI  = s_pml[D_P0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P1]._dims;
    _py_z_VI  = s_pml[D_P1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    dims      = s_pml[D_P2]._dims;
    _pz_z_VI  = s_pml[D_P2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _pxend = _px + nx; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_III) = ((*_px_x_III) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_III) = ((*_px_y_III) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_px_z_III) =  (*_px_z_III) + dfdz*(*_mp01);
            *_px         = *_px_x_III + *_px_y_III + *_px_z_III;
                
            (*_py_x_III) = ((*_py_x_III) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_III) = ((*_py_y_III) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);
            (*_py_z_III) =  (*_py_z_III) + dfdz*(*_mp01);
            *_py         = *_py_x_III + *_py_y_III + *_py_z_III;
            
            (*_pz_x_III) = ((*_pz_x_III) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_III) = ((*_pz_y_III) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_pz_z_III) =  (*_pz_z_III) + dfdz*(*_mp00);
            *_pz         = *_pz_x_III + *_pz_y_III + *_pz_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;
          _px_x_III++; _px_y_III++; _px_z_III++;
          _py_x_III++; _py_y_III++; _py_z_III++;
          _pz_x_III++; _pz_y_III++; _pz_z_III++;
	
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        _px += px_a; _py += py_a; _pz += pz_a; 
        _mp00 += mp00_a; _mp01 += mp01_a;
        _vx3 += vx_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a;
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _epx -= nx;
        _px_x_III += px_pml_III_a; _py_x_III += py_pml_III_a; _pz_x_III += pz_pml_III_a;
        _px_y_III += px_pml_III_a; _py_y_III += py_pml_III_a; _pz_y_III += pz_pml_III_a;
        _px_z_III += px_pml_III_a; _py_z_III += py_pml_III_a; _pz_z_III += pz_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
        _epy ++;
        for ( _pxend = _px + gxe_pml_V-gxs_pml_V+1; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_V) = ((*_px_x_V) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_V) =  (*_px_y_V) + dfdy*(*_mp01);
            (*_px_z_V) =  (*_px_z_V) + dfdz*(*_mp01);
            *_px       = *_px_x_V + *_px_y_V + *_px_z_V;
            
            (*_py_x_V) = ((*_py_x_V) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_V) =  (*_py_y_V) + dfdy*(*_mp00);
            (*_py_z_V) =  (*_py_z_V) + dfdz*(*_mp01);
            *_py       = *_py_x_V + *_py_y_V + *_py_z_V;
            
            (*_pz_x_V) = ((*_pz_x_V) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_V) =  (*_pz_y_V) + dfdy*(*_mp01);
            (*_pz_z_V) =  (*_pz_z_V) + dfdz*(*_mp00);
            *_pz       = *_pz_x_V + *_pz_y_V + *_pz_z_V; 
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;	
          _px_x_V++; _px_y_V++; _px_z_V++;
          _py_x_V++; _py_y_V++; _py_z_V++;
          _pz_x_V++; _pz_y_V++; _pz_z_V++;
	
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        
        /** physical region */
        for ( _pxend = _px + gxs_pml_VI-gxe_pml_V-1; _px < _pxend;) {
          vx3 = *_vx3++;
          _epx ++;
          
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            *_px = *_px + dfdx*(*_mp00) + dfdy*(*_mp01) + dfdz*(*_mp01);
            *_py = *_py + dfdx*(*_mp01) + dfdy*(*_mp00) + dfdz*(*_mp01);
            *_pz = *_pz + dfdx*(*_mp01) + dfdy*(*_mp01) + dfdz*(*_mp00);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;	
          vx0 = vx1; vx1 = vx2; vx2 = vx3;          
        }
        
        /** pml region VI */
        for ( _pxend = _px + gxe_pml_VI-gxs_pml_VI+1; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_VI) = ((*_px_x_VI) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_VI) =  (*_px_y_VI) + dfdy*(*_mp01);
            (*_px_z_VI) =  (*_px_z_VI) + dfdz*(*_mp01);
            *_px        = *_px_x_VI + *_px_y_VI + *_px_z_VI;
            
            (*_py_x_VI) = ((*_py_x_VI) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_VI) =  (*_py_y_VI) + dfdy*(*_mp00);
            (*_py_z_VI) =  (*_py_z_VI) + dfdz*(*_mp01);
            *_py        = *_py_x_VI + *_py_y_VI + *_py_z_VI;
            
            (*_pz_x_VI) = ((*_pz_x_VI) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_VI) =  (*_pz_y_VI) + dfdy*(*_mp01);
            (*_pz_z_VI) =  (*_pz_z_VI) + dfdz*(*_mp00);
            *_pz        = *_pz_x_VI + *_pz_y_VI + *_pz_z_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;
          _px_x_VI++; _px_y_VI++; _px_z_VI++;
          _py_x_VI++; _py_y_VI++; _py_z_VI++;
          _pz_x_VI++; _pz_y_VI++; _pz_z_VI++;
	        
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _px += px_a; _py += py_a; _pz += pz_a; 
        _mp00 += mp00_a; _mp01 += mp01_a; 
        _vx3 += vx_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _epx -= nx;
        _px_x_V += px_pml_V_a; _py_x_V += py_pml_V_a; _pz_x_V += pz_pml_V_a;
        _px_y_V += px_pml_V_a; _py_y_V += py_pml_V_a; _pz_y_V += pz_pml_V_a;
        _px_z_V += px_pml_V_a; _py_z_V += py_pml_V_a; _pz_z_V += pz_pml_V_a;
        _px_x_VI += px_pml_VI_a; _py_x_VI += py_pml_VI_a; _pz_x_VI += pz_pml_VI_a;
        _px_y_VI += px_pml_VI_a; _py_y_VI += py_pml_VI_a; _pz_y_VI += pz_pml_VI_a;
        _px_z_VI += px_pml_VI_a; _py_z_VI += py_pml_VI_a; _pz_z_VI += pz_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _pxend = _px + nx; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_IV) = ((*_px_x_IV) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_IV) = ((*_px_y_IV) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_px_z_IV) =  (*_px_z_IV) + dfdz*(*_mp01);
            *_px        = *_px_x_IV + *_px_y_IV + *_px_z_IV;
                
            (*_py_x_IV) = ((*_py_x_IV) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_IV) = ((*_py_y_IV) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);
            (*_py_z_IV) =  (*_py_z_IV) + dfdz*(*_mp01);
            *_py        = *_py_x_IV + *_py_y_IV + *_py_z_IV;
            
            (*_pz_x_IV) = ((*_pz_x_IV) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_IV) = ((*_pz_y_IV) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_pz_z_IV) =  (*_pz_z_IV) + dfdz*(*_mp00);
            *_pz        = *_pz_x_IV + *_pz_y_IV + *_pz_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;
          _px_x_IV++; _px_y_IV++; _px_z_IV++;
          _py_x_IV++; _py_y_IV++; _py_z_IV++;
          _pz_x_IV++; _pz_y_IV++; _pz_z_IV++;
	
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        _px += px_a; _py += py_a; _pz += pz_a; 
        _mp00 += mp00_a; _mp01 += mp01_a;
        _vx3 += vx_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a;
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _epx -= nx;
        _px_x_IV += px_pml_IV_a; _py_x_IV += py_pml_IV_a; _pz_x_IV += pz_pml_IV_a;
        _px_y_IV += px_pml_IV_a; _py_y_IV += py_pml_IV_a; _pz_y_IV += pz_pml_IV_a;
        _px_z_IV += px_pml_IV_a; _py_z_IV += py_pml_IV_a; _pz_z_IV += pz_pml_IV_a;
      }
      _epy -= ny;
      _px += px_aa; _py += py_aa; _pz += pz_aa; 
      _mp00 += mp00_aa; _mp01 += mp01_aa;
      _vx3 += vx_aa;
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa;
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
     
      _px_x_III += px_pml_III_aa; _py_x_III += py_pml_III_aa; _pz_x_III += pz_pml_III_aa;
      _px_y_III += px_pml_III_aa; _py_y_III += py_pml_III_aa; _pz_y_III += pz_pml_III_aa;
      _px_z_III += px_pml_III_aa; _py_z_III += py_pml_III_aa; _pz_z_III += pz_pml_III_aa;

      _px_x_IV += px_pml_IV_aa; _py_x_IV += py_pml_IV_aa; _pz_x_IV += pz_pml_IV_aa;
      _px_y_IV += px_pml_IV_aa; _py_y_IV += py_pml_IV_aa; _pz_y_IV += pz_pml_IV_aa;
      _px_z_IV += px_pml_IV_aa; _py_z_IV += py_pml_IV_aa; _pz_z_IV += pz_pml_IV_aa;

      _px_x_V += px_pml_V_aa; _py_x_V += py_pml_V_aa; _pz_x_V += pz_pml_V_aa;
      _px_y_V += px_pml_V_aa; _py_y_V += py_pml_V_aa; _pz_y_V += pz_pml_V_aa;
      _px_z_V += px_pml_V_aa; _py_z_V += py_pml_V_aa; _pz_z_V += pz_pml_V_aa;

      _px_x_VI += px_pml_VI_aa; _py_x_VI += py_pml_VI_aa; _pz_x_VI += pz_pml_VI_aa;
      _px_y_VI += px_pml_VI_aa; _py_y_VI += py_pml_VI_aa; _pz_y_VI += pz_pml_VI_aa;
      _px_z_VI += px_pml_VI_aa; _py_z_VI += py_pml_VI_aa; _pz_z_VI += pz_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _px   += (gzs_pml_II + tid - iz) * s[D_P0]._dims[0].n0 * s[D_P0]._dims[1].n0;
    _py   += (gzs_pml_II + tid - iz) * s[D_P1]._dims[0].n0 * s[D_P1]._dims[1].n0;
    _pz   += (gzs_pml_II + tid - iz) * s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _mp00 += (gzs_pml_II + tid - iz) * s[D_MP00]._dims[0].n0 * s[D_MP00]._dims[1].n0;
    _mp01 += (gzs_pml_II + tid - iz) * s[D_MP01]._dims[0].n0 * s[D_MP01]._dims[1].n0;
    _vx3  += (gzs_pml_II + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vy2  += (gzs_pml_II + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vz2  += (gzs_pml_II + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;

    _vz3 = _vz2 + s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vz1 = _vz2 - s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vz0 = _vz1 - s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
  
    s_pml    = ld_pml[3]._s;
    dims     = s_pml[D_P0]._dims;
    _px_x_II = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P1]._dims;
    _py_x_II = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P2]._dims;
    _pz_x_II = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[4]._s;
    dims     = s_pml[D_P0]._dims;
    _px_y_II = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P1]._dims;
    _py_y_II = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P2]._dims;
    _pz_y_II = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[5]._s;
    dims     = s_pml[D_P0]._dims;
    _px_z_II = s_pml[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P1]._dims;
    _py_z_II = s_pml[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims     = s_pml[D_P2]._dims;
    _pz_z_II = s_pml[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gzs_pml_II + tid - s[D_EP[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _pxend = _px + nx; _px < _pxend; ) {
          vx3 = *_vx3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
          dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
          dfdz = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_px_x_II) = ((*_px_x_II) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
            (*_px_y_II) = ((*_px_y_II) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_px_z_II) = ((*_px_z_II) * (1.0f - etazdt) + dfdz*(*_mp01))/(1.0f + etazdt);
            *_px        = *_px_x_II + *_px_y_II + * _px_z_II;

            (*_py_x_II) = ((*_py_x_II) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_py_y_II) = ((*_py_y_II) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);
            (*_py_z_II) = ((*_py_z_II) * (1.0f - etazdt) + dfdz*(*_mp01))/(1.0f + etazdt);
            *_py        = *_py_x_II + *_py_y_II + * _py_z_II;
            
            (*_pz_x_II) = ((*_pz_x_II) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
            (*_pz_y_II) = ((*_pz_y_II) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
            (*_pz_z_II) = ((*_pz_z_II) * (1.0f - etazdt) + dfdz*(*_mp00))/(1.0f + etazdt);
            *_pz        = *_pz_x_II + *_pz_y_II + * _pz_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _px++; _py++; _pz++; _mp00++; _mp01++;
          _px_x_II++; _px_y_II++; _px_z_II++; 
          _py_x_II++; _py_y_II++; _py_z_II++;
          _pz_x_II++; _pz_y_II++; _pz_z_II++;
          
          vx0 = vx1; vx1 = vx2; vx2 = vx3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _px += px_a; _py += py_a; _pz += pz_a; 
        _mp00 += mp00_a; _mp01 += mp01_a; 
        _vx3 += vx_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _epx -= nx;
        _px_x_II += px_pml_II_a; _py_x_II += py_pml_II_a; _pz_x_II += pz_pml_II_a;
        _px_y_II += px_pml_II_a; _py_y_II += py_pml_II_a; _pz_y_II += pz_pml_II_a;
        _px_z_II += px_pml_II_a; _py_z_II += py_pml_II_a; _pz_z_II += pz_pml_II_a; 
      }
      _px += px_aa; _py += py_aa; _pz += pz_aa; 
      _mp00 += mp00_aa; _mp01 += mp01_aa; 
      _vx3 += vx_aa; 
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa;
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
      _epy -= ny; _epz += tsz;
      _px_x_II += px_pml_II_aa; _py_x_II += py_pml_II_aa; _pz_x_II += pz_pml_II_aa;
      _px_y_II += px_pml_II_aa; _py_y_II += py_pml_II_aa; _pz_y_II += pz_pml_II_aa;
      _px_z_II += px_pml_II_aa; _py_z_II += py_pml_II_aa; _pz_z_II += pz_pml_II_aa; 
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** shear stress s0 (sxy) */
int esgn_gts3d_24ss0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, sxy_a, ms0_a, vx_a, vy_a;
  int sxy_aa, ms0_aa, vx_aa, vy_aa, iy, iz, tsz, tid;
  int sxy_pml_I_a,   sxy_pml_I_aa;
  int sxy_pml_II_a,  sxy_pml_II_aa;
  int sxy_pml_III_a, sxy_pml_III_aa;
  int sxy_pml_IV_a,  sxy_pml_IV_aa;
  int sxy_pml_V_a,   sxy_pml_V_aa;
  int sxy_pml_VI_a,  sxy_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _sxy, * restrict _sxyend;
  register ireal * restrict _sxy_x_I,   * restrict _sxy_y_I; 
  register ireal * restrict _sxy_x_II,  * restrict _sxy_y_II;
  register ireal * restrict _sxy_x_III, * restrict _sxy_y_III;
  register ireal * restrict _sxy_x_IV,  * restrict _sxy_y_IV;
  register ireal * restrict _sxy_x_V,   * restrict _sxy_y_V;
  register ireal * restrict _sxy_x_VI,  * restrict _sxy_y_VI;
  register ireal * restrict _ms0;
  register ireal * restrict _vy3;
  register ireal * restrict _vx3, * restrict _vx2, * restrict _vx1, * restrict _vx0;
  register ireal * restrict _epx, * restrict _epy;
  register ireal lax, lay, dt2, vy3, vy2, vy1, vy0, dfdx, dfdy, etaxdt, etaydt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_S0]._dims[0].n;
  ny = s[D_S0]._dims[1].n;
  nz = s[D_S0]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_S0]._dims[0].gs;
  gxe = s[D_S0]._dims[0].ge;
  gys = s[D_S0]._dims[1].gs;
  gye = s[D_S0]._dims[1].ge;
  gzs = s[D_S0]._dims[2].gs;
  gze = s[D_S0]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
 
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_S0,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_S0]._dims[2].gs;
    gze_pml_I = s_pml[D_S0]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_S0,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_S0]._dims[2].gs;
    gze_pml_II = s_pml[D_S0]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_S0,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_S0]._dims[1].gs;
    gye_pml_III = s_pml[D_S0]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_S0,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_S0]._dims[1].gs;
    gye_pml_IV = s_pml[D_S0]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_S0,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_S0]._dims[0].gs;
    gxe_pml_V = s_pml[D_S0]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_S0,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_S0]._dims[0].gs;
    gxe_pml_VI = s_pml[D_S0]._dims[0].ge;
  }
  
  sxy_a = s[D_S0]._dims[0].n0 - nx;
 
  ms0_a = s[D_MS0]._dims[0].n0 - nx;
    
  vx_a = s[D_V0]._dims[0].n0 - nx;
  vy_a = s[D_V1]._dims[0].n0 - nx;
 
  sxy_pml_I_a   = ld_pml[0 ]._s[D_S0]._dims[0].n0 - nx;   
  sxy_pml_II_a  = ld_pml[3 ]._s[D_S0]._dims[0].n0 - nx;
  sxy_pml_III_a = ld_pml[6 ]._s[D_S0]._dims[0].n0 - nx;
  sxy_pml_IV_a  = ld_pml[9 ]._s[D_S0]._dims[0].n0 - nx;
  sxy_pml_V_a   = ld_pml[12]._s[D_S0]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  sxy_pml_VI_a  = ld_pml[15]._s[D_S0]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);

#pragma omp parallel private(\
  tsz,tid,iy,_sxy,_sxyend,_ms0,_epx,_epy,_vy3,_vx3,_vx2,_vx1,_vx0,vy3,vy2,vy1,vy0, \
  _sxy_x_I,  _sxy_y_I,  _sxy_x_II, _sxy_y_II, _sxy_x_III, _sxy_y_III, \
  _sxy_x_IV, _sxy_y_IV, _sxy_x_V,  _sxy_y_V,  _sxy_x_VI,  _sxy_y_VI,  \
  dfdx,dfdy,etaxdt,etaydt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      sxy_aa = (tsz * s[D_S0]._dims[1].n0 - ny) * s[D_S0]._dims[0].n0;
      
      ms0_aa = (tsz * s[D_MS0]._dims[1].n0 - ny) * s[D_MS0]._dims[0].n0;
       
      vx_aa = (tsz * s[D_V0]._dims[1].n0 - ny) * s[D_V0]._dims[0].n0;
      vy_aa = (tsz * s[D_V1]._dims[1].n0 - ny) * s[D_V1]._dims[0].n0;
            
      sxy_pml_I_aa   = (tsz * ld_pml[0]._s[D_S0]._dims[1].n0 - ny) * ld_pml[0]._s[D_S0]._dims[0].n0;
      sxy_pml_II_aa  = (tsz * ld_pml[3]._s[D_S0]._dims[1].n0 - ny) * ld_pml[3]._s[D_S0]._dims[0].n0;
      sxy_pml_III_aa = (tsz * ld_pml[6]._s[D_S0]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_S0]._dims[0].n0;
      sxy_pml_IV_aa  = (tsz * ld_pml[9]._s[D_S0]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_S0]._dims[0].n0;
      sxy_pml_V_aa   = (tsz * ld_pml[12]._s[D_S0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_S0]._dims[0].n0;
      sxy_pml_VI_aa  = (tsz * ld_pml[15]._s[D_S0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_S0]._dims[0].n0;
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims     = s[D_S0]._dims;
    _sxy     = s[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    s_pml    = ld_pml[0]._s;
    dims     = s_pml[D_S0]._dims;
    _sxy_x_I = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml    = ld_pml[1]._s;
    dims     = s_pml[D_S0]._dims;
    _sxy_y_I = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
            
    dims     = s[D_MS0]._dims;
    _ms0     = cs[D_MS0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_V1]._dims;
    _vy3    = rs[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 2;

    dims    = s[D_V0]._dims;
    _vx1    = rs[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vx2    = _vx1 + dims[0].n0;
    _vx3    = _vx2 + dims[0].n0;
    _vx0    = _vx1 - dims[0].n0;

    _epx    = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);              /* 1D */
       
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      for ( iy = 0; iy < ny; iy ++ ) {
        vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxyend = _sxy + nx; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_I) = ((*_sxy_x_I) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_I) = ((*_sxy_y_I) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
            *_sxy       = *_sxy_x_I + *_sxy_y_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          _sxy_x_I++; _sxy_y_I++;
                  
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxy += sxy_a;
        _ms0 += ms0_a;
        _vy3 += vy_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; 
        _epx -= nx;
        _sxy_x_I += sxy_pml_I_a; 
        _sxy_y_I += sxy_pml_I_a;
      }
      _sxy += sxy_aa;
      _ms0 += ms0_aa; 
      _vy3 += vy_aa; 
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa; 
      _epy -= ny;
      _sxy_x_I += sxy_pml_I_aa; 
      _sxy_y_I += sxy_pml_I_aa; 
    }
    /*************************************************************************************/

    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _sxy  += (gze_pml_I + 1 + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _ms0  += (gze_pml_I + 1 + tid - iz) * s[D_MS0]._dims[0].n0 * s[D_MS0]._dims[1].n0;
    _vy3  += (gze_pml_I + 1 + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vx2  += (gze_pml_I + 1 + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
       
    _vx3 = _vx2 + s[D_V0]._dims[0].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0;

    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
        
    s_pml     = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims       = s_pml[D_S0]._dims;
    _sxy_x_III = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml      = ld_pml[7]._s;
    dims       = s_pml[D_S0]._dims;
    _sxy_y_III = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S0]._dims;
    _sxy_x_IV = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[10]._s;
    dims      = s_pml[D_S0]._dims;
    _sxy_y_IV = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S0]._dims;
    _sxy_x_V  = s_pml[D_S0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[13]._s;
    dims      = s_pml[D_S0]._dims;
    _sxy_y_V  = s_pml[D_S0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S0]._dims;
    _sxy_x_VI = s_pml[D_S0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[16]._s;
    dims      = s_pml[D_S0]._dims;
    _sxy_y_VI = s_pml[D_S0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
            
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxyend = _sxy + nx; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_III) = ((*_sxy_x_III) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_III) = ((*_sxy_y_III) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
            *_sxy         = *_sxy_x_III + *_sxy_y_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          _sxy_x_III++; _sxy_y_III++;
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        _sxy += sxy_a;
        _ms0 += ms0_a;
        _vy3 += vy_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a;
        _epx -= nx;
        _sxy_x_III += sxy_pml_III_a;
        _sxy_y_III += sxy_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
        _epy ++;
        for ( _sxyend = _sxy + gxe_pml_V-gxs_pml_V+1; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_V) = ((*_sxy_x_V) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_V) =  (*_sxy_y_V) + dfdy*(*_ms0);
            *_sxy       = *_sxy_x_V + *_sxy_y_V;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;	
          _sxy_x_V++; _sxy_y_V++;
          
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        
        /** physical region */
        for ( _sxyend = _sxy + gxs_pml_VI-gxe_pml_V-1; _sxy < _sxyend;) {
          vy3 = *_vy3++;
          _epx ++;
          
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            *_sxy = *_sxy + (dfdx + dfdy) * (*_ms0);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          vy0 = vy1; vy1 = vy2; vy2 = vy3;          
        }
        
        /** pml region VI */
        for ( _sxyend = _sxy + gxe_pml_VI-gxs_pml_VI+1; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_VI) = ((*_sxy_x_VI) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_VI) =  (*_sxy_y_VI) + dfdy*(*_ms0);
            *_sxy        = *_sxy_x_VI + *_sxy_y_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          _sxy_x_VI++; _sxy_y_VI++;
         	        
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxy += sxy_a;
        _ms0 += ms0_a;
        _vy3 += vy_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; 
        _epx -= nx;
        _sxy_x_V  += sxy_pml_V_a; 
        _sxy_y_V  += sxy_pml_V_a; 
        _sxy_x_VI += sxy_pml_VI_a;
        _sxy_y_VI += sxy_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxyend = _sxy + nx; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_IV) = ((*_sxy_x_IV) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_IV) = ((*_sxy_y_IV) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
            *_sxy        = *_sxy_x_IV + *_sxy_y_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          _sxy_x_IV++; _sxy_y_IV++;
          	
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        _sxy += sxy_a;
        _ms0 += ms0_a;
        _vy3 += vy_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a;
        _epx -= nx;
        _sxy_x_IV += sxy_pml_IV_a;
        _sxy_y_IV += sxy_pml_IV_a;
      }
      _epy -= ny;
      _sxy += sxy_aa;
      _ms0 += ms0_aa;
      _vy3 += vy_aa;
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa;
           
      _sxy_x_III += sxy_pml_III_aa; 
      _sxy_y_III += sxy_pml_III_aa; 
      
      _sxy_x_IV += sxy_pml_IV_aa; 
      _sxy_y_IV += sxy_pml_IV_aa; 
      
      _sxy_x_V += sxy_pml_V_aa;
      _sxy_y_V += sxy_pml_V_aa;
      
      _sxy_x_VI += sxy_pml_VI_aa;
      _sxy_y_VI += sxy_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _sxy  += (gzs_pml_II + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _ms0  += (gzs_pml_II + tid - iz) * s[D_MS0]._dims[0].n0 * s[D_MS0]._dims[1].n0;
    _vy3  += (gzs_pml_II + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vx2  += (gzs_pml_II + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
        
    _vx3 = _vx2 + s[D_V0]._dims[0].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0;

    s_pml     = ld_pml[3]._s;
    dims      = s_pml[D_S0]._dims;
    _sxy_x_II = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[4]._s;
    dims      = s_pml[D_S0]._dims;
    _sxy_y_II = s_pml[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
     
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
   
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) {
      for ( iy = 0; iy < ny; iy ++) {
        vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxyend = _sxy + nx; _sxy < _sxyend; ) {
          vy3 = *_vy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
          dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
          if(_fwd) {
            (*_sxy_x_II) = ((*_sxy_x_II) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
            (*_sxy_y_II) = ((*_sxy_y_II) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
            *_sxy        = *_sxy_x_II + *_sxy_y_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxy++; _ms0++;
          _sxy_x_II++; _sxy_y_II++;
          
          vy0 = vy1; vy1 = vy2; vy2 = vy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxy += sxy_a;
        _ms0 += ms0_a;
        _vy3 += vy_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; 
        _epx -= nx;
        _sxy_x_II += sxy_pml_II_a;
        _sxy_y_II += sxy_pml_II_a;
      }
      _sxy += sxy_aa;
      _ms0 += ms0_aa;
      _vy3 += vy_aa; 
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa;
      _epy -= ny;
      _sxy_x_II += sxy_pml_II_aa;
      _sxy_y_II += sxy_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** shear stress s1 (syz) */
int esgn_gts3d_24ss1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, syz_a, ms1_a, vy_a, vz_a;
  int syz_aa, ms1_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
  int syz_pml_I_a,   syz_pml_I_aa;
  int syz_pml_II_a,  syz_pml_II_aa;
  int syz_pml_III_a, syz_pml_III_aa;
  int syz_pml_IV_a,  syz_pml_IV_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  register ireal * restrict _syz, * restrict _syzend;
  register ireal * restrict _syz_y_I,   * restrict _syz_z_I;
  register ireal * restrict _syz_y_II,  * restrict _syz_z_II;
  register ireal * restrict _syz_y_III, * restrict _syz_z_III;
  register ireal * restrict _syz_y_IV,  * restrict _syz_z_IV;
  register ireal * restrict _ms1;
  register ireal * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
  register ireal * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  register ireal * restrict _epy, * restrict _epz;
  register ireal lay, laz, dt2, dfdy, dfdz, etaydt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_S1]._dims[0].n;
  ny = s[D_S1]._dims[1].n;
  nz = s[D_S1]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_S1]._dims[0].gs;
  gxe = s[D_S1]._dims[0].ge;
  gys = s[D_S1]._dims[1].gs;
  gye = s[D_S1]._dims[1].ge;
  gzs = s[D_S1]._dims[2].gs;
  gze = s[D_S1]._dims[2].ge;
  
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_S1,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_S1]._dims[2].gs;
    gze_pml_I = s_pml[D_S1]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_S1,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_S1]._dims[2].gs;
    gze_pml_II = s_pml[D_S1]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_S1,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_S1]._dims[1].gs;
    gye_pml_III = s_pml[D_S1]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_S1,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_S1]._dims[1].gs;
    gye_pml_IV = s_pml[D_S1]._dims[1].ge;
  }
    
  syz_a = s[D_S1]._dims[0].n0 - nx;
  
  ms1_a = s[D_MS1]._dims[0].n0 - nx;
    
  vy_a = s[D_V1]._dims[0].n0 - nx;
  vz_a = s[D_V2]._dims[0].n0 - nx;

  syz_pml_I_a   = ld_pml[0]._s[D_S1]._dims[0].n0 - nx;
  syz_pml_II_a  = ld_pml[3]._s[D_S1]._dims[0].n0 - nx;
  syz_pml_III_a = ld_pml[6]._s[D_S1]._dims[0].n0 - nx;
  syz_pml_IV_a  = ld_pml[9]._s[D_S1]._dims[0].n0 - nx;
  
#pragma omp parallel private(\
  tsz,tid,iy,_syz,_syzend,_ms1,_epy,_epz,_vy3,_vy2,_vy1,_vy0,_vz3,_vz2,_vz1,_vz0, \
  _syz_y_I,   _syz_z_I,  \
  _syz_y_II,  _syz_z_II, \
  _syz_y_III, _syz_z_III,\
  _syz_y_IV,  _syz_z_IV, \
  dfdy,dfdz,etaydt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      syz_aa = (tsz * s[D_S1]._dims[1].n0 - ny) * s[D_S1]._dims[0].n0;
          
      ms1_aa = (tsz * s[D_MS1]._dims[1].n0 - ny) * s[D_MS1]._dims[0].n0;
      
      vy_aa = (tsz * s[D_V1]._dims[1].n0 - ny) * s[D_V1]._dims[0].n0;
      vz_aa = (tsz * s[D_V2]._dims[1].n0 - ny) * s[D_V2]._dims[0].n0;
      
      syz_pml_I_aa   = (tsz * ld_pml[0]._s[D_S1]._dims[1].n0 - ny) * ld_pml[0]._s[D_S1]._dims[0].n0;
      syz_pml_II_aa  = (tsz * ld_pml[3]._s[D_S1]._dims[1].n0 - ny) * ld_pml[3]._s[D_S1]._dims[0].n0;
      syz_pml_III_aa = (tsz * ld_pml[6]._s[D_S1]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_S1]._dims[0].n0;
      syz_pml_IV_aa  = (tsz * ld_pml[9]._s[D_S1]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_S1]._dims[0].n0;
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_S1]._dims;
    _syz    = s[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
             
    s_pml    = ld_pml[1]._s;
    dims     = s_pml[D_S1]._dims;
    _syz_y_I = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
   
    s_pml    = ld_pml[2]._s;
    dims     = s_pml[D_S1]._dims;
    _syz_z_I = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_MS1]._dims;
    _ms1    = cs[D_MS1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
   
    dims    = s[D_V1]._dims;
    _vy1    = rs[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy2    = _vy1 + dims[0].n0 * dims[1].n0;
    _vy3    = _vy2 + dims[0].n0 * dims[1].n0;
    _vy0    = _vy1 - dims[0].n0 * dims[1].n0;

    dims    = s[D_V2]._dims;
    _vz1    = rs[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz2    = _vz1 + dims[0].n0;
    _vz3    = _vz2 + dims[0].n0;
    _vz0    = _vz1 - dims[0].n0;

    _epy    = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EV[2]]._s + (gzs_pml_I + tid - s[D_EV[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _syzend = _syz + nx; _syz < _syzend; ) {
          dfdy = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * lay;
          dfdz = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * laz;
          if(_fwd) {
            (*_syz_y_I) = ((*_syz_y_I) * (1.0f - etaydt) + dfdy*(*_ms1))/(1.0f + etaydt);
            (*_syz_z_I) = ((*_syz_z_I) * (1.0f - etazdt) + dfdz*(*_ms1))/(1.0f + etazdt);
            *_syz       = *_syz_y_I + * _syz_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _syz++; _ms1++;
          _syz_y_I++; _syz_z_I++; 
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _syz += syz_a; 
        _ms1 += ms1_a;
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _syz_y_I += syz_pml_I_a;
        _syz_z_I += syz_pml_I_a;
      }
      _syz += syz_aa; 
      _ms1 += ms1_aa;
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; 
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
      _epy -= ny; _epz += tsz;
      _syz_y_I += syz_pml_I_aa;
      _syz_z_I += syz_pml_I_aa;
    }
    /*************************************************************************************/

    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _syz  += (gze_pml_I + 1 + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _ms1  += (gze_pml_I + 1 + tid - iz) * s[D_MS1]._dims[0].n0 * s[D_MS1]._dims[1].n0;
    _vy2  += (gze_pml_I + 1 + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vz2  += (gze_pml_I + 1 + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;

    _vz3 = _vz2 + s[D_V2]._dims[0].n0;
    _vz1 = _vz2 - s[D_V2]._dims[0].n0;
    _vz0 = _vz1 - s[D_V2]._dims[0].n0;

    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gze_pml_I + 1 + tid - s[D_EV[2]]._dims[0].gs);        /* 1D */
    
    s_pml      = ld_pml[7]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims       = s_pml[D_S1]._dims;
    _syz_y_III = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml      = ld_pml[8]._s;
    dims       = s_pml[D_S1]._dims;
    _syz_z_III = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
          
    s_pml     = ld_pml[10]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S1]._dims;
    _syz_y_IV = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_S1]._dims;
    _syz_z_IV = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _syzend = _syz + nx; _syz < _syzend; ) {
          dfdy = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * lay;
          dfdz = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * laz;
          if(_fwd) {
            (*_syz_y_III) = ((*_syz_y_III) * (1.0f - etaydt) + dfdy*(*_ms1))/(1.0f + etaydt);
            (*_syz_z_III) =  (*_syz_z_III) + dfdz*(*_ms1);
            *_syz         = *_syz_y_III + *_syz_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _syz++; _ms1++;
          _syz_y_III++; _syz_z_III++;
        }
        _syz += syz_a; 
        _ms1 += ms1_a;
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a;
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _syz_y_III += syz_pml_III_a;
        _syz_z_III += syz_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      /** for shear stress s1, the formulars in these 3 regions are the same */ 
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        _epy ++;
        for ( _syzend = _syz + nx; _syz < _syzend;) {
          dfdy = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * lay;
          dfdz = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * laz;
          if(_fwd) {
            *_syz = *_syz + (dfdy + dfdz) * (*_ms1);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _syz++; _ms1++;	
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _syz += syz_a; 
        _ms1 += ms1_a; 
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _syzend = _syz + nx; _syz < _syzend; ) {
          dfdy = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * lay;
          dfdz = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * laz;
          if(_fwd) {
            (*_syz_y_IV) = ((*_syz_y_IV) * (1.0f - etaydt) + dfdy*(*_ms1))/(1.0f + etaydt);
            (*_syz_z_IV) =  (*_syz_z_IV) + dfdz*(*_ms1);
            *_syz        = *_syz_y_IV + *_syz_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _syz++; _ms1++;
          _syz_y_IV++; _syz_z_IV++;
        }
        _syz += syz_a; 
        _ms1 += ms1_a;
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a;
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _syz_y_IV += syz_pml_IV_a;
        _syz_z_IV += syz_pml_IV_a;
      }
      _epy -= ny;
      _syz += syz_aa; 
      _ms1 += ms1_aa;
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa;
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
     
      _syz_y_III += syz_pml_III_aa;
      _syz_z_III += syz_pml_III_aa;

      _syz_y_IV += syz_pml_IV_aa;
      _syz_z_IV += syz_pml_IV_aa;
    }
    /* redtom */
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _syz  += (gzs_pml_II + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _ms1  += (gzs_pml_II + tid - iz) * s[D_MS1]._dims[0].n0 * s[D_MS1]._dims[1].n0;
    _vy2  += (gzs_pml_II + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vz2  += (gzs_pml_II + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;

    _vz3 = _vz2 + s[D_V2]._dims[0].n0;
    _vz1 = _vz2 - s[D_V2]._dims[0].n0;
    _vz0 = _vz1 - s[D_V2]._dims[0].n0;
  
    s_pml     = ld_pml[4]._s;
    dims      = s_pml[D_S1]._dims;
    _syz_y_II = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[5]._s;
    dims      = s_pml[D_S1]._dims;
    _syz_z_II = s_pml[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gzs_pml_II + tid - s[D_EV[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _syzend = _syz + nx; _syz < _syzend; ) {
          dfdy = ((*_vz3++) - (*_vz0++) + ((*_vz1++) - (*_vz2++)) * 27.0) * lay;
          dfdz = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * laz;
          if(_fwd) {
            (*_syz_y_II) = ((*_syz_y_II) * (1.0f - etaydt) + dfdy*(*_ms1))/(1.0f + etaydt);
            (*_syz_z_II) = ((*_syz_z_II) * (1.0f - etazdt) + dfdz*(*_ms1))/(1.0f + etazdt);
            *_syz        = *_syz_y_II + * _syz_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _syz++; _ms1++;
          _syz_y_II++; _syz_z_II++; 
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _syz += syz_a; 
        _ms1 += ms1_a;
        _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; 
        _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a;
        _syz_y_II += syz_pml_II_a;
        _syz_z_II += syz_pml_II_a;
      }
      _syz += syz_aa; 
      _ms1 += ms1_aa;
      _vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa;
      _vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa;
      _epy -= ny; _epz += tsz;
      _syz_y_II += syz_pml_II_aa;
      _syz_z_II += syz_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** shear stress s2 (sxz) */
int esgn_gts3d_24ss2(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, sxz_a, ms2_a, vx_a, vz_a;
  int sxz_aa, ms2_aa, vx_aa, vz_aa, iy, iz, tsz, tid;
  int sxz_pml_I_a,   sxz_pml_I_aa;
  int sxz_pml_II_a,  sxz_pml_II_aa;
  int sxz_pml_III_a, sxz_pml_III_aa;
  int sxz_pml_IV_a,  sxz_pml_IV_aa;
  int sxz_pml_V_a,   sxz_pml_V_aa;
  int sxz_pml_VI_a,  sxz_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _sxz, * restrict _sxzend;
  register ireal * restrict _sxz_x_I,   * restrict _sxz_z_I; 
  register ireal * restrict _sxz_x_II,  * restrict _sxz_z_II;
  register ireal * restrict _sxz_x_III, * restrict _sxz_z_III;
  register ireal * restrict _sxz_x_IV,  * restrict _sxz_z_IV;
  register ireal * restrict _sxz_x_V,   * restrict _sxz_z_V;
  register ireal * restrict _sxz_x_VI,  * restrict _sxz_z_VI;
  register ireal * restrict _ms2;
  register ireal * restrict _vz3;
  register ireal * restrict _vx3, * restrict _vx2, * restrict _vx1, * restrict _vx0;
  register ireal * restrict _epx, * restrict _epz;
  register ireal lax, laz, dt2, vz3, vz2, vz1, vz0, dfdx, dfdz, etaxdt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_S2]._dims[0].n;
  ny = s[D_S2]._dims[1].n;
  nz = s[D_S2]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_S2]._dims[0].gs;
  gxe = s[D_S2]._dims[0].ge;
  gys = s[D_S2]._dims[1].gs;
  gye = s[D_S2]._dims[1].ge;
  gzs = s[D_S2]._dims[2].gs;
  gze = s[D_S2]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_S2,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_S2]._dims[2].gs;
    gze_pml_I = s_pml[D_S2]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_S2,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_S2]._dims[2].gs;
    gze_pml_II = s_pml[D_S2]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_S2,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_S2]._dims[1].gs;
    gye_pml_III = s_pml[D_S2]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_S2,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_S2]._dims[1].gs;
    gye_pml_IV = s_pml[D_S2]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_S2,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_S2]._dims[0].gs;
    gxe_pml_V = s_pml[D_S2]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_S2,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_S2]._dims[0].gs;
    gxe_pml_VI = s_pml[D_S2]._dims[0].ge;
  }
  
  sxz_a = s[D_S2]._dims[0].n0 - nx;

  ms2_a = s[D_MS2]._dims[0].n0 - nx;
  
  vx_a = s[D_V0]._dims[0].n0 - nx;
  vz_a = s[D_V2]._dims[0].n0 - nx;

  sxz_pml_I_a   = ld_pml[0]._s[D_S2]._dims[0].n0 - nx;
  sxz_pml_II_a  = ld_pml[3]._s[D_S2]._dims[0].n0 - nx;
  sxz_pml_III_a = ld_pml[6]._s[D_S2]._dims[0].n0 - nx;
  sxz_pml_IV_a  = ld_pml[9]._s[D_S2]._dims[0].n0 - nx;
  sxz_pml_V_a   = ld_pml[12]._s[D_S2]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  sxz_pml_VI_a  = ld_pml[15]._s[D_S2]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);  

#pragma omp parallel private( \
  tsz,tid,iy,_sxz,_sxzend,_ms2,_epx,_epz,_vz3,_vx3,_vx2,_vx1,_vx0,vz3,vz2,vz1,vz0, \
  _sxz_x_I,   _sxz_z_I,   \
  _sxz_x_II,  _sxz_z_II,  \
  _sxz_x_III, _sxz_z_III, \
  _sxz_x_IV,  _sxz_z_IV,  \
  _sxz_x_V,   _sxz_z_V,   \
  _sxz_x_VI,  _sxz_z_VI,  \
  dfdx,dfdz,etaxdt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      sxz_aa = (tsz * s[D_S2]._dims[1].n0 - ny) * s[D_S2]._dims[0].n0;
      
      ms2_aa = (tsz * s[D_MS2]._dims[1].n0 - ny) * s[D_MS2]._dims[0].n0;
       
      vx_aa = (tsz * s[D_V0]._dims[1].n0 - ny) * s[D_V0]._dims[0].n0;
      vz_aa = (tsz * s[D_V2]._dims[1].n0 - ny) * s[D_V2]._dims[0].n0;
      
      sxz_pml_I_aa   = (tsz * ld_pml[0]._s[D_S2]._dims[1].n0 - ny) * ld_pml[0]._s[D_S2]._dims[0].n0;
      sxz_pml_II_aa  = (tsz * ld_pml[3]._s[D_S2]._dims[1].n0 - ny) * ld_pml[3]._s[D_S2]._dims[0].n0;
      sxz_pml_III_aa = (tsz * ld_pml[6]._s[D_S2]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_S2]._dims[0].n0;
      sxz_pml_IV_aa  = (tsz * ld_pml[9]._s[D_S2]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_S2]._dims[0].n0;
      sxz_pml_V_aa   = (tsz * ld_pml[12]._s[D_S2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_S2]._dims[0].n0;
      sxz_pml_VI_aa  = (tsz * ld_pml[15]._s[D_S2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_S2]._dims[0].n0;  
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_S2]._dims;
    _sxz     = s[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[0]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_x_I  = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[2]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_z_I  = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
   
    dims    = s[D_MS2]._dims;
    _ms2    = cs[D_MS2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    dims    = s[D_V2]._dims;
    _vz3    = rs[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 2;
    
    dims    = s[D_V0]._dims;
    _vx1    = rs[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vx2    = _vx1 + dims[1].n0 * dims[0].n0;
    _vx3    = _vx2 + dims[1].n0 * dims[0].n0;
    _vx0    = _vx1 - dims[1].n0 * dims[0].n0;
    
    _epx    = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EV[2]]._s + (gzs_pml_I + tid - s[D_EV[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        vz2 = _vz3[-1]; vz1 = _vz3[-2]; vz0 = _vz3[-3];
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxzend = _sxz + nx; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_I) = ((*_sxz_x_I) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_I) = ((*_sxz_z_I) * (1.0f - etazdt) + dfdz*(*_ms2))/(1.0f + etazdt);
            *_sxz       = *_sxz_x_I + * _sxz_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_I++; _sxz_z_I++;
          
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxz += sxz_a; 
        _ms2 += ms2_a;
        _vz3 += vz_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a;
        _epx -= nx;
        _sxz_x_I += sxz_pml_I_a;
        _sxz_z_I += sxz_pml_I_a;
      }
      _sxz += sxz_aa;
      _ms2 += ms2_aa;
      _vz3 += vz_aa; 
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa;
      _epz += tsz;
      _sxz_x_I += sxz_pml_I_aa; 
      _sxz_z_I += sxz_pml_I_aa; 
    }
    /*************************************************************************************/

    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _sxz  += (gze_pml_I + 1 + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _ms2  += (gze_pml_I + 1 + tid - iz) * s[D_MP00]._dims[0].n0 * s[D_MP00]._dims[1].n0;
    _vz3  += (gze_pml_I + 1 + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vx2  += (gze_pml_I + 1 + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
      
    _vx3 = _vx2 + s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;

    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gze_pml_I + 1 + tid - s[D_EV[2]]._dims[0].gs);        /* 1D */
    
    s_pml      = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims       = s_pml[D_S2]._dims;
    _sxz_x_III = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml      = ld_pml[8]._s;
    dims       = s_pml[D_S2]._dims;
    _sxz_z_III = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S2]._dims;
    _sxz_x_IV = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_z_IV = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S2]._dims;
    _sxz_x_V  = s_pml[D_S2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
   
    s_pml     = ld_pml[14]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_z_V  = s_pml[D_S2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_S2]._dims;
    _sxz_x_VI = s_pml[D_S2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[17]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_z_VI = s_pml[D_S2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        vz2 = _vz3[-1]; vz1 = _vz3[-2]; vz0 = _vz3[-3];
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxzend = _sxz + nx; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_III) = ((*_sxz_x_III) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_III) =  (*_sxz_z_III) + dfdz*(*_ms2);
            *_sxz         = *_sxz_x_III + *_sxz_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_III++; _sxz_z_III++;
          
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        _sxz += sxz_a; 
        _ms2 += ms2_a;
        _vz3 += vz_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a;
        _epx -= nx;
        _sxz_x_III += sxz_pml_III_a;
        _sxz_z_III += sxz_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        vz2 = _vz3[-1]; vz1 = _vz3[-2]; vz0 = _vz3[-3];
        for ( _sxzend = _sxz + gxe_pml_V-gxs_pml_V+1; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_V) = ((*_sxz_x_V) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_V) =  (*_sxz_z_V) + dfdz*(*_ms2);
            *_sxz       = *_sxz_x_V + *_sxz_z_V;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_V++; _sxz_z_V++;
       
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        
        /** physical region */
        for ( _sxzend = _sxz + gxs_pml_VI-gxe_pml_V-1; _sxz < _sxzend;) {
          vz3 = *_vz3++;
          _epx ++;
          
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            *_sxz = *_sxz + (dfdx + dfdz) * (*_ms2);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          vz0 = vz1; vz1 = vz2; vz2 = vz3;          
        }
        
        /** pml region VI */
        for ( _sxzend = _sxz + gxe_pml_VI-gxs_pml_VI+1; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_VI) = ((*_sxz_x_VI) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_VI) =  (*_sxz_z_VI) + dfdz*(*_ms2);
            *_sxz        = *_sxz_x_VI + *_sxz_z_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_VI++; _sxz_z_VI++;
          
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxz += sxz_a; 
        _ms2 += ms2_a;
        _vz3 += vz_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; 
        _epx -= nx;
        _sxz_x_V += sxz_pml_V_a; 
        _sxz_z_V += sxz_pml_V_a; 
        _sxz_x_VI += sxz_pml_VI_a;
        _sxz_z_VI += sxz_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        vz2 = _vz3[-1]; vz1 = _vz3[-2]; vz0 = _vz3[-3];
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxzend = _sxz + nx; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_IV) = ((*_sxz_x_IV) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_IV) =  (*_sxz_z_IV) + dfdz*(*_ms2);
            *_sxz        = *_sxz_x_IV + *_sxz_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_IV++; _sxz_z_IV++;
          
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        _sxz += sxz_a; 
        _ms2 += ms2_a;
        _vz3 += vz_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a;
        _epx -= nx;
        _sxz_x_IV += sxz_pml_IV_a;
        _sxz_z_IV += sxz_pml_IV_a;
      }
      _sxz += sxz_aa;
      _ms2 += ms2_aa;
      _vz3 += vz_aa;
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa;
      
      _sxz_x_III += sxz_pml_III_aa;
      _sxz_z_III += sxz_pml_III_aa;

      _sxz_x_IV += sxz_pml_IV_aa;
      _sxz_z_IV += sxz_pml_IV_aa;

      _sxz_x_V += sxz_pml_V_aa;
      _sxz_z_V += sxz_pml_V_aa;

      _sxz_x_VI += sxz_pml_VI_aa;
      _sxz_z_VI += sxz_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _sxz  += (gzs_pml_II + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _ms2  += (gzs_pml_II + tid - iz) * s[D_MS2]._dims[0].n0 * s[D_MS2]._dims[1].n0;
    _vz3  += (gzs_pml_II + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _vx2  += (gzs_pml_II + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    
    _vx3 = _vx2 + s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
   
    s_pml     = ld_pml[3]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_x_II = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[5]._s;
    dims      = s_pml[D_S2]._dims;
    _sxz_z_II = s_pml[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gzs_pml_II + tid - s[D_EV[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        vz2 = _vz3[-1]; vz1 = _vz3[-2]; vz0 = _vz3[-3];
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _sxzend = _sxz + nx; _sxz < _sxzend; ) {
          vz3 = *_vz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (vz3 - vz0 + (vz1 - vz2) * 27.0) * lax; 
          dfdz = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * laz;
          if(_fwd) {
            (*_sxz_x_II) = ((*_sxz_x_II) * (1.0f - etaxdt) + dfdx*(*_ms2))/(1.0f + etaxdt);
            (*_sxz_z_II) = ((*_sxz_z_II) * (1.0f - etazdt) + dfdz*(*_ms2))/(1.0f + etazdt);
            *_sxz        = *_sxz_x_II + *_sxz_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _sxz++; _ms2++;
          _sxz_x_II++; _sxz_z_II++;
          
          vz0 = vz1; vz1 = vz2; vz2 = vz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _sxz += sxz_a; 
        _ms2 += ms2_a;
        _vz3 += vz_a; 
        _vx0 += vx_a; _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; 
        _epx -= nx;
        _sxz_x_II += sxz_pml_II_a;
        _sxz_z_II += sxz_pml_II_a;
      }
      _sxz += sxz_aa; 
      _ms2 += ms2_aa;
      _vz3 += vz_aa; 
      _vx0 += vx_aa; _vx1 += vx_aa; _vx2 += vx_aa; _vx3 += vx_aa;
      _epz += tsz;
      _sxz_x_II += sxz_pml_II_aa;
      _sxz_z_II += sxz_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** velocity component v0 (vx) */
int esgn_gts3d_24v0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, vx_a, mvx_a, px_a, sxy_a, sxz_a;
  int vx_aa, px_aa, sxy_aa, sxz_aa, mvx_aa, iy, iz, tsz, tid;
  int vx_pml_I_a,   vx_pml_I_aa;
  int vx_pml_II_a,  vx_pml_II_aa;
  int vx_pml_III_a, vx_pml_III_aa; 
  int vx_pml_IV_a,  vx_pml_IV_aa;
  int vx_pml_V_a,   vx_pml_V_aa;
  int vx_pml_VI_a,  vx_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _vx, * restrict _vxend;
  register ireal * restrict _vx_x_I,   * restrict _vx_y_I,   * restrict _vx_z_I; 
  register ireal * restrict _vx_x_II,  * restrict _vx_y_II,  * restrict _vx_z_II;
  register ireal * restrict _vx_x_III, * restrict _vx_y_III, * restrict _vx_z_III;
  register ireal * restrict _vx_x_IV,  * restrict _vx_y_IV,  * restrict _vx_z_IV;
  register ireal * restrict _vx_x_V,   * restrict _vx_y_V,   * restrict _vx_z_V;
  register ireal * restrict _vx_x_VI,  * restrict _vx_y_VI,  * restrict _vx_z_VI;
  register ireal * restrict _mvx;
  register ireal * restrict _px3;
  register ireal * restrict _sxy3, * restrict _sxy2, * restrict _sxy1, * restrict _sxy0;
  register ireal * restrict _sxz3, * restrict _sxz2, * restrict _sxz1, * restrict _sxz0;
  register ireal * restrict _epx, * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, dt2, px3, px2, px1, px0, dfdx, dfdy, dfdz, etaxdt, etaydt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  nz = s[D_V0]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_V0]._dims[0].gs;
  gxe = s[D_V0]._dims[0].ge;
  gys = s[D_V0]._dims[1].gs;
  gye = s[D_V0]._dims[1].ge;
  gzs = s[D_V0]._dims[2].gs;
  gze = s[D_V0]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_V0,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_V0]._dims[2].gs;
    gze_pml_I = s_pml[D_V0]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_V0,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_V0]._dims[2].gs;
    gze_pml_II = s_pml[D_V0]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_V0,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_V0]._dims[1].gs;
    gye_pml_III = s_pml[D_V0]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_V0,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_V0]._dims[1].gs;
    gye_pml_IV = s_pml[D_V0]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_V0,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_V0]._dims[0].gs;
    gxe_pml_V = s_pml[D_V0]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_V0,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_V0]._dims[0].gs;
    gxe_pml_VI = s_pml[D_V0]._dims[0].ge;
  }
  
  vx_a = s[D_V0]._dims[0].n0 - nx;
  
  mvx_a = s[D_MV0]._dims[0].n0 - nx;
  
  px_a  = s[D_P0]._dims[0].n0 - nx;
  sxy_a = s[D_S0]._dims[0].n0 - nx;
  sxz_a = s[D_S2]._dims[0].n0 - nx;
  
  vx_pml_I_a   = ld_pml[0]._s[D_V0]._dims[0].n0 - nx;
  vx_pml_II_a  = ld_pml[3]._s[D_V0]._dims[0].n0 - nx;
  vx_pml_III_a = ld_pml[6]._s[D_V0]._dims[0].n0 - nx;
  vx_pml_IV_a  = ld_pml[9]._s[D_V0]._dims[0].n0 - nx;
  vx_pml_V_a   = ld_pml[12]._s[D_V0]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  vx_pml_VI_a  = ld_pml[15]._s[D_V0]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);
  
#pragma omp parallel private                                            \
  (                                                                     \
   tsz,tid,iy,_vx,_vxend,_mvx,_epx,_epy,_epz,_px3,_sxy3,_sxy2,_sxy1,_sxy0,_sxz3,_sxz2,_sxz1,_sxz0,px3,px2,px1,px0, \
   _vx_x_I,  _vx_y_I,  _vx_z_I,                                         \
   _vx_x_II, _vx_y_II, _vx_z_II,                                        \
   _vx_x_III,_vx_y_III,_vx_z_III,                                       \
   _vx_x_IV, _vx_y_IV, _vx_z_IV,                                        \
   _vx_x_V,  _vx_y_V,  _vx_z_V,                                         \
   _vx_x_VI, _vx_y_VI, _vx_z_VI,                                        \
   dfdx,dfdy,dfdz,etaxdt,etaydt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      vx_aa  = (tsz * s[D_V0]._dims[1].n0 - ny) * s[D_V0]._dims[0].n0;
      mvx_aa = (tsz * s[D_MV0]._dims[1].n0 - ny) * s[D_MV0]._dims[0].n0;
      px_aa  = (tsz * s[D_P0]._dims[1].n0 - ny) * s[D_P0]._dims[0].n0;
      sxy_aa = (tsz * s[D_S0]._dims[1].n0 - ny) * s[D_S0]._dims[0].n0;
      sxz_aa = (tsz * s[D_S2]._dims[1].n0 - ny) * s[D_S2]._dims[0].n0;
      
      vx_pml_I_aa   = (tsz * ld_pml[0]._s[D_V0]._dims[1].n0 - ny) * ld_pml[0]._s[D_V0]._dims[0].n0;
      vx_pml_II_aa  = (tsz * ld_pml[3]._s[D_V0]._dims[1].n0 - ny) * ld_pml[3]._s[D_V0]._dims[0].n0;
      vx_pml_III_aa = (tsz * ld_pml[6]._s[D_V0]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_V0]._dims[0].n0;
      vx_pml_IV_aa  = (tsz * ld_pml[9]._s[D_V0]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_V0]._dims[0].n0;
      vx_pml_V_aa   = (tsz * ld_pml[12]._s[D_V0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_V0]._dims[0].n0;
      vx_pml_VI_aa  = (tsz * ld_pml[15]._s[D_V0]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_V0]._dims[0].n0;
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_V0]._dims;
    _vx     = s[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    dims    = s_pml[D_V0]._dims;
    _vx_x_I = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml   = ld_pml[1]._s;
    dims    = s_pml[D_V0]._dims;
    _vx_y_I = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[2]._s;
    dims    = s_pml[D_V0]._dims;
    _vx_z_I = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_MV0]._dims;
    _mvx    = cs[D_MV0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_P0]._dims;
    _px3    = rs[D_P0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 2;

    dims     = s[D_S0]._dims;
    _sxy2    = rs[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _sxy3    = _sxy2 + dims[0].n0;
    _sxy1    = _sxy2 - dims[0].n0;
    _sxy0    = _sxy1 - dims[0].n0;

    dims     = s[D_S2]._dims;
    _sxz2    = rs[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _sxz3    = _sxz2 + dims[1].n0 * dims[0].n0;
    _sxz1    = _sxz2 - dims[1].n0 * dims[0].n0;
    _sxz0    = _sxz1 - dims[1].n0 * dims[0].n0;


    _epx    = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EP[2]]._s + (gzs_pml_I + tid - s[D_EP[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vxend = _vx + nx; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_I) = ((*_vx_x_I) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_I) = ((*_vx_y_I) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
            (*_vx_z_I) = ((*_vx_z_I) * (1.0f - etazdt) + dfdz*(*_mvx))/(1.0f + etazdt);
            *_vx       = *_vx_x_I + *_vx_y_I + * _vx_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          _vx_x_I++; _vx_y_I++; _vx_z_I++; 
          
          px0 = px1; px1 = px2; px2 = px3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vx += vx_a;
        _mvx += mvx_a;
        _px3 += px_a; 
        _sxy0 += sxy_a; _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; 
        _sxz0 += sxz_a; _sxz1 += sxz_a; _sxz2 += sxz_a; _sxz3 += sxz_a;
        _epx -= nx;
        _vx_x_I += vx_pml_I_a;
        _vx_y_I += vx_pml_I_a;
        _vx_z_I += vx_pml_I_a;
      }
      _vx += vx_aa;
      _mvx += mvx_aa;
      _px3 += px_aa; 
      _sxy0 += sxy_aa; _sxy1 += sxy_aa; _sxy2 += sxy_aa; _sxy3 += sxy_aa; 
      _sxz0 += sxz_aa; _sxz1 += sxz_aa; _sxz2 += sxz_aa; _sxz3 += sxz_aa;
      _epy -= ny; _epz += tsz;
      _vx_x_I += vx_pml_I_aa;
      _vx_y_I += vx_pml_I_aa;
      _vx_z_I += vx_pml_I_aa;
    }
    /*************************************************************************************/
    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _vx   += (gze_pml_I + 1 + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _mvx  += (gze_pml_I + 1 + tid - iz) * s[D_MV0]._dims[0].n0 * s[D_MV0]._dims[1].n0;
    _px3  += (gze_pml_I + 1 + tid - iz) * s[D_P0]._dims[0].n0 * s[D_P0]._dims[1].n0;
    _sxy2 += (gze_pml_I + 1 + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _sxz2 += (gze_pml_I + 1 + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    
    _sxy3 = _sxy2 + s[D_S0]._dims[0].n0;
    _sxy1 = _sxy2 - s[D_S0]._dims[0].n0;
    _sxy0 = _sxy1 - s[D_S0]._dims[0].n0;

    _sxz3 = _sxz2 + s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _sxz1 = _sxz2 - s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _sxz0 = _sxz1 - s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;

    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gze_pml_I + 1 + tid - s[D_EP[2]]._dims[0].gs);        /* 1D */
    
    s_pml     = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V0]._dims;
    _vx_x_III = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[7]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_y_III = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[8]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_z_III = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V0]._dims;
    _vx_x_IV  = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[10]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_y_IV  = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_z_IV  = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V0]._dims;
    _vx_x_V   = s_pml[D_V0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[13]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_y_V   = s_pml[D_V0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[14]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_z_V   = s_pml[D_V0]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V0]._dims;
    _vx_x_VI  = s_pml[D_V0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[16]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_y_VI  = s_pml[D_V0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[17]._s;
    dims      = s_pml[D_V0]._dims;
    _vx_z_VI  = s_pml[D_V0]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vxend = _vx + nx; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_III) = ((*_vx_x_III) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_III) = ((*_vx_y_III) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
            (*_vx_z_III) =  (*_vx_z_III) + dfdz*(*_mvx);
            *_vx         = *_vx_x_III + *_vx_y_III + *_vx_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          _vx_x_III++; _vx_y_III++; _vx_z_III++;
          
          px0 = px1; px1 = px2; px2 = px3;
        }
        _vx += vx_a;
        _mvx += mvx_a;
        _px3 += px_a; 
        _sxy0 += sxy_a; _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a;
        _sxz0 += sxz_a; _sxz1 += sxz_a; _sxz2 += sxz_a; _sxz3 += sxz_a;
        _epx -= nx;
        _vx_x_III += vx_pml_III_a;
        _vx_y_III += vx_pml_III_a;
        _vx_z_III += vx_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
        _epy ++;
        for ( _vxend = _vx + gxe_pml_V-gxs_pml_V+1; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_V) = ((*_vx_x_V) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_V) =  (*_vx_y_V) + dfdy*(*_mvx);
            (*_vx_z_V) =  (*_vx_z_V) + dfdz*(*_mvx);
            *_vx       = *_vx_x_V + *_vx_y_V + *_vx_z_V;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++; 
          _vx_x_V++; _vx_y_V++; _vx_z_V++;
          
          px0 = px1; px1 = px2; px2 = px3;
        }
        
        /** physical region */
        for ( _vxend = _vx + gxs_pml_VI-gxe_pml_V-1; _vx < _vxend;) {
          px3 = *_px3++;
          _epx ++;
          
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            *_vx = *_vx + (dfdx + dfdy + dfdz) * (*_mvx);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          px0 = px1; px1 = px2; px2 = px3;          
        }
        
        /** pml region VI */
        for ( _vxend = _vx + gxe_pml_VI-gxs_pml_VI+1; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_VI) = ((*_vx_x_VI) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_VI) =  (*_vx_y_VI) + dfdy*(*_mvx);
            (*_vx_z_VI) =  (*_vx_z_VI) + dfdz*(*_mvx);
            *_vx        = *_vx_x_VI + *_vx_y_VI + *_vx_z_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          _vx_x_VI++; _vx_y_VI++; _vx_z_VI++;
        	        
          px0 = px1; px1 = px2; px2 = px3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vx += vx_a;
        _mvx += mvx_a;
        _px3 += px_a; 
        _sxy0 += sxy_a; _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; 
        _sxz0 += sxz_a; _sxz1 += sxz_a; _sxz2 += sxz_a; _sxz3 += sxz_a;
        _epx -= nx;
        _vx_x_V += vx_pml_V_a;
        _vx_y_V += vx_pml_V_a;
        _vx_z_V += vx_pml_V_a;
        _vx_x_VI += vx_pml_VI_a;
        _vx_y_VI += vx_pml_VI_a;
        _vx_z_VI += vx_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vxend = _vx + nx; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_IV) = ((*_vx_x_IV) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_IV) = ((*_vx_y_IV) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
            (*_vx_z_IV) =  (*_vx_z_IV) + dfdz*(*_mvx);
            *_vx        = *_vx_x_IV + *_vx_y_IV + *_vx_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          _vx_x_IV++; _vx_y_IV++; _vx_z_IV++;
       	
          px0 = px1; px1 = px2; px2 = px3;
        }
        _vx += vx_a;
        _mvx += mvx_a;
        _px3 += px_a; 
        _sxy0 += sxy_a; _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a;
        _sxz0 += sxz_a; _sxz1 += sxz_a; _sxz2 += sxz_a; _sxz3 += sxz_a;
        _epx -= nx;
        _vx_x_IV += vx_pml_IV_a;
        _vx_y_IV += vx_pml_IV_a;
        _vx_z_IV += vx_pml_IV_a;
      }
      _epy -= ny;
      _vx += vx_aa;
      _mvx += mvx_aa;
      _px3 += px_aa;
      _sxy0 += sxy_aa; _sxy1 += sxy_aa; _sxy2 += sxy_aa; _sxy3 += sxy_aa;
      _sxz0 += sxz_aa; _sxz1 += sxz_aa; _sxz2 += sxz_aa; _sxz3 += sxz_aa;
     
      _vx_x_III += vx_pml_III_aa;
      _vx_y_III += vx_pml_III_aa;
      _vx_z_III += vx_pml_III_aa;

      _vx_x_IV += vx_pml_IV_aa; 
      _vx_y_IV += vx_pml_IV_aa;
      _vx_z_IV += vx_pml_IV_aa;

      _vx_x_V += vx_pml_V_aa;
      _vx_y_V += vx_pml_V_aa;
      _vx_z_V += vx_pml_V_aa;

      _vx_x_VI += vx_pml_VI_aa;
      _vx_y_VI += vx_pml_VI_aa;
      _vx_z_VI += vx_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _vx   += (gzs_pml_II + tid - iz) * s[D_V0]._dims[0].n0 * s[D_V0]._dims[1].n0;
    _mvx  += (gzs_pml_II + tid - iz) * s[D_MV0]._dims[0].n0 * s[D_MV0]._dims[1].n0;
    _px3  += (gzs_pml_II + tid - iz) * s[D_P0]._dims[0].n0 * s[D_P0]._dims[1].n0;
    _sxy2 += (gzs_pml_II + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _sxz2 += (gzs_pml_II + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    
    _sxy3 = _sxy2 + s[D_S0]._dims[0].n0;
    _sxy1 = _sxy2 - s[D_S0]._dims[0].n0;
    _sxy0 = _sxy1 - s[D_S0]._dims[0].n0;

    _sxz3 = _sxz2 + s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _sxz1 = _sxz2 - s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _sxz0 = _sxz1 - s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
  
    s_pml    = ld_pml[3]._s;
    dims     = s_pml[D_V0]._dims;
    _vx_x_II = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[4]._s;
    dims     = s_pml[D_V0]._dims;
    _vx_y_II = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml    = ld_pml[5]._s;
    dims     = s_pml[D_V0]._dims;
    _vx_z_II = s_pml[D_V0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
      
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gzs_pml_II + tid - s[D_EP[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vxend = _vx + nx; _vx < _vxend; ) {
          px3 = *_px3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
          dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
          dfdz = ((*_sxz3++) - (*_sxz0++) + ((*_sxz1++) - (*_sxz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vx_x_II) = ((*_vx_x_II) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
            (*_vx_y_II) = ((*_vx_y_II) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
            (*_vx_z_II) = ((*_vx_z_II) * (1.0f - etazdt) + dfdz*(*_mvx))/(1.0f + etazdt);
            *_vx        = *_vx_x_II + *_vx_y_II + * _vx_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vx++; _mvx++;
          _vx_x_II++; _vx_y_II++; _vx_z_II++; 
                 
          px0 = px1; px1 = px2; px2 = px3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vx += vx_a;
        _mvx += mvx_a;
        _px3 += px_a; 
        _sxy0 += sxy_a; _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; 
        _sxz0 += sxz_a; _sxz1 += sxz_a; _sxz2 += sxz_a; _sxz3 += sxz_a;
        _epx -= nx;
        _vx_x_II += vx_pml_II_a;
        _vx_y_II += vx_pml_II_a;
        _vx_z_II += vx_pml_II_a;
      }
      _vx += vx_aa;
      _mvx += mvx_aa; 
      _px3 += px_aa; 
      _sxy0 += sxy_aa; _sxy1 += sxy_aa; _sxy2 += sxy_aa; _sxy3 += sxy_aa;
      _sxz0 += sxz_aa; _sxz1 += sxz_aa; _sxz2 += sxz_aa; _sxz3 += sxz_aa;
      _epy -= ny; _epz += tsz;
      _vx_x_II += vx_pml_II_aa;
      _vx_y_II += vx_pml_II_aa;
      _vx_z_II += vx_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** velocity component v1 (vy) */
int esgn_gts3d_24v1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, vy_a, mvy_a, sxy_a, py_a, syz_a;
  int vy_aa, py_aa, sxy_aa, syz_aa, mvy_aa, iy, iz, tsz, tid;
  int vy_pml_I_a,   vy_pml_I_aa;
  int vy_pml_II_a,  vy_pml_II_aa;
  int vy_pml_III_a, vy_pml_III_aa; 
  int vy_pml_IV_a,  vy_pml_IV_aa;
  int vy_pml_V_a,   vy_pml_V_aa;
  int vy_pml_VI_a,  vy_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _vy, * restrict _vyend;
  register ireal * restrict _vy_x_I,   * restrict _vy_y_I,   * restrict _vy_z_I; 
  register ireal * restrict _vy_x_II,  * restrict _vy_y_II,  * restrict _vy_z_II;
  register ireal * restrict _vy_x_III, * restrict _vy_y_III, * restrict _vy_z_III;
  register ireal * restrict _vy_x_IV,  * restrict _vy_y_IV,  * restrict _vy_z_IV;
  register ireal * restrict _vy_x_V,   * restrict _vy_y_V,   * restrict _vy_z_V;
  register ireal * restrict _vy_x_VI,  * restrict _vy_y_VI,  * restrict _vy_z_VI;
  register ireal * restrict _mvy;
  register ireal * restrict _sxy3;
  register ireal * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0;
  register ireal * restrict _syz3, * restrict _syz2, * restrict _syz1, * restrict _syz0;
  register ireal * restrict _epx, * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, dt2, sxy3, sxy2, sxy1, sxy0, dfdx, dfdy, dfdz, etaxdt, etaydt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  nz = s[D_V1]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_V1]._dims[0].gs;
  gxe = s[D_V1]._dims[0].ge;
  gys = s[D_V1]._dims[1].gs;
  gye = s[D_V1]._dims[1].ge;
  gzs = s[D_V1]._dims[2].gs;
  gze = s[D_V1]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_V1,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_V1]._dims[2].gs;
    gze_pml_I = s_pml[D_V1]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_V1,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_V1]._dims[2].gs;
    gze_pml_II = s_pml[D_V1]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_V1,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_V1]._dims[1].gs;
    gye_pml_III = s_pml[D_V1]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_V1,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_V1]._dims[1].gs;
    gye_pml_IV = s_pml[D_V1]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_V1,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_V1]._dims[0].gs;
    gxe_pml_V = s_pml[D_V1]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_V1,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_V1]._dims[0].gs;
    gxe_pml_VI = s_pml[D_V1]._dims[0].ge;
  }
  
  vy_a = s[D_V1]._dims[0].n0 - nx;
  
  mvy_a = s[D_MV1]._dims[0].n0 - nx;
  
  sxy_a = s[D_S0]._dims[0].n0 - nx;
  py_a  = s[D_P1]._dims[0].n0 - nx;
  syz_a = s[D_S1]._dims[0].n0 - nx;
  
  vy_pml_I_a   = ld_pml[0]._s[D_V1]._dims[0].n0 - nx;
  vy_pml_II_a  = ld_pml[3]._s[D_V1]._dims[0].n0 - nx;
  vy_pml_III_a = ld_pml[6]._s[D_V1]._dims[0].n0 - nx;
  vy_pml_IV_a  = ld_pml[9]._s[D_V1]._dims[0].n0 - nx;
  vy_pml_V_a   = ld_pml[12]._s[D_V1]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  vy_pml_VI_a  = ld_pml[15]._s[D_V1]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);
  
#pragma omp parallel private                                            \
  (                                                                     \
   tsz,tid,iy,_vy,_vyend,_mvy,_epx,_epy,_epz,_sxy3,_py3,_py2,_py1,_py0,_syz3,_syz2,_syz1,_syz0,sxy3,sxy2,sxy1,sxy0, \
   _vy_x_I,  _vy_y_I,  _vy_z_I,                                         \
   _vy_x_II, _vy_y_II, _vy_z_II,                                        \
   _vy_x_III,_vy_y_III,_vy_z_III,                                       \
   _vy_x_IV, _vy_y_IV, _vy_z_IV,                                        \
   _vy_x_V,  _vy_y_V,  _vy_z_V,                                         \
   _vy_x_VI, _vy_y_VI, _vy_z_VI,                                        \
   dfdx,dfdy,dfdz,etaxdt,etaydt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      vy_aa  = (tsz * s[D_V1]._dims[1].n0 - ny) * s[D_V1]._dims[0].n0;
      mvy_aa = (tsz * s[D_MV1]._dims[1].n0 - ny) * s[D_MV1]._dims[0].n0;
      sxy_aa = (tsz * s[D_S0]._dims[1].n0 - ny) * s[D_S0]._dims[0].n0;
      py_aa  = (tsz * s[D_P1]._dims[1].n0 - ny) * s[D_P1]._dims[0].n0;
      syz_aa = (tsz * s[D_S1]._dims[1].n0 - ny) * s[D_S1]._dims[0].n0;
      
      vy_pml_I_aa   = (tsz * ld_pml[0]._s[D_V1]._dims[1].n0 - ny) * ld_pml[0]._s[D_V1]._dims[0].n0;
      vy_pml_II_aa  = (tsz * ld_pml[3]._s[D_V1]._dims[1].n0 - ny) * ld_pml[3]._s[D_V1]._dims[0].n0;
      vy_pml_III_aa = (tsz * ld_pml[6]._s[D_V1]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_V1]._dims[0].n0;
      vy_pml_IV_aa  = (tsz * ld_pml[9]._s[D_V1]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_V1]._dims[0].n0;
      vy_pml_V_aa   = (tsz * ld_pml[12]._s[D_V1]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_V1]._dims[0].n0;
      vy_pml_VI_aa  = (tsz * ld_pml[15]._s[D_V1]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_V1]._dims[0].n0;
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_V1]._dims;
    _vy     = s[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    dims    = s_pml[D_V1]._dims;
    _vy_x_I = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml   = ld_pml[1]._s;
    dims    = s_pml[D_V1]._dims;
    _vy_y_I = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[2]._s;
    dims    = s_pml[D_V1]._dims;
    _vy_z_I = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_MV1]._dims;
    _mvy    = cs[D_MV1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_S0]._dims;
    _sxy3   = rs[D_S0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 1;

    dims     = s[D_P1]._dims;
    _py1     = rs[D_P1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _py2     = _py1 + dims[0].n0;
    _py3     = _py2 + dims[0].n0;
    _py0     = _py1 - dims[0].n0;

    dims     = s[D_S1]._dims;
    _syz2    = rs[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _syz3    = _syz2 + dims[1].n0 * dims[0].n0;
    _syz1    = _syz2 - dims[1].n0 * dims[0].n0;
    _syz0    = _syz1 - dims[1].n0 * dims[0].n0;


    _epx    = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EP[2]]._s + (gzs_pml_I + tid - s[D_EP[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vyend = _vy + nx; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_I) = ((*_vy_x_I) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_I) = ((*_vy_y_I) * (1.0f - etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);
            (*_vy_z_I) = ((*_vy_z_I) * (1.0f - etazdt) + dfdz*(*_mvy))/(1.0f + etazdt);
            *_vy       = *_vy_x_I + *_vy_y_I + * _vy_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          _vy_x_I++; _vy_y_I++; _vy_z_I++; 
          
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vy += vy_a;
        _mvy += mvy_a;
        _sxy3 += sxy_a; 
        _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _epx -= nx;
        _vy_x_I += vy_pml_I_a;
        _vy_y_I += vy_pml_I_a;
        _vy_z_I += vy_pml_I_a;
      }
      _vy += vy_aa;
      _mvy += mvy_aa;
      _sxy3 += sxy_aa; 
      _py0 += py_aa; _py1 += py_aa; _py2 += py_aa; _py3 += py_aa; 
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa;
      _epy -= ny; _epz += tsz;
      _vy_x_I += vy_pml_I_aa;
      _vy_y_I += vy_pml_I_aa;
      _vy_z_I += vy_pml_I_aa;
    }
    /*************************************************************************************/
    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _vy   += (gze_pml_I + 1 + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _mvy  += (gze_pml_I + 1 + tid - iz) * s[D_MV1]._dims[0].n0 * s[D_MV1]._dims[1].n0;
    _sxy3 += (gze_pml_I + 1 + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _py2  += (gze_pml_I + 1 + tid - iz) * s[D_P1]._dims[0].n0 * s[D_P1]._dims[1].n0;
    _syz2 += (gze_pml_I + 1 + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    
    _py3  = _py2 + s[D_P1]._dims[0].n0;
    _py1  = _py2 - s[D_P1]._dims[0].n0;
    _py0  = _py1 - s[D_P1]._dims[0].n0;

    _syz3 = _syz2 + s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _syz1 = _syz2 - s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _syz0 = _syz1 - s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;

    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gze_pml_I + 1 + tid - s[D_EP[2]]._dims[0].gs);        /* 1D */
    
    s_pml     = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V1]._dims;
    _vy_x_III = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[7]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_y_III = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[8]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_z_III = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V1]._dims;
    _vy_x_IV  = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[10]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_y_IV  = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_z_IV  = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V1]._dims;
    _vy_x_V   = s_pml[D_V1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[13]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_y_V   = s_pml[D_V1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[14]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_z_V   = s_pml[D_V1]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V1]._dims;
    _vy_x_VI  = s_pml[D_V1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[16]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_y_VI  = s_pml[D_V1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[17]._s;
    dims      = s_pml[D_V1]._dims;
    _vy_z_VI  = s_pml[D_V1]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vyend = _vy + nx; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_III) = ((*_vy_x_III) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_III) = ((*_vy_y_III) * (1.0f - etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);
            (*_vy_z_III) =  (*_vy_z_III) + dfdz*(*_mvy);
            *_vy         = *_vy_x_III + *_vy_y_III + *_vy_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          _vy_x_III++; _vy_y_III++; _vy_z_III++;
          
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        _vy += vy_a;
        _mvy += mvy_a;
        _sxy3 += sxy_a; 
        _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a;
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _epx -= nx;
        _vy_x_III += vy_pml_III_a;
        _vy_y_III += vy_pml_III_a;
        _vy_z_III += vy_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
        _epy ++;
        for ( _vyend = _vy + gxe_pml_V-gxs_pml_V+1; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_V) = ((*_vy_x_V) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_V) =  (*_vy_y_V) + dfdy*(*_mvy);
            (*_vy_z_V) =  (*_vy_z_V) + dfdz*(*_mvy);
            *_vy       = *_vy_x_V + *_vy_y_V + *_vy_z_V;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++; 
          _vy_x_V++; _vy_y_V++; _vy_z_V++;
          
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        
        /** physical region */
        for ( _vyend = _vy + gxs_pml_VI-gxe_pml_V-1; _vy < _vyend;) {
          sxy3 = *_sxy3++;
          _epx ++;
          
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            *_vy = *_vy + (dfdx + dfdy + dfdz) * (*_mvy);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;          
        }
        
        /** pml region VI */
        for ( _vyend = _vy + gxe_pml_VI-gxs_pml_VI+1; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_VI) = ((*_vy_x_VI) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_VI) =  (*_vy_y_VI) + dfdy*(*_mvy);
            (*_vy_z_VI) =  (*_vy_z_VI) + dfdz*(*_mvy);
            *_vy        = *_vy_x_VI + *_vy_y_VI + *_vy_z_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          _vy_x_VI++; _vy_y_VI++; _vy_z_VI++;
        	        
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vy += vy_a;
        _mvy += mvy_a;
        _sxy3 += sxy_a; 
        _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _epx -= nx;
        _vy_x_V += vy_pml_V_a;
        _vy_y_V += vy_pml_V_a;
        _vy_z_V += vy_pml_V_a;
        _vy_x_VI += vy_pml_VI_a;
        _vy_y_VI += vy_pml_VI_a;
        _vy_z_VI += vy_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vyend = _vy + nx; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_IV) = ((*_vy_x_IV) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_IV) = ((*_vy_y_IV) * (1.0f - etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);
            (*_vy_z_IV) =  (*_vy_z_IV) + dfdz*(*_mvy);
            *_vy        = *_vy_x_IV + *_vy_y_IV + *_vy_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          _vy_x_IV++; _vy_y_IV++; _vy_z_IV++;
       	
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        _vy += vy_a;
        _mvy += mvy_a;
        _sxy3 += sxy_a; 
        _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a;
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _epx -= nx;
        _vy_x_IV += vy_pml_IV_a;
        _vy_y_IV += vy_pml_IV_a;
        _vy_z_IV += vy_pml_IV_a;
      }
      _epy -= ny;
      _vy += vy_aa;
      _mvy += mvy_aa;
      _sxy3 += sxy_aa;
      _py0 += py_aa; _py1 += py_aa; _py2 += py_aa; _py3 += py_aa;
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa;
     
      _vy_x_III += vy_pml_III_aa;
      _vy_y_III += vy_pml_III_aa;
      _vy_z_III += vy_pml_III_aa;

      _vy_x_IV += vy_pml_IV_aa; 
      _vy_y_IV += vy_pml_IV_aa;
      _vy_z_IV += vy_pml_IV_aa;

      _vy_x_V += vy_pml_V_aa;
      _vy_y_V += vy_pml_V_aa;
      _vy_z_V += vy_pml_V_aa;

      _vy_x_VI += vy_pml_VI_aa;
      _vy_y_VI += vy_pml_VI_aa;
      _vy_z_VI += vy_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _vy   += (gzs_pml_II + tid - iz) * s[D_V1]._dims[0].n0 * s[D_V1]._dims[1].n0;
    _mvy  += (gzs_pml_II + tid - iz) * s[D_MV1]._dims[0].n0 * s[D_MV1]._dims[1].n0;
    _sxy3 += (gzs_pml_II + tid - iz) * s[D_S0]._dims[0].n0 * s[D_S0]._dims[1].n0;
    _py2  += (gzs_pml_II + tid - iz) * s[D_P1]._dims[0].n0 * s[D_P1]._dims[1].n0;
    _syz2 += (gzs_pml_II + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    
    _py3 = _py2 + s[D_P1]._dims[0].n0;
    _py1 = _py2 - s[D_P1]._dims[0].n0;
    _py0 = _py1 - s[D_P1]._dims[0].n0;

    _syz3 = _syz2 + s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _syz1 = _syz2 - s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _syz0 = _syz1 - s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
  
    s_pml    = ld_pml[3]._s;
    dims     = s_pml[D_V1]._dims;
    _vy_x_II = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[4]._s;
    dims     = s_pml[D_V1]._dims;
    _vy_y_II = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml    = ld_pml[5]._s;
    dims     = s_pml[D_V1]._dims;
    _vy_z_II = s_pml[D_V1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
      
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EP[2]]._s + (gzs_pml_II + tid - s[D_EP[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vyend = _vy + nx; _vy < _vyend; ) {
          sxy3 = *_sxy3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
          dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
          dfdz = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vy_x_II) = ((*_vy_x_II) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
            (*_vy_y_II) = ((*_vy_y_II) * (1.0f - etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);
            (*_vy_z_II) = ((*_vy_z_II) * (1.0f - etazdt) + dfdz*(*_mvy))/(1.0f + etazdt);
            *_vy        = *_vy_x_II + *_vy_y_II + * _vy_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vy++; _mvy++;
          _vy_x_II++; _vy_y_II++; _vy_z_II++; 
                 
          sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vy += vy_a;
        _mvy += mvy_a;
        _sxy3 += sxy_a; 
        _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _epx -= nx;
        _vy_x_II += vy_pml_II_a;
        _vy_y_II += vy_pml_II_a;
        _vy_z_II += vy_pml_II_a;
      }
      _vy += vy_aa;
      _mvy += mvy_aa; 
      _sxy3 += sxy_aa; 
      _py0 += py_aa; _py1 += py_aa; _py2 += py_aa; _py3 += py_aa;
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa;
      _epy -= ny; _epz += tsz;
      _vy_x_II += vy_pml_II_aa;
      _vy_y_II += vy_pml_II_aa;
      _vy_z_II += vy_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** velocity component v2 (vz) */
int esgn_gts3d_24v2(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, nz, gxs, gys, gzs, gxe, gye, gze, vz_a, mvz_a, pz_a, sxz_a, syz_a;
  int vz_aa, pz_aa, sxz_aa, syz_aa, mvz_aa, iy, iz, tsz, tid;
  int vz_pml_I_a,   vz_pml_I_aa;
  int vz_pml_II_a,  vz_pml_II_aa;
  int vz_pml_III_a, vz_pml_III_aa; 
  int vz_pml_IV_a,  vz_pml_IV_aa;
  int vz_pml_V_a,   vz_pml_V_aa;
  int vz_pml_VI_a,  vz_pml_VI_aa;
  int gzs_pml_I, gze_pml_I, gzs_pml_II, gze_pml_II;
  int gys_pml_III, gye_pml_III, gys_pml_IV, gye_pml_IV;
  int gxs_pml_V, gxe_pml_V, gxs_pml_VI, gxe_pml_VI;
  register ireal * restrict _vz, * restrict _vzend;
  register ireal * restrict _vz_x_I,   * restrict _vz_y_I,   * restrict _vz_z_I; 
  register ireal * restrict _vz_x_II,  * restrict _vz_y_II,  * restrict _vz_z_II;
  register ireal * restrict _vz_x_III, * restrict _vz_y_III, * restrict _vz_z_III;
  register ireal * restrict _vz_x_IV,  * restrict _vz_y_IV,  * restrict _vz_z_IV;
  register ireal * restrict _vz_x_V,   * restrict _vz_y_V,   * restrict _vz_z_V;
  register ireal * restrict _vz_x_VI,  * restrict _vz_y_VI,  * restrict _vz_z_VI;
  register ireal * restrict _mvz;
  register ireal * restrict _sxz3;
  register ireal * restrict _syz3, * restrict _syz2, * restrict _syz1, * restrict _syz0;
  register ireal * restrict _pz3, * restrict _pz2, * restrict _pz1, * restrict _pz0;
  register ireal * restrict _epx, * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, dt2, sxz3, sxz2, sxz1, sxz0, dfdx, dfdy, dfdz, etaxdt, etaydt, etazdt;

  /*  register ireal *tmp;  pointer for swapping p with mp during imaging accumulation*/ 
  /* register ireal *_rmpx;  pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;
  INFODIM * dims;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_V2]._dims[0].n;
  ny = s[D_V2]._dims[1].n;
  nz = s[D_V2]._dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
  
  gxs = s[D_V2]._dims[0].gs;
  gxe = s[D_V2]._dims[0].ge;
  gys = s[D_V2]._dims[1].gs;
  gye = s[D_V2]._dims[1].ge;
  gzs = s[D_V2]._dims[2].gs;
  gze = s[D_V2]._dims[2].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  laz = ((ESGN_TS_PARS*)pars)->lam[2] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_V2,&empty);
  if (empty) {
    gzs_pml_I = gzs;
    gze_pml_I = gzs-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gzs_pml_I = s_pml[D_V2]._dims[2].gs;
    gze_pml_I = s_pml[D_V2]._dims[2].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+3,D_V2,&empty);
  if (empty) {
    gzs_pml_II = gze+1;
    gze_pml_II = gze;
  }
  else {
    s_pml = ld_pml[3]._s;
    gzs_pml_II = s_pml[D_V2]._dims[2].gs;
    gze_pml_II = s_pml[D_V2]._dims[2].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+6,D_V2,&empty);
  if (empty) {
    gys_pml_III = gys;
    gye_pml_III = gys-1;
  }
  else {
    s_pml = ld_pml[6]._s;
    gys_pml_III = s_pml[D_V2]._dims[1].gs;
    gye_pml_III = s_pml[D_V2]._dims[1].ge;
  }
  /** pml region IV */
  rd_empty(ld_pml+9,D_V2,&empty);
  if (empty) {
    gys_pml_IV = gye+1;
    gye_pml_IV = gye;
  }
  else {
    s_pml = ld_pml[9]._s;
    gys_pml_IV = s_pml[D_V2]._dims[1].gs;
    gye_pml_IV = s_pml[D_V2]._dims[1].ge;
  }
  /** pml region V */
  rd_empty(ld_pml+12,D_V2,&empty);
  if (empty) {
    gxs_pml_V = gxs;
    gxe_pml_V = gxs-1;
  }
  else {
    s_pml = ld_pml[12]._s;
    gxs_pml_V = s_pml[D_V2]._dims[0].gs;
    gxe_pml_V = s_pml[D_V2]._dims[0].ge;
  }
  /** pml region VI */
  rd_empty(ld_pml+15,D_V2,&empty);
  if (empty) {
    gxs_pml_VI = gxe+1;
    gxe_pml_VI = gxe;
  }
  else {
    s_pml = ld_pml[15]._s;
    gxs_pml_VI = s_pml[D_V2]._dims[0].gs;
    gxe_pml_VI = s_pml[D_V2]._dims[0].ge;
  }
  
  vz_a = s[D_V2]._dims[0].n0 - nx;
  
  mvz_a = s[D_MV2]._dims[0].n0 - nx;
  
  sxz_a = s[D_S2]._dims[0].n0 - nx;
  syz_a = s[D_S1]._dims[0].n0 - nx;
  pz_a  = s[D_P2]._dims[0].n0 - nx;
  
  vz_pml_I_a   = ld_pml[0]._s[D_V2]._dims[0].n0 - nx;
  vz_pml_II_a  = ld_pml[3]._s[D_V2]._dims[0].n0 - nx;
  vz_pml_III_a = ld_pml[6]._s[D_V2]._dims[0].n0 - nx;
  vz_pml_IV_a  = ld_pml[9]._s[D_V2]._dims[0].n0 - nx;
  vz_pml_V_a   = ld_pml[12]._s[D_V2]._dims[0].n0 - (gxe_pml_V - gxs_pml_V + 1);
  vz_pml_VI_a  = ld_pml[15]._s[D_V2]._dims[0].n0 - (gxe_pml_VI - gxs_pml_VI + 1);
  
#pragma omp parallel private                                            \
  (                                                                     \
   tsz,tid,iy,_vz,_vzend,_mvz,_epx,_epy,_epz,_sxz3,_syz3,_syz2,_syz1,_syz0,_pz3,_pz2,_pz1,_pz0,sxz3,sxz2,sxz1,sxz0, \
   _vz_x_I,  _vz_y_I,  _vz_z_I,                                         \
   _vz_x_II, _vz_y_II, _vz_z_II,                                        \
   _vz_x_III,_vz_y_III,_vz_z_III,                                       \
   _vz_x_IV, _vz_y_IV, _vz_z_IV,                                        \
   _vz_x_V,  _vz_y_V,  _vz_z_V,                                         \
   _vz_x_VI, _vz_y_VI, _vz_z_VI,                                        \
   dfdx,dfdy,dfdz,etaxdt,etaydt,etazdt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

#pragma omp single
    {    
      vz_aa  = (tsz * s[D_V2]._dims[1].n0 - ny) * s[D_V2]._dims[0].n0;
      mvz_aa = (tsz * s[D_MV2]._dims[1].n0 - ny) * s[D_MV2]._dims[0].n0;
      sxz_aa = (tsz * s[D_S2]._dims[1].n0 - ny) * s[D_S2]._dims[0].n0;
      syz_aa = (tsz * s[D_S1]._dims[1].n0 - ny) * s[D_S1]._dims[0].n0;
      pz_aa  = (tsz * s[D_P2]._dims[1].n0 - ny) * s[D_P2]._dims[0].n0;
      
      vz_pml_I_aa   = (tsz * ld_pml[0]._s[D_V2]._dims[1].n0 - ny) * ld_pml[0]._s[D_V2]._dims[0].n0;
      vz_pml_II_aa  = (tsz * ld_pml[3]._s[D_V2]._dims[1].n0 - ny) * ld_pml[3]._s[D_V2]._dims[0].n0;
      vz_pml_III_aa = (tsz * ld_pml[6]._s[D_V2]._dims[1].n0 - (gye_pml_III - gys_pml_III + 1)) * ld_pml[6]._s[D_V2]._dims[0].n0;
      vz_pml_IV_aa  = (tsz * ld_pml[9]._s[D_V2]._dims[1].n0 - (gye_pml_IV - gys_pml_IV + 1)) * ld_pml[9]._s[D_V2]._dims[0].n0;
      vz_pml_V_aa   = (tsz * ld_pml[12]._s[D_V2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[12]._s[D_V2]._dims[0].n0;
      vz_pml_VI_aa  = (tsz * ld_pml[15]._s[D_V2]._dims[1].n0 - (gys_pml_IV - gye_pml_III - 1)) * ld_pml[15]._s[D_V2]._dims[0].n0;
    }
#pragma omp barrier
    /** pml region I *********************************************************************/
    /** gzs == gzs_pml_I */
    dims    = s[D_V2]._dims;
    _vz     = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    dims    = s_pml[D_V2]._dims;
    _vz_x_I = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml   = ld_pml[1]._s;
    dims    = s_pml[D_V2]._dims;
    _vz_y_I = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml   = ld_pml[2]._s;
    dims    = s_pml[D_V2]._dims;
    _vz_z_I = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims    = s[D_MV2]._dims;
    _mvz    = cs[D_MV2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
       
    dims     = s[D_S2]._dims;
    _sxz3    = rs[D_S2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 1;

    dims     = s[D_S1]._dims;
    _syz2    = rs[D_S1]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _syz3    = _syz2 + dims[0].n0;
    _syz1    = _syz2 - dims[0].n0;
    _syz0    = _syz1 - dims[0].n0;

    dims     = s[D_P2]._dims;
    _pz1    = rs[D_P2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_I + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _pz2    = _pz1 + dims[1].n0 * dims[0].n0;
    _pz3    = _pz2 + dims[1].n0 * dims[0].n0;
    _pz0    = _pz1 - dims[1].n0 * dims[0].n0;


    _epx    = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);              /* 1D */
    _epz    = s[D_EV[2]]._s + (gzs_pml_I + tid - s[D_EV[2]]._dims[0].gs);  /* 1D */
    
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_I+tid; iz < gze_pml_I+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++ ) {
        sxz2 = _sxz3[-1]; sxz1 = _sxz3[-2]; sxz0 = _sxz3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vzend = _vz + nx; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_I) = ((*_vz_x_I) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_I) = ((*_vz_y_I) * (1.0f - etaydt) + dfdy*(*_mvz))/(1.0f + etaydt);
            (*_vz_z_I) = ((*_vz_z_I) * (1.0f - etazdt) + dfdz*(*_mvz))/(1.0f + etazdt);
            *_vz       = *_vz_x_I + *_vz_y_I + * _vz_z_I;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          _vz_x_I++; _vz_y_I++; _vz_z_I++; 
          
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vz += vz_a;
        _mvz += mvz_a;
        _sxz3 += sxz_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a; 
        _pz0 += pz_a; _pz1 += pz_a; _pz2 += pz_a; _pz3 += pz_a;
        _epx -= nx;
        _vz_x_I += vz_pml_I_a;
        _vz_y_I += vz_pml_I_a;
        _vz_z_I += vz_pml_I_a;
      }
      _vz += vz_aa;
      _mvz += mvz_aa;
      _sxz3 += sxz_aa; 
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa; 
      _pz0 += pz_aa; _pz1 += pz_aa; _pz2 += pz_aa; _pz3 += pz_aa;
      _epy -= ny; _epz += tsz;
      _vz_x_I += vz_pml_I_aa;
      _vz_y_I += vz_pml_I_aa;
      _vz_z_I += vz_pml_I_aa;
    }
    /*************************************************************************************/
    /** pml region III, IV, V, VI and physical region *******************************************/
    /** adjust pointers */
    _vz   += (gze_pml_I + 1 + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _mvz  += (gze_pml_I + 1 + tid - iz) * s[D_MV2]._dims[0].n0 * s[D_MV2]._dims[1].n0;
    _sxz3 += (gze_pml_I + 1 + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _syz2 += (gze_pml_I + 1 + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _pz2  += (gze_pml_I + 1 + tid - iz) * s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    
    _syz3 = _syz2 + s[D_S1]._dims[0].n0;
    _syz1 = _syz2 - s[D_S1]._dims[0].n0;
    _syz0 = _syz1 - s[D_S1]._dims[0].n0;

    _pz3 = _pz2 + s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _pz1 = _pz2 - s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _pz0 = _pz1 - s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;

    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gze_pml_I + 1 + tid - s[D_EV[2]]._dims[0].gs);        /* 1D */
    
    s_pml     = ld_pml[6]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V2]._dims;
    _vz_x_III = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[7]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_y_III = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[8]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_z_III = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_III - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[9]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V2]._dims;
    _vz_x_IV  = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[10]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_y_IV  = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[11]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_z_IV  = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys_pml_IV - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
      
    s_pml     = ld_pml[12]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V2]._dims;
    _vz_x_V   = s_pml[D_V2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
     
    s_pml     = ld_pml[13]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_y_V   = s_pml[D_V2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
        
    s_pml     = ld_pml[14]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_z_V   = s_pml[D_V2]._s + (gxs_pml_V - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[15]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    dims      = s_pml[D_V2]._dims;
    _vz_x_VI  = s_pml[D_V2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
       
    s_pml     = ld_pml[16]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_y_VI  = s_pml[D_V2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    s_pml     = ld_pml[17]._s;
    dims      = s_pml[D_V2]._dims;
    _vz_z_VI  = s_pml[D_V2]._s + (gxs_pml_VI - dims[0].gs) + (gye_pml_III + 1 - dims[1].gs + (gze_pml_I + 1 - dims[2].gs + tid) * dims[1].n0) * dims[0].n0;
    
    /** adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    
    for ( iz = gze_pml_I+1+tid; iz < gzs_pml_II; iz += tsz) { 
      /** pml region III */
      for ( iy = gys_pml_III; iy < gye_pml_III + 1; iy ++) {
        sxz2 = _sxz3[-1]; sxz1 = _sxz3[-2]; sxz0 = _sxz3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vzend = _vz + nx; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_III) = ((*_vz_x_III) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_III) = ((*_vz_y_III) * (1.0f - etaydt) + dfdy*(*_mvz))/(1.0f + etaydt);
            (*_vz_z_III) =  (*_vz_z_III) + dfdz*(*_mvz);
            *_vz         = *_vz_x_III + *_vz_y_III + *_vz_z_III;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          _vz_x_III++; _vz_y_III++; _vz_z_III++;
          
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        _vz += vz_a;
        _mvz += mvz_a;
        _sxz3 += sxz_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _pz0 += pz_a; _pz1 += pz_a; _pz2 += pz_a; _pz3 += pz_a;
        _epx -= nx;
        _vz_x_III += vz_pml_III_a;
        _vz_y_III += vz_pml_III_a;
        _vz_z_III += vz_pml_III_a;
      }
      
      /** pml region V, physical region and pml region VI*/
      for (iy = gye_pml_III+1; iy < gys_pml_IV; iy ++) {
        /** pml region V */
        sxz2 = _sxz3[-1]; sxz1 = _sxz3[-2]; sxz0 = _sxz3[-3];
        _epy ++;
        for ( _vzend = _vz + gxe_pml_V-gxs_pml_V+1; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_V) = ((*_vz_x_V) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_V) =  (*_vz_y_V) + dfdy*(*_mvz);
            (*_vz_z_V) =  (*_vz_z_V) + dfdz*(*_mvz);
            *_vz       = *_vz_x_V + *_vz_y_V + *_vz_z_V;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++; 
          _vz_x_V++; _vz_y_V++; _vz_z_V++;
          
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        
        /** physical region */
        for ( _vzend = _vz + gxs_pml_VI-gxe_pml_V-1; _vz < _vzend;) {
          sxz3 = *_sxz3++;
          _epx ++;
          
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            *_vz = *_vz + (dfdx + dfdy + dfdz) * (*_mvz);
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;          
        }
        
        /** pml region VI */
        for ( _vzend = _vz + gxe_pml_VI-gxs_pml_VI+1; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_VI) = ((*_vz_x_VI) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_VI) =  (*_vz_y_VI) + dfdy*(*_mvz);
            (*_vz_z_VI) =  (*_vz_z_VI) + dfdz*(*_mvz);
            *_vz        = *_vz_x_VI + *_vz_y_VI + *_vz_z_VI;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          _vz_x_VI++; _vz_y_VI++; _vz_z_VI++;
        	        
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vz += vz_a;
        _mvz += mvz_a;
        _sxz3 += sxz_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a; 
        _pz0 += pz_a; _pz1 += pz_a; _pz2 += pz_a; _pz3 += pz_a;
        _epx -= nx;
        _vz_x_V += vz_pml_V_a;
        _vz_y_V += vz_pml_V_a;
        _vz_z_V += vz_pml_V_a;
        _vz_x_VI += vz_pml_VI_a;
        _vz_y_VI += vz_pml_VI_a;
        _vz_z_VI += vz_pml_VI_a;
      }
      /** pml region IV */
      for ( iy = gys_pml_IV; iy < gye_pml_IV + 1; iy ++) {
        sxz2 = _sxz3[-1]; sxz1 = _sxz3[-2]; sxz0 = _sxz3[-3];
        etaydt = (*_epy ++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vzend = _vz + nx; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_IV) = ((*_vz_x_IV) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_IV) = ((*_vz_y_IV) * (1.0f - etaydt) + dfdy*(*_mvz))/(1.0f + etaydt);
            (*_vz_z_IV) =  (*_vz_z_IV) + dfdz*(*_mvz);
            *_vz        = *_vz_x_IV + *_vz_y_IV + *_vz_z_IV;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          _vz_x_IV++; _vz_y_IV++; _vz_z_IV++;
       	
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        _vz += vz_a;
        _mvz += mvz_a;
        _sxz3 += sxz_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a;
        _pz0 += pz_a; _pz1 += pz_a; _pz2 += pz_a; _pz3 += pz_a;
        _epx -= nx;
        _vz_x_IV += vz_pml_IV_a;
        _vz_y_IV += vz_pml_IV_a;
        _vz_z_IV += vz_pml_IV_a;
      }
      _epy -= ny;
      _vz += vz_aa;
      _mvz += mvz_aa;
      _sxz3 += sxz_aa;
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa;
      _pz0 += pz_aa; _pz1 += pz_aa; _pz2 += pz_aa; _pz3 += pz_aa;
     
      _vz_x_III += vz_pml_III_aa;
      _vz_y_III += vz_pml_III_aa;
      _vz_z_III += vz_pml_III_aa;

      _vz_x_IV += vz_pml_IV_aa; 
      _vz_y_IV += vz_pml_IV_aa;
      _vz_z_IV += vz_pml_IV_aa;

      _vz_x_V += vz_pml_V_aa;
      _vz_y_V += vz_pml_V_aa;
      _vz_z_V += vz_pml_V_aa;

      _vz_x_VI += vz_pml_VI_aa;
      _vz_y_VI += vz_pml_VI_aa;
      _vz_z_VI += vz_pml_VI_aa;
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _vz   += (gzs_pml_II + tid - iz) * s[D_V2]._dims[0].n0 * s[D_V2]._dims[1].n0;
    _mvz  += (gzs_pml_II + tid - iz) * s[D_MV2]._dims[0].n0 * s[D_MV2]._dims[1].n0;
    _sxz3 += (gzs_pml_II + tid - iz) * s[D_S2]._dims[0].n0 * s[D_S2]._dims[1].n0;
    _syz2 += (gzs_pml_II + tid - iz) * s[D_S1]._dims[0].n0 * s[D_S1]._dims[1].n0;
    _pz2  += (gzs_pml_II + tid - iz) * s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    
    _syz3 = _syz2 + s[D_S1]._dims[0].n0;
    _syz1 = _syz2 - s[D_S1]._dims[0].n0;
    _syz0 = _syz1 - s[D_S1]._dims[0].n0;

    _pz3 = _pz2 + s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _pz1 = _pz2 - s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
    _pz0 = _pz1 - s[D_P2]._dims[0].n0 * s[D_P2]._dims[1].n0;
  
    s_pml    = ld_pml[3]._s;
    dims     = s_pml[D_V2]._dims;
    _vz_x_II = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    
    s_pml    = ld_pml[4]._s;
    dims     = s_pml[D_V2]._dims;
    _vz_y_II = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
    s_pml    = ld_pml[5]._s;
    dims     = s_pml[D_V2]._dims;
    _vz_z_II = s_pml[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (gzs_pml_II + tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
      
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs);                        /* 1D */
    _epz = s[D_EV[2]]._s + (gzs_pml_II + tid - s[D_EV[2]]._dims[0].gs);           /* 1D */
 
    /* adjoint formula not derived */
    /* _rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0; */
    for ( iz = gzs_pml_II+tid; iz < gze_pml_II+1; iz += tsz) { 
      etazdt = (*_epz) * dt2;
      for ( iy = 0; iy < ny; iy ++) {
        sxz2 = _sxz3[-1]; sxz1 = _sxz3[-2]; sxz0 = _sxz3[-3];
        etaydt = (*_epy++) * dt2;
        /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
        /*
          if(!_fwd) { 
          tmp = _px;
          _px=_mpx; 
          _mpx=tmp;
          }
        */
        for ( _vzend = _vz + nx; _vz < _vzend; ) {
          sxz3 = *_sxz3++;
          etaxdt = (*_epx++) * dt2;
        
          dfdx = (sxz3 - sxz0 + (sxz1 - sxz2) * 27.0) * lax; 
          dfdy = ((*_syz3++) - (*_syz0++) + ((*_syz1++) - (*_syz2++)) * 27.0) * lay;
          dfdz = ((*_pz3++) - (*_pz0++) + ((*_pz1++) - (*_pz2++)) * 27.0) * laz;
          if(_fwd) {
            (*_vz_x_II) = ((*_vz_x_II) * (1.0f - etaxdt) + dfdx*(*_mvz))/(1.0f + etaxdt);
            (*_vz_y_II) = ((*_vz_y_II) * (1.0f - etaydt) + dfdy*(*_mvz))/(1.0f + etaydt);
            (*_vz_z_II) = ((*_vz_z_II) * (1.0f - etazdt) + dfdz*(*_mvz))/(1.0f + etazdt);
            *_vz        = *_vz_x_II + *_vz_y_II + * _vz_z_II;
          }
          else {
            /*
              (*_px) = (*_px) + delta/(*_rmpx++);
            */
            /* (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt); */
          }
          _vz++; _mvz++;
          _vz_x_II++; _vz_y_II++; _vz_z_II++; 
                 
          sxz0 = sxz1; sxz1 = sxz2; sxz2 = sxz3;
        }
        /* swap pointers back */
        /*
          if(!_fwd){ 
          tmp = _mpx;
          _mpx=_px; 
          _px=tmp;
          }
        */
        _vz += vz_a;
        _mvz += mvz_a;
        _sxz3 += sxz_a; 
        _syz0 += syz_a; _syz1 += syz_a; _syz2 += syz_a; _syz3 += syz_a; 
        _pz0 += pz_a; _pz1 += pz_a; _pz2 += pz_a; _pz3 += pz_a;
        _epx -= nx;
        _vz_x_II += vz_pml_II_a;
        _vz_y_II += vz_pml_II_a;
        _vz_z_II += vz_pml_II_a;
      }
      _vz += vz_aa;
      _mvz += mvz_aa; 
      _sxz3 += sxz_aa; 
      _syz0 += syz_aa; _syz1 += syz_aa; _syz2 += syz_aa; _syz3 += syz_aa;
      _pz0 += pz_aa; _pz1 += pz_aa; _pz2 += pz_aa; _pz3 += pz_aa;
      _epy -= ny; _epz += tsz;
      _vz_x_II += vz_pml_II_aa;
      _vz_y_II += vz_pml_II_aa;
      _vz_z_II += vz_pml_II_aa;
    }
    /*************************************************************************************/
  }/* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

/*---- END POINTER BRANCH ----------------------------------------------------*/
#endif


/*********************************** original sgn functions ****************************************/
