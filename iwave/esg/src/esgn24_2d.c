/**
   time step functions
   XW 09/24/2011  
   2-4 stagger finite difference in 2D for isotropic elastic wave
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
 * Normal stress updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts2d_24ns (RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Shear stress updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts2d_24ss0(RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Velocity component vz updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts2d_24v0 (RDOM *, RDOM *, RDOM *, void *, int);
/**
 * Velocity component vx updater of 2D stagger-grid 2-4 scheme for EWE  
 */
int esgn_gts2d_24v1 (RDOM *, RDOM *, RDOM *, void *, int);

int esgn_gts2d_24(RDOM *dom, RDOM *rdom, RDOM *cdom, int iarr, void *pars,int _fwd) {

  ESGN_TS_PARS * esgnpars = (ESGN_TS_PARS *) pars;  
  /* sanity test - storage mode */
#ifdef STORE
  if ((dom->_s)[D_P0]._dims[0].n > NSTORE) return E_NOTIMESTEP;
#endif
  if ( iarr == D_P1 || iarr == D_P2)  return 0;
  if ( iarr == D_P0 ) return esgn_gts2d_24ns (dom, rdom, cdom, esgnpars, _fwd);
  if ( iarr == D_S0 ) return esgn_gts2d_24ss0(dom, rdom ,cdom, esgnpars, _fwd);
  if ( iarr == D_V0 ) return esgn_gts2d_24v0 (dom, rdom, cdom, esgnpars, _fwd);
  if ( iarr == D_V1 ) return esgn_gts2d_24v1( dom, rdom, cdom, esgnpars, _fwd);
  
  return E_NOTIMESTEP;
}

/*----------------------------------------------------------------------------*/
#ifdef VEC
int esgn_gts2d_24ns (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  
  fprintf(stderr, "Warning: VEC version 2D normal stress updater not implemented\n");
  return 0;
}

int esgn_gts2d_24ss0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {

  fprintf(stderr, "Warning: VEC version 2D shear stress updater not implemented\n");
  return 0;
}

int esgn_gts2d_24v0 (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  
  fprintf(stderr, "Warning: VEC version 2D v0 updater not implemented\n");
  return 0;
}

int esgn_gts2d_24v1 (RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  fprintf(stderr, "Warning: VEC version 2D v1 updater not implemented\n");
  return 0;
}

#else
/*----------------------------------------------------------------------------*/

int esgn_gts2d_24ns(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, gxe, gye, px_a, py_a, mp00_a, mp01_a, vx_a, vy_a, iy, tsz, tid;
  int px_pml_I_a, py_pml_I_a, px_pml_II_a, py_pml_II_a, px_pml_III_a, py_pml_III_a, px_pml_IV_a, py_pml_IV_a;
  int gxs_pml_III, gxe_pml_III, gxs_pml_IV, gxe_pml_IV, gys_pml_I, gye_pml_I, gys_pml_II, gye_pml_II;
  register ireal * restrict _px, * restrict _py, * restrict _pxend;
  register ireal * restrict _px_x_I, * restrict _px_y_I, * restrict _py_x_I, * restrict _py_y_I;
  register ireal * restrict _px_x_II, * restrict _px_y_II, * restrict _py_x_II, * restrict _py_y_II;
  register ireal * restrict _px_x_III, * restrict _px_y_III, * restrict _py_x_III, * restrict _py_y_III;
  register ireal * restrict _px_x_IV, * restrict _px_y_IV, * restrict _py_x_IV, * restrict _py_y_IV;
  register ireal * restrict _mp00, * restrict _mp01;
  register ireal * restrict _vx3, * restrict _vy3, * restrict _vy2;
  register ireal * restrict _vy1, * restrict _vy0;
  register ireal * restrict _epx, * restrict _epy;
  register ireal lax, lay, dt2, vx3, vx2, vx1, vx0, dfdx, dfdy, etaxdt, etaydt;

  //  register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  // register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_P0]._dims[0].gs;
  gxe = s[D_P0]._dims[0].ge;
  gys = s[D_P0]._dims[1].gs;
  gye = s[D_P0]._dims[1].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_P0,&empty);
  if (empty) {
    gys_pml_I = gys;
    gye_pml_I = gys-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gys_pml_I = s_pml[D_P0]._dims[1].gs;
    gye_pml_I = s_pml[D_P0]._dims[1].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+2,D_P0,&empty);
  if (empty) {
    gys_pml_II = gye+1;
    gye_pml_II = gye;
  }
  else {
    s_pml = ld_pml[2]._s;
    gys_pml_II = s_pml[D_P0]._dims[1].gs;
    gye_pml_II = s_pml[D_P0]._dims[1].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+4,D_P0,&empty);
  if (empty) {
    gxs_pml_III = gxs;
    gxe_pml_III = gxs-1;
  }
  else {
    s_pml = ld_pml[4]._s;
    gxs_pml_III = s_pml[D_P0]._dims[0].gs;
    gxe_pml_III = s_pml[D_P0]._dims[0].ge;
  }
 /** pml region IV */
  rd_empty(ld_pml+6,D_P0,&empty);
  if (empty) {
    gxs_pml_IV = gxe+1;
    gxe_pml_IV = gxe;
  }
  else {
    s_pml = ld_pml[6]._s;
    gxs_pml_IV = s_pml[D_P0]._dims[0].gs;
    gxe_pml_IV = s_pml[D_P0]._dims[0].ge;
  }
  
#pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mp00,_mp01,_epx,_epy,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,_px_x_I,_px_y_I,_py_x_I,_py_y_I,_px_x_II,_px_y_II,_py_x_II,_py_y_II,_px_x_III,_px_y_III,_py_x_III,_py_y_III,_px_x_IV,_px_y_IV,_py_x_IV,_py_y_IV,dfdx,dfdy,etaxdt,etaydt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

    /** pml region I *********************************************************************/
#pragma omp single
    {            
      px_a         = tsz * s[D_P0]._dims[0].n0 - nx;
      py_a         = tsz * s[D_P1]._dims[0].n0 - nx;
      
      px_pml_I_a   = tsz * ld_pml[0]._s[D_P0]._dims[0].n0-nx;
      py_pml_I_a   = tsz * ld_pml[0]._s[D_P1]._dims[0].n0-nx;

      px_pml_II_a  = tsz * ld_pml[2]._s[D_P0]._dims[0].n0-nx;
      py_pml_II_a  = tsz * ld_pml[2]._s[D_P1]._dims[0].n0-nx;

      px_pml_III_a = tsz * ld_pml[4]._s[D_P0]._dims[0].n0-(gxe_pml_III-gxs_pml_III+1);
      py_pml_III_a = tsz * ld_pml[4]._s[D_P1]._dims[0].n0-(gxe_pml_III-gxs_pml_III+1);

      px_pml_IV_a  = tsz * ld_pml[6]._s[D_P0]._dims[0].n0-(gxe_pml_IV-gxs_pml_IV+1);
      py_pml_IV_a  = tsz * ld_pml[6]._s[D_P1]._dims[0].n0-(gxe_pml_IV-gxs_pml_IV+1);

      mp00_a       = tsz * s[D_MP00]._dims[0].n0 - nx;
      mp01_a       = tsz * s[D_MP01]._dims[0].n0 - nx;

      vx_a         = tsz * s[D_V0]._dims[0].n0 - nx;
      vy_a         = tsz * s[D_V1]._dims[0].n0 - nx;
    }
#pragma omp barrier
    /** gys == gys_pml_I */
    _px     = s[D_P0]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py     = s[D_P1]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    _px_x_I = s_pml[D_P0]._s + (gxs - s_pml[D_P0]._dims[0].gs) + (gys - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_x_I = s_pml[D_P1]._s + (gxs - s_pml[D_P1]._dims[0].gs) + (gys - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
    
    s_pml   = ld_pml[1]._s;
    _px_y_I = s_pml[D_P0]._s + (gxs - s_pml[D_P0]._dims[0].gs) + (gys - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_y_I = s_pml[D_P1]._s + (gxs - s_pml[D_P1]._dims[0].gs) + (gys - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
   
    _mp00   = cs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    _mp01   = cs[D_MP01]._s + (gxs - s[D_MP01]._dims[0].gs) + (gys - s[D_MP01]._dims[1].gs + tid) * s[D_MP01]._dims[0].n0;
    
    _vx3    = rs[D_V0]._s + (gxs - s[D_V0]._dims[0].gs) + (gys - s[D_V0]._dims[1].gs + tid) * s[D_V0]._dims[0].n0 + 1;
    _vy2    = rs[D_V1]._s + (gxs - s[D_V1]._dims[0].gs) + (gys - s[D_V1]._dims[1].gs + tid) * s[D_V1]._dims[0].n0;
    _vy3    = _vy2 + s[D_V1]._dims[0].n0;
    _vy1    = _vy2 - s[D_V1]._dims[0].n0;
    _vy0    = _vy1 - s[D_V1]._dims[0].n0;
    _epx    = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);              /* 1D */
    _epy    = s[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs + tid);        /* 1D */
    
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
   
    for ( iy = gys_pml_I+tid; iy < gye_pml_I+1; iy += tsz )
    {
      vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _pxend = _px + nx; _px < _pxend; )
      {
        vx3 = *_vx3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
        dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_px_x_I) = ((*_px_x_I) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
          (*_py_x_I) = ((*_py_x_I) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
          (*_px_y_I) = ((*_px_y_I) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
          (*_py_y_I) = ((*_py_y_I) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);

          *_px = *_px_x_I + *_px_y_I;
          *_py = *_py_x_I + *_py_y_I;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _px++; _py++; _mp00++; _mp01++;
        _px_x_I++; _px_y_I++; _py_x_I++; _py_y_I++;
	
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
      _px += px_a; _py += py_a; _mp00 += mp00_a; _mp01 += mp01_a; _vx3 += vx_a; _vy0 += vy_a;
      _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx; _epy += tsz;
      _px_x_I += px_pml_I_a; _py_x_I += py_pml_I_a;
      _px_y_I += px_pml_I_a; _py_y_I += py_pml_I_a;
    }
    /*************************************************************************************/

    /** pml region III, IV and physical region *******************************************/
    /** adjust pointers */
    _px   += (gye_pml_I + 1 + tid - iy) * s[D_P0 ]._dims[0].n0;
    _py   += (gye_pml_I + 1 + tid - iy) * s[D_P1 ]._dims[0].n0;
    _mp00 += (gye_pml_I + 1 + tid - iy) * s[D_MP00]._dims[0].n0;
    _mp01 += (gye_pml_I + 1 + tid - iy) * s[D_MP01]._dims[0].n0;
    _vx3  += (gye_pml_I + 1 + tid - iy) * s[D_V0 ]._dims[0].n0;
    _vy2  += (gye_pml_I + 1 + tid - iy) * s[D_V1 ]._dims[0].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
  
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gye_pml_I + 1 + tid - s[D_EP[1]]._dims[0].gs);        /* 1D */
    
    s_pml = ld_pml[4]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    _px_x_III = s_pml[D_P0]._s + (gxs_pml_III - s_pml[D_P0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_x_III = s_pml[D_P1]._s + (gxs_pml_III - s_pml[D_P1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
   
    s_pml = ld_pml[5]._s;
    _px_y_III = s_pml[D_P0]._s + (gxs_pml_III - s_pml[D_P0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_y_III = s_pml[D_P1]._s + (gxs_pml_III - s_pml[D_P1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
    
    s_pml = ld_pml[6]._s;
    _px_x_IV = s_pml[D_P0]._s + (gxs_pml_IV - s_pml[D_P0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_x_IV = s_pml[D_P1]._s + (gxs_pml_IV - s_pml[D_P1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
    
    s_pml = ld_pml[7]._s;
    _px_y_IV = s_pml[D_P0]._s + (gxs_pml_IV - s_pml[D_P0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_y_IV = s_pml[D_P1]._s + (gxs_pml_IV - s_pml[D_P1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
    
    /** adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gye_pml_I+1+tid; iy < gys_pml_II; iy += tsz )
    {
      vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
      //etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      /** pml region III */
      for ( _pxend = _px + gxe_pml_III-gxs_pml_III+1; _px < _pxend; )
      {
        vx3 = *_vx3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
        dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_px_x_III) = ((*_px_x_III) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
          (*_py_x_III) = ((*_py_x_III) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
          (*_px_y_III) = (*_px_y_III) + dfdy*(*_mp01);
          (*_py_y_III) = (*_py_y_III) + dfdy*(*_mp00);

          *_px = *_px_x_III + *_px_y_III;
          *_py = *_py_x_III + *_py_y_III;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _px++; _py++; _mp00++; _mp01++;
        _px_x_III++; _px_y_III++; _py_x_III++; _py_y_III++;
	
        vx0 = vx1; vx1 = vx2; vx2 = vx3;
      }
      
      /** physical region */
      for ( _pxend = _px + gxs_pml_IV-gxe_pml_III-1; _px < _pxend; )
      {
        vx3 = *_vx3++;
        _epx++;
        
        dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
        dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
        if(_fwd) {
          *_px = *_px + dfdx*(*_mp00) + dfdy*(*_mp01);
          *_py = *_py + dfdx*(*_mp01) + dfdy*(*_mp00);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _px++; _py++; _mp00++; _mp01++;	

        vx0 = vx1; vx1 = vx2; vx2 = vx3;
      }
      
      /** pml region IV */
      for ( _pxend = _px + gxe_pml_IV-gxs_pml_IV+1; _px < _pxend; )
      {
        vx3 = *_vx3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
        dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_px_x_IV) = ((*_px_x_IV) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
          (*_py_x_IV) = ((*_py_x_IV) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
          (*_px_y_IV) = (*_px_y_IV) + dfdy*(*_mp01);
          (*_py_y_IV) = (*_py_y_IV) + dfdy*(*_mp00);

          *_px = *_px_x_IV + *_px_y_IV;
          *_py = *_py_x_IV + *_py_y_IV;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _px++; _py++; _mp00++; _mp01++;	
        _px_x_IV++; _px_y_IV++; _py_x_IV++; _py_y_IV++;
        
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
      _px += px_a; _py += py_a; _mp00 += mp00_a; _mp01 += mp01_a; _vx3 += vx_a; _vy0 += vy_a;
      _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx;
      _px_x_III += px_pml_III_a; _py_x_III += py_pml_III_a;
      _px_y_III += px_pml_III_a; _py_y_III += py_pml_III_a;
      _px_x_IV += px_pml_IV_a; _py_x_IV += py_pml_IV_a;
      _px_y_IV += px_pml_IV_a; _py_y_IV += py_pml_IV_a;
    
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _px   += (gys_pml_II + tid - iy) * s[D_P0 ]._dims[0].n0;
    _py   += (gys_pml_II + tid - iy) * s[D_P1 ]._dims[0].n0;
    _mp00 += (gys_pml_II + tid - iy) * s[D_MP00]._dims[0].n0;
    _mp01 += (gys_pml_II + tid - iy) * s[D_MP01]._dims[0].n0;
    _vx3  += (gys_pml_II + tid - iy) * s[D_V0 ]._dims[0].n0;
    _vy2  += (gys_pml_II + tid - iy) * s[D_V1 ]._dims[0].n0;
    
    _vy3 = _vy2 + s[D_V1]._dims[0].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
  
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys_pml_II - s[D_EP[1]]._dims[0].gs + tid);           /* 1D */
    
    s_pml = ld_pml[2]._s;
    _px_x_II = s_pml[D_P0]._s + (gxs - s_pml[D_P0]._dims[0].gs) + (gys_pml_II - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_x_II = s_pml[D_P1]._s + (gxs - s_pml[D_P1]._dims[0].gs) + (gys_pml_II - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
    
    s_pml = ld_pml[3]._s;
    _px_y_II = s_pml[D_P0]._s + (gxs - s_pml[D_P0]._dims[0].gs) + (gys_pml_II - s_pml[D_P0]._dims[1].gs + tid) * s_pml[D_P0]._dims[0].n0;
    _py_y_II = s_pml[D_P1]._s + (gxs - s_pml[D_P1]._dims[0].gs) + (gys_pml_II - s_pml[D_P1]._dims[1].gs + tid) * s_pml[D_P1]._dims[0].n0;
 
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gys_pml_II+tid; iy < gye_pml_II+1; iy += tsz )
    {
      vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _pxend = _px + nx; _px < _pxend; )
      {
        vx3 = *_vx3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vx3 - vx0 + (vx1 - vx2) * 27.0) * lax; 
        dfdy = ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_px_x_II) = ((*_px_x_II) * (1.0f - etaxdt) + dfdx*(*_mp00))/(1.0f + etaxdt);
          (*_py_x_II) = ((*_py_x_II) * (1.0f - etaxdt) + dfdx*(*_mp01))/(1.0f + etaxdt);
          (*_px_y_II) = ((*_px_y_II) * (1.0f - etaydt) + dfdy*(*_mp01))/(1.0f + etaydt);
          (*_py_y_II) = ((*_py_y_II) * (1.0f - etaydt) + dfdy*(*_mp00))/(1.0f + etaydt);

          *_px = *_px_x_II + *_px_y_II;
          *_py = *_py_x_II + *_py_y_II;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _px++; _py++; _mp00++; _mp01++;
        _px_x_II++; _px_y_II++; _py_x_II++; _py_y_II++;
	
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
      _px += px_a; _py += py_a; _mp00 += mp00_a; _mp01 += mp01_a; _vx3 += vx_a; _vy0 += vy_a;
      _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx; _epy += tsz;
      _px_x_II += px_pml_II_a; _py_x_II += py_pml_II_a;
      _px_y_II += px_pml_II_a; _py_y_II += py_pml_II_a;
    }
    /*************************************************************************************/

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/
int esgn_gts2d_24ss0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, gxe, gye, iy, tsz, tid, sxy_a, ms0_a, vx_a, vy_a;
  int sxy_pml_I_a, sxy_pml_II_a, sxy_pml_III_a, sxy_pml_IV_a;
  int gxs_pml_III, gxe_pml_III, gxs_pml_IV, gxe_pml_IV, gys_pml_I, gye_pml_I, gys_pml_II, gye_pml_II;
  register ireal * restrict _sxy, * restrict _sxyend;
  register ireal * restrict _sxy_x_I, * restrict _sxy_y_I;
  register ireal * restrict _sxy_x_II, * restrict _sxy_y_II;
  register ireal * restrict _sxy_x_III, * restrict _sxy_y_III;
  register ireal * restrict _sxy_x_IV, * restrict _sxy_y_IV;
  register ireal * restrict _ms0;
  register ireal * restrict _vx3, * restrict _vy3, * restrict _vx2;
  register ireal * restrict _vx1, * restrict _vx0;
  register ireal * restrict _epx, * restrict _epy;
  register ireal lax, lay, dt2, vy3, vy2, vy1, vy0, dfdx, dfdy, etaxdt, etaydt;

  //register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  //register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_S0]._dims[0].n;
  ny = s[D_S0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_S0]._dims[0].gs;
  gxe = s[D_S0]._dims[0].ge;
  gys = s[D_S0]._dims[1].gs;
  gye = s[D_S0]._dims[1].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_S0,&empty);
  if (empty) {
    gys_pml_I = gys;
    gye_pml_I = gys-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gys_pml_I = s_pml[D_S0]._dims[1].gs;
    gye_pml_I = s_pml[D_S0]._dims[1].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+2,D_S0,&empty);
  if (empty) {
    gys_pml_II = gye+1;
    gye_pml_II = gye;
  }
  else {
    s_pml = ld_pml[2]._s;
    gys_pml_II = s_pml[D_S0]._dims[1].gs;
    gye_pml_II = s_pml[D_S0]._dims[1].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+4,D_S0,&empty);
  if (empty) {
    gxs_pml_III = gxs;
    gxe_pml_III = gxs-1;
  }
  else {
    s_pml = ld_pml[4]._s;
    gxs_pml_III = s_pml[D_S0]._dims[0].gs;
    gxe_pml_III = s_pml[D_S0]._dims[0].ge;
  }
 /** pml region IV */
  rd_empty(ld_pml+6,D_S0,&empty);
  if (empty) {
    gxs_pml_IV = gxe+1;
    gxe_pml_IV = gxe;
  }
  else {
    s_pml = ld_pml[6]._s;
    gxs_pml_IV = s_pml[D_S0]._dims[0].gs;
    gxe_pml_IV = s_pml[D_S0]._dims[0].ge;
  }
  
#pragma omp parallel private(tsz,tid,iy,_sxy,_sxyend,_ms0,_epx,_epy,_vy3,_vx3,_vx2,_vx1,_vx0,vy3,vy2,vy1,vy0,_sxy_x_I,_sxy_y_I,_sxy_x_II,_sxy_y_II,_sxy_x_III,_sxy_y_III,_sxy_x_IV,_sxy_y_IV,dfdx,dfdy,etaxdt,etaydt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

    /** pml region I *********************************************************************/
#pragma omp single
    {            
      sxy_a         = tsz * s[D_S0]._dims[0].n0 - nx;
      ms0_a         = tsz * s[D_MS0]._dims[0].n0 - nx;
      vx_a          = tsz * s[D_V0]._dims[0].n0 - nx;
      vy_a          = tsz * s[D_V1]._dims[0].n0 - nx;
   
      sxy_pml_I_a   = tsz * ld_pml[0]._s[D_S0]._dims[0].n0-nx;
      sxy_pml_II_a  = tsz * ld_pml[2]._s[D_S0]._dims[0].n0-nx;   
      sxy_pml_III_a = tsz * ld_pml[4]._s[D_S0]._dims[0].n0-(gxe_pml_III-gxs_pml_III+1);  
      sxy_pml_IV_a  = tsz * ld_pml[6]._s[D_S0]._dims[0].n0-(gxe_pml_IV-gxs_pml_IV+1);
    }
#pragma omp barrier
    /** gys == gys_pml_I */
    _sxy     = s[D_S0]._s + (gxs - s[D_S0 ]._dims[0].gs) + (gys - s[D_S0 ]._dims[1].gs + tid) * s[D_S0 ]._dims[0].n0;
      
    s_pml    = ld_pml[0]._s;
    _sxy_x_I = s_pml[D_S0]._s + (gxs - s_pml[D_S0]._dims[0].gs) + (gys - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
      
    s_pml    = ld_pml[1]._s;
    _sxy_y_I = s_pml[D_S0]._s + (gxs - s_pml[D_S0]._dims[0].gs) + (gys - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
      
    _ms0     = cs[D_MS0]._s + (gxs - s[D_MS0]._dims[0].gs) + (gys - s[D_MS0]._dims[1].gs + tid) * s[D_MS0]._dims[0].n0;
        
    _vy3     = rs[D_V1]._s + (gxs - s[D_V1]._dims[0].gs) + (gys - s[D_V1]._dims[1].gs + tid) * s[D_V1]._dims[0].n0 + 2;
    _vx1     = rs[D_V0]._s + (gxs - s[D_V0]._dims[0].gs) + (gys - s[D_V0]._dims[1].gs + tid) * s[D_V0]._dims[0].n0;
    _vx2     = _vx1 + s[D_V0]._dims[0].n0;
    _vx3     = _vx2 + s[D_V0]._dims[0].n0;
    _vx0     = _vx1 - s[D_V0]._dims[0].n0;

    _epx     = rs[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);              /* 1D */
    _epy     = rs[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs + tid);        /* 1D */
    
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
   
    for ( iy = gys_pml_I+tid; iy < gye_pml_I+1; iy += tsz )
    {
      vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _sxyend = _sxy + nx; _sxy < _sxyend; )
      {
        vy3 = *_vy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
        dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
        if(_fwd) {
          (*_sxy_x_I) = ((*_sxy_x_I) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
          (*_sxy_y_I) = ((*_sxy_y_I) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
        
          *_sxy = *_sxy_x_I + *_sxy_y_I;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
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
      _sxy += sxy_a; _ms0 += ms0_a; _vy3 += vy_a; _vx0 += vx_a;
      _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; _epx -= nx; _epy += tsz;
      _sxy_x_I += sxy_pml_I_a; 
      _sxy_y_I += sxy_pml_I_a;
    }
    /*************************************************************************************/

    /** pml region III, IV and physical region *******************************************/
    /** adjust pointers */
    _sxy  += (gye_pml_I + 1 + tid - iy) * s[D_S0 ]._dims[0].n0;
    _ms0  += (gye_pml_I + 1 + tid - iy) * s[D_MS0]._dims[0].n0;
    _vy3  += (gye_pml_I + 1 + tid - iy) * s[D_V1 ]._dims[0].n0;
    _vx2  += (gye_pml_I + 1 + tid - iy) * s[D_V0 ]._dims[0].n0;
    
    _vx3 = _vx2 + s[D_V0]._dims[0].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0;
  
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gye_pml_I + 1 + tid - s[D_EV[1]]._dims[0].gs);        /* 1D */
    
    s_pml = ld_pml[4]._s;
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    _sxy_x_III = s_pml[D_S0]._s + (gxs_pml_III - s_pml[D_S0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
      
    s_pml = ld_pml[5]._s;
    _sxy_y_III = s_pml[D_S0]._s + (gxs_pml_III - s_pml[D_S0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
       
    s_pml = ld_pml[6]._s;
    _sxy_x_IV = s_pml[D_S0]._s + (gxs_pml_IV - s_pml[D_S0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
      
    s_pml = ld_pml[7]._s;
    _sxy_y_IV = s_pml[D_S0]._s + (gxs_pml_IV - s_pml[D_S0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
      
    /** adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gye_pml_I+1+tid; iy < gys_pml_II; iy += tsz )
    {
      vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
      //etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      /** pml region III */
      for ( _sxyend = _sxy + gxe_pml_III-gxs_pml_III+1; _sxy < _sxyend; )
      {
        vy3 = *_vy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
        dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
        if(_fwd) {
          (*_sxy_x_III) = ((*_sxy_x_III) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
          (*_sxy_y_III) = (*_sxy_y_III) + dfdy*(*_ms0);
          
          *_sxy = *_sxy_x_III + *_sxy_y_III;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _sxy++; _ms0++;
        _sxy_x_III++; _sxy_y_III++;
	
        vy0 = vy1; vy1 = vy2; vy2 = vy3;
      }
      
      /** physical region */
      for ( _sxyend = _sxy + gxs_pml_IV-gxe_pml_III-1; _sxy < _sxyend; )
      {
        vy3 = *_vy3++;
        _epx++;
        
        dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
        dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
        if(_fwd) {
          *_sxy = *_sxy + (dfdx + dfdy) * (*_ms0);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _sxy++; _ms0++;

        vy0 = vy1; vy1 = vy2; vy2 = vy3;
      }
      
      /** pml region IV */
      for ( _sxyend = _sxy + gxe_pml_IV-gxs_pml_IV+1; _sxy < _sxyend; )
      {
        vy3 = *_vy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
        dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
        if(_fwd) {
          (*_sxy_x_IV) = ((*_sxy_x_IV) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
          (*_sxy_y_IV) = (*_sxy_y_IV) + dfdy*(*_ms0);
          
          *_sxy = *_sxy_x_IV + *_sxy_y_IV;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _sxy++; _ms0++;
        _sxy_x_IV++; _sxy_y_IV++; 
        
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
      _sxy += sxy_a; _ms0 += ms0_a; _vy3 += vy_a; _vx0 += vx_a;
      _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; _epx -= nx;
      _sxy_x_III += sxy_pml_III_a; 
      _sxy_y_III += sxy_pml_III_a;
      _sxy_x_IV  += sxy_pml_IV_a; 
      _sxy_y_IV  += sxy_pml_IV_a; 
    
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _sxy  += (gys_pml_II + tid - iy) * s[D_S0 ]._dims[0].n0;
    _ms0  += (gys_pml_II + tid - iy) * s[D_MS0]._dims[0].n0;
    _vy3  += (gys_pml_II + tid - iy) * s[D_V1 ]._dims[0].n0;
    _vx2  += (gys_pml_II + tid - iy) * s[D_V0 ]._dims[0].n0;
    
    _vx3 = _vx2 + s[D_V0]._dims[0].n0;
    _vx1 = _vx2 - s[D_V0]._dims[0].n0;
    _vx0 = _vx1 - s[D_V0]._dims[0].n0;
  
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys_pml_II - s[D_EV[1]]._dims[0].gs + tid);           /* 1D */
    
    s_pml = ld_pml[2]._s;
    _sxy_x_II = s_pml[D_S0]._s + (gxs - s_pml[D_S0]._dims[0].gs) + (gys_pml_II - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
        
    s_pml = ld_pml[3]._s;
    _sxy_y_II = s_pml[D_S0]._s + (gxs - s_pml[D_S0]._dims[0].gs) + (gys_pml_II - s_pml[D_S0]._dims[1].gs + tid) * s_pml[D_S0]._dims[0].n0;
    
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gys_pml_II+tid; iy < gye_pml_II+1; iy += tsz )
    {
      vy2 = _vy3[-1]; vy1 = _vy3[-2]; vy0 = _vy3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _sxyend = _sxy + nx; _sxy < _sxyend; )
      {
        vy3 = *_vy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (vy3 - vy0 + (vy1 - vy2) * 27.0) * lax; 
        dfdy = ((*_vx3++) - (*_vx0++) + ((*_vx1++) - (*_vx2++)) * 27.0) * lay;
        if(_fwd) {
          (*_sxy_x_II) = ((*_sxy_x_II) * (1.0f - etaxdt) + dfdx*(*_ms0))/(1.0f + etaxdt);
          (*_sxy_y_II) = ((*_sxy_y_II) * (1.0f - etaydt) + dfdy*(*_ms0))/(1.0f + etaydt);
          
          *_sxy = *_sxy_x_II + *_sxy_y_II;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
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
      _sxy += sxy_a; _ms0 += ms0_a; _vy3 += vy_a; _vx0 += vx_a;
      _vx1 += vx_a; _vx2 += vx_a; _vx3 += vx_a; _epx -= nx; _epy += tsz;
      _sxy_x_II += sxy_pml_II_a;
      _sxy_y_II += sxy_pml_II_a;
    }
    /*************************************************************************************/

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/
int esgn_gts2d_24v0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, gxe, gye, iy, tsz, tid, vx_a, mvx_a, px_a, sxy_a;
  int vx_pml_I_a, vx_pml_II_a, vx_pml_III_a, vx_pml_IV_a;
  int gxs_pml_III, gxe_pml_III, gxs_pml_IV, gxe_pml_IV, gys_pml_I, gye_pml_I, gys_pml_II, gye_pml_II;
  register ireal * restrict _vx, * restrict _mvx, * restrict _vxend;
  register ireal * restrict _vx_x_I, * restrict _vx_y_I;
  register ireal * restrict _vx_x_II, * restrict _vx_y_II;
  register ireal * restrict _vx_x_III, * restrict _vx_y_III;
  register ireal * restrict _vx_x_IV, * restrict _vx_y_IV;
  register ireal * restrict _px3, * restrict _sxy3, * restrict _sxy2;
  register ireal * restrict _sxy1, * restrict _sxy0;
  register ireal * restrict _epx, * restrict _epy;
  register ireal lax, lay, dt2, px3, px2, px1, px0, dfdx, dfdy, etaxdt, etaydt;

  //register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  //register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_V0]._dims[0].gs;
  gxe = s[D_V0]._dims[0].ge;
  gys = s[D_V0]._dims[1].gs;
  gye = s[D_V0]._dims[1].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_V0,&empty);
  if (empty) {
    gys_pml_I = gys;
    gye_pml_I = gys-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gys_pml_I = s_pml[D_V0]._dims[1].gs;
    gye_pml_I = s_pml[D_V0]._dims[1].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+2,D_V0,&empty);
  if (empty) {
    gys_pml_II = gye+1;
    gye_pml_II = gye;
  }
  else {
    s_pml = ld_pml[2]._s;
    gys_pml_II = s_pml[D_V0]._dims[1].gs;
    gye_pml_II = s_pml[D_V0]._dims[1].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+4,D_V0,&empty);
  if (empty) {
    gxs_pml_III = gxs;
    gxe_pml_III = gxs-1;
  }
  else {
    s_pml = ld_pml[4]._s;
    gxs_pml_III = s_pml[D_V0]._dims[0].gs;
    gxe_pml_III = s_pml[D_V0]._dims[0].ge;
  }
 /** pml region IV */
  rd_empty(ld_pml+6,D_V0,&empty);
  if (empty) {
    gxs_pml_IV = gxe+1;
    gxe_pml_IV = gxe;
  }
  else {
    s_pml = ld_pml[6]._s;
    gxs_pml_IV = s_pml[D_V0]._dims[0].gs;
    gxe_pml_IV = s_pml[D_V0]._dims[0].ge;
  }
  
#pragma omp parallel private(tsz,tid,iy,_vx,_vxend,_mvx,_epx,_epy,_px3,_sxy3,_sxy2,_sxy1,_sxy0,px3,px2,px1,px0,_vx_x_I,_vx_y_I,_vx_x_II,_vx_y_II,_vx_x_III,_vx_y_III,_vx_x_IV,_vx_y_IV,dfdx,dfdy,etaxdt,etaydt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

    /** pml region I *********************************************************************/
#pragma omp single
    {            
      vx_a         = tsz * s[D_V0]._dims[0].n0 - nx;
        
      vx_pml_I_a   = tsz * ld_pml[0]._s[D_V0]._dims[0].n0-nx;
      vx_pml_II_a  = tsz * ld_pml[2]._s[D_V0]._dims[0].n0-nx; 
      vx_pml_III_a = tsz * ld_pml[4]._s[D_V0]._dims[0].n0-(gxe_pml_III-gxs_pml_III+1);
      vx_pml_IV_a  = tsz * ld_pml[6]._s[D_V0]._dims[0].n0-(gxe_pml_IV-gxs_pml_IV+1);
      
      mvx_a        = tsz * s[D_MV0]._dims[0].n0 - nx;
 
      px_a         = tsz * s[D_P0]._dims[0].n0 - nx;
      sxy_a        = tsz * s[D_S0]._dims[0].n0 - nx;
    }
#pragma omp barrier
    /** gys == gys_pml_I */
    _vx     = s[D_V0]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
       
    s_pml   = ld_pml[0]._s;
    _vx_x_I = s_pml[D_V0]._s + (gxs - s_pml[D_V0]._dims[0].gs) + (gys - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
       
    s_pml   = ld_pml[1]._s;
    _vx_y_I = s_pml[D_V0]._s + (gxs - s_pml[D_V0]._dims[0].gs) + (gys - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
      
    _mvx    = cs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs + tid) * s[D_MV0]._dims[0].n0;
       
    _px3    = rs[D_P0]._s + (gxs - s[D_P0]._dims[0].gs) + (gys - s[D_P0]._dims[1].gs + tid) * s[D_P0]._dims[0].n0 + 2;
    _sxy2   = rs[D_S0]._s + (gxs - s[D_S0]._dims[0].gs) + (gys - s[D_S0]._dims[1].gs + tid) * s[D_S0]._dims[0].n0;
    _sxy3   = _sxy2 + s[D_S0]._dims[0].n0;
    _sxy1   = _sxy2 - s[D_S0]._dims[0].n0;
    _sxy0   = _sxy1 - s[D_S0]._dims[0].n0;
    _epx    = rs[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);              /* 1D */
    _epy    = rs[D_EP[1]]._s + (gys - s[D_EP[1]]._dims[0].gs + tid);        /* 1D */
    
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
   
    for ( iy = gys_pml_I+tid; iy < gye_pml_I+1; iy += tsz )
    {
      px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _vxend = _vx + nx; _vx < _vxend; )
      {
        px3 = *_px3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
        dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vx_x_I) = ((*_vx_x_I) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
          (*_vx_y_I) = ((*_vx_y_I) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
          
          *_vx = *_vx_x_I + *_vx_y_I;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vx++; _mvx++;
        _vx_x_I++; _vx_y_I++;
	
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
      _vx += vx_a; _mvx += mvx_a; _px3 += px_a; _sxy0 += sxy_a;
      _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; _epx -= nx; _epy += tsz;
      _vx_x_I += vx_pml_I_a;
      _vx_y_I += vx_pml_I_a;
    }
    /*************************************************************************************/

    /** pml region III, IV and physical region *******************************************/
    /** adjust pointers */
    _vx   += (gye_pml_I + 1 + tid - iy) * s[D_V0 ]._dims[0].n0;
    _mvx  += (gye_pml_I + 1 + tid - iy) * s[D_MV0]._dims[0].n0;
   
    _px3  += (gye_pml_I + 1 + tid - iy) * s[D_P0 ]._dims[0].n0;
    _sxy2 += (gye_pml_I + 1 + tid - iy) * s[D_S0 ]._dims[0].n0;
    
    _sxy3 = _sxy2 + s[D_S0]._dims[0].n0;
    _sxy1 = _sxy2 - s[D_S0]._dims[0].n0;
    _sxy0 = _sxy1 - s[D_S0]._dims[0].n0;
 
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gye_pml_I + 1 + tid - s[D_EP[1]]._dims[0].gs);        /* 1D */
    
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    s_pml = ld_pml[4]._s;
    _vx_x_III = s_pml[D_V0]._s + (gxs_pml_III - s_pml[D_V0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
      
    s_pml = ld_pml[5]._s;
    _vx_y_III = s_pml[D_V0]._s + (gxs_pml_III - s_pml[D_V0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
       
    s_pml = ld_pml[6]._s;
    _vx_x_IV = s_pml[D_V0]._s + (gxs_pml_IV - s_pml[D_V0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
   
    s_pml = ld_pml[7]._s;
    _vx_y_IV = s_pml[D_V0]._s + (gxs_pml_IV - s_pml[D_V0]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
        
    /** adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gye_pml_I+1+tid; iy < gys_pml_II; iy += tsz )
    {
      px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
      //etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      /** pml region III */
      for ( _vxend = _vx + gxe_pml_III-gxs_pml_III+1; _vx < _vxend; )
      {
        px3 = *_px3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
        dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vx_x_III) = ((*_vx_x_III) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
          (*_vx_y_III) = (*_vx_y_III) + dfdy*(*_mvx);
          
          *_vx = *_vx_x_III + *_vx_y_III;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vx++; _mvx++;
        _vx_x_III++; _vx_y_III++;
	
        px0 = px1; px1 = px2; px2 = px3;
      }
      
      /** physical region */
      for ( _vxend = _vx + gxs_pml_IV-gxe_pml_III-1; _vx < _vxend; )
      {
        px3 = *_px3++;
        _epx++;
        
        dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
        dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
        if(_fwd) {
          *_vx = *_vx + (dfdx + dfdy) * (*_mvx);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vx++; _mvx++;	

        px0 = px1; px1 = px2; px2 = px3;
      }
      
      /** pml region IV */
      for ( _vxend = _vx + gxe_pml_IV-gxs_pml_IV+1; _vx < _vxend; )
      {
        px3 = *_px3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
        dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vx_x_IV) = ((*_vx_x_IV) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
          (*_vx_y_IV) = (*_vx_y_IV) + dfdy*(*_mvx);
          
          *_vx = *_vx_x_IV + *_vx_y_IV;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vx++; _mvx++;
        _vx_x_IV++; _vx_y_IV++;
        
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
      _vx += vx_a; _mvx += mvx_a; _px3 += px_a; _sxy0 += sxy_a;
      _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; _epx -= nx;
      _vx_x_III += vx_pml_III_a;
      _vx_y_III += vx_pml_III_a;
      _vx_x_IV += vx_pml_IV_a;
      _vx_y_IV += vx_pml_IV_a;
    
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _vx   += (gys_pml_II + tid - iy) * s[D_V0 ]._dims[0].n0;
    _mvx  += (gys_pml_II + tid - iy) * s[D_MV0]._dims[0].n0;

    _px3  += (gys_pml_II + tid - iy) * s[D_P0 ]._dims[0].n0;
    _sxy2 += (gys_pml_II + tid - iy) * s[D_S0 ]._dims[0].n0;
    
    _sxy3 = _sxy2 + s[D_S0]._dims[0].n0;
    _sxy1 = _sxy2 - s[D_S0]._dims[0].n0;
    _sxy0 = _sxy1 - s[D_S0]._dims[0].n0;
  
    _epx = s[D_EV[0]]._s + (gxs - s[D_EV[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EP[1]]._s + (gys_pml_II - s[D_EP[1]]._dims[0].gs + tid);           /* 1D */
    
    s_pml = ld_pml[2]._s;
    _vx_x_II = s_pml[D_V0]._s + (gxs - s_pml[D_V0]._dims[0].gs) + (gys_pml_II - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
        
    s_pml = ld_pml[3]._s;
    _vx_y_II = s_pml[D_V0]._s + (gxs - s_pml[D_V0]._dims[0].gs) + (gys_pml_II - s_pml[D_V0]._dims[1].gs + tid) * s_pml[D_V0]._dims[0].n0;
     
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gys_pml_II+tid; iy < gye_pml_II+1; iy += tsz )
    {
      px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _vxend = _vx + nx; _vx < _vxend; )
      {
        px3 = *_px3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (px3 - px0 + (px1 - px2) * 27.0) * lax; 
        dfdy = ((*_sxy3++) - (*_sxy0++) + ((*_sxy1++) - (*_sxy2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vx_x_II) = ((*_vx_x_II) * (1.0f - etaxdt) + dfdx*(*_mvx))/(1.0f + etaxdt);
          (*_vx_y_II) = ((*_vx_y_II) * (1.0f - etaydt) + dfdy*(*_mvx))/(1.0f + etaydt);
          
          *_vx = *_vx_x_II + *_vx_y_II;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vx++; _mvx++;
        _vx_x_II++; _vx_y_II++;
	
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
      _vx += vx_a; _mvx += mvx_a; _px3 += px_a; _sxy0 += sxy_a;
      _sxy1 += sxy_a; _sxy2 += sxy_a; _sxy3 += sxy_a; _epx -= nx; _epy += tsz;
      _vx_x_II += vx_pml_II_a;
      _vx_y_II += vx_pml_II_a;
    }
    /*************************************************************************************/

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/
int esgn_gts2d_24v1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, gxe, gye, iy, tsz, tid, vy_a, mvy_a, py_a, sxy_a;
  int vy_pml_I_a, vy_pml_II_a, vy_pml_III_a, vy_pml_IV_a;
  int gxs_pml_III, gxe_pml_III, gxs_pml_IV, gxe_pml_IV, gys_pml_I, gye_pml_I, gys_pml_II, gye_pml_II;
  register ireal * restrict _vy, * restrict _mvy, * restrict _vyend;
  register ireal * restrict _vy_x_I, * restrict _vy_y_I;
  register ireal * restrict _vy_x_II, * restrict _vy_y_II;
  register ireal * restrict _vy_x_III, * restrict _vy_y_III;
  register ireal * restrict _vy_x_IV, * restrict _vy_y_IV;
  register ireal * restrict _sxy3, * restrict _py3, * restrict _py2;
  register ireal * restrict _py1, * restrict _py0;
  register ireal * restrict _epx, * restrict _epy;
  register ireal lax, lay, dt2, sxy3, sxy2, sxy1, sxy0, dfdx, dfdy, etaxdt, etaydt;

  //register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  //register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  RARR *s_pml;
  RDOM *ld_pml;
  int empty;

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  ld_pml = ((ESGN_TS_PARS*)pars)->ld_pml;
  
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_V1]._dims[0].gs;
  gxe = s[D_V1]._dims[0].ge;
  gys = s[D_V1]._dims[1].gs;
  gye = s[D_V1]._dims[1].ge;
  
  lax = ((ESGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((ESGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((ESGN_TS_PARS*)pars)->dt / 2.0;

  /** pml region I */
  rd_empty(ld_pml+0,D_V1,&empty);
  if (empty) {
    gys_pml_I = gys;
    gye_pml_I = gys-1;
  }
  else {
    s_pml = ld_pml[0]._s;
    gys_pml_I = s_pml[D_V1]._dims[1].gs;
    gye_pml_I = s_pml[D_V1]._dims[1].ge;
  }
  /** pml region II */
  rd_empty(ld_pml+2,D_V1,&empty);
  if (empty) {
    gys_pml_II = gye+1;
    gye_pml_II = gye;
  }
  else {
    s_pml = ld_pml[2]._s;
    gys_pml_II = s_pml[D_V1]._dims[1].gs;
    gye_pml_II = s_pml[D_V1]._dims[1].ge;
  }
  /** pml region III */
  rd_empty(ld_pml+4,D_V1,&empty);
  if (empty) {
    gxs_pml_III = gxs;
    gxe_pml_III = gxs-1;
  }
  else {
    s_pml = ld_pml[4]._s;
    gxs_pml_III = s_pml[D_V1]._dims[0].gs;
    gxe_pml_III = s_pml[D_V1]._dims[0].ge;
  }
 /** pml region IV */
  rd_empty(ld_pml+6,D_V1,&empty);
  if (empty) {
    gxs_pml_IV = gxe+1;
    gxe_pml_IV = gxe;
  }
  else {
    s_pml = ld_pml[6]._s;
    gxs_pml_IV = s_pml[D_V1]._dims[0].gs;
    gxe_pml_IV = s_pml[D_V1]._dims[0].ge;
  }
  
#pragma omp parallel private(tsz,tid,iy,_vy,_vyend,_mvy,_epx,_epy,_sxy3,_py3,_py2,_py1,_py0,sxy3,sxy2,sxy1,sxy0,_vy_x_I,_vy_y_I,_vy_x_II,_vy_y_II,_vy_x_III,_vy_y_III,_vy_x_IV,_vy_y_IV,dfdx,dfdy,etaxdt,etaydt)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif

    /** pml region I *********************************************************************/
#pragma omp single
    {            
      vy_a         = tsz * s[D_V1]._dims[0].n0 - nx;
           
      vy_pml_I_a   = tsz * ld_pml[0]._s[D_V1]._dims[0].n0-nx;
      vy_pml_II_a  = tsz * ld_pml[2]._s[D_V1]._dims[0].n0-nx;
      vy_pml_III_a = tsz * ld_pml[4]._s[D_V1]._dims[0].n0-(gxe_pml_III-gxs_pml_III+1);
      vy_pml_IV_a  = tsz * ld_pml[6]._s[D_V1]._dims[0].n0-(gxe_pml_IV-gxs_pml_IV+1);
      
      mvy_a        = tsz * s[D_MV1]._dims[0].n0 - nx;
      
      py_a         = tsz * s[D_P1]._dims[0].n0 - nx;
      sxy_a        = tsz * s[D_S0]._dims[0].n0 - nx;
    }
#pragma omp barrier
    /** gys == gys_pml_I */
    _vy     = s[D_V1]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    
    s_pml   = ld_pml[0]._s;
    _vy_x_I = s_pml[D_V1]._s + (gxs - s_pml[D_V1]._dims[0].gs) + (gys - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
    
    s_pml   = ld_pml[1]._s;
    _vy_y_I = s_pml[D_V1]._s + (gxs - s_pml[D_V1]._dims[0].gs) + (gys - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
   
    _mvy    = cs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs + tid) * s[D_MV1]._dims[0].n0;
        
    _sxy3   = rs[D_S0]._s + (gxs - s[D_S0]._dims[0].gs) + (gys - s[D_S0]._dims[1].gs + tid) * s[D_S0]._dims[0].n0 + 1;
    _py1    = rs[D_P1]._s + (gxs - s[D_P1]._dims[0].gs) + (gys - s[D_P1]._dims[1].gs + tid) * s[D_P1]._dims[0].n0;
    _py2    = _py1 + s[D_P1]._dims[0].n0;
    _py3    = _py2 + s[D_P1]._dims[0].n0;
    _py0    = _py1 - s[D_P1]._dims[0].n0;
    _epx    = rs[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);              /* 1D */
    _epy    = rs[D_EV[1]]._s + (gys - s[D_EV[1]]._dims[0].gs + tid);        /* 1D */
    
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
   
    for ( iy = gys_pml_I+tid; iy < gye_pml_I+1; iy += tsz )
    {
      sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _vyend = _vy + nx; _vy < _vyend; )
      {
        sxy3 = *_sxy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
        dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vy_x_I) = ((*_vy_x_I) * (1.0f -etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
          (*_vy_y_I) = ((*_vy_y_I) * (1.0f -etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);
          
          *_vy = *_vy_x_I + *_vy_y_I;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vy++; _mvy++;
        _vy_x_I++; _vy_y_I++;
	
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
      _vy += vy_a; _mvy += mvy_a; _sxy3 += sxy_a; _py0 += py_a;
      _py1 += py_a; _py2 += py_a; _py3 += py_a; _epx -= nx; _epy += tsz;
      _vy_x_I += vy_pml_I_a;
      _vy_y_I += vy_pml_I_a;
    }
    /*************************************************************************************/

    /** pml region III, IV and physical region *******************************************/
    /** adjust pointers */
    _vy   += (gye_pml_I + 1 + tid - iy) * s[D_V1 ]._dims[0].n0;
    _mvy  += (gye_pml_I + 1 + tid - iy) * s[D_MV1]._dims[0].n0;
    _sxy3 += (gye_pml_I + 1 + tid - iy) * s[D_S0 ]._dims[0].n0;
    _py2  += (gye_pml_I + 1 + tid - iy) * s[D_P1 ]._dims[0].n0;
    
    _py3 = _py2 + s[D_P1]._dims[0].n0;
    _py1 = _py2 - s[D_P1]._dims[0].n0;
    _py0 = _py1 - s[D_P1]._dims[0].n0;
  
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gye_pml_I + 1 + tid - s[D_EV[1]]._dims[0].gs);        /* 1D */
    
    /** Assumption: pml is thin, so the integer resulted from sum and multiple is not too large */
    s_pml = ld_pml[4]._s;
    _vy_x_III = s_pml[D_V1]._s + (gxs_pml_III - s_pml[D_V1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
   
    s_pml = ld_pml[5]._s;
    _vy_y_III = s_pml[D_V1]._s + (gxs_pml_III - s_pml[D_V1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
    
    s_pml = ld_pml[6]._s;
    _vy_x_IV = s_pml[D_V1]._s + (gxs_pml_IV - s_pml[D_V1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
    
    s_pml = ld_pml[7]._s;
    _vy_y_IV = s_pml[D_V1]._s + (gxs_pml_IV - s_pml[D_V1]._dims[0].gs) + (gye_pml_I + 1 - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
    
    /** adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gye_pml_I+1+tid; iy < gys_pml_II; iy += tsz )
    {
      sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
      //etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      /** pml region III */
      for ( _vyend = _vy + gxe_pml_III-gxs_pml_III+1; _vy < _vyend; )
      {
        sxy3 = *_sxy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
        dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vy_x_III) = ((*_vy_x_III) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
          (*_vy_y_III) = (*_vy_y_III) + dfdy*(*_mvy);

          *_vy = *_vy_x_III + *_vy_y_III;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vy++; _mvy++;
        _vy_x_III++; _vy_y_III++;
	
        sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
      }
      
      /** physical region */
      for ( _vyend = _vy + gxs_pml_IV-gxe_pml_III-1; _vy < _vyend; )
      {
        sxy3 = *_sxy3++;
        _epx++;
        
        dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
        dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
        if(_fwd) {
          *_vy = *_vy + (dfdx + dfdy) * (*_mvy);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vy++; _mvy++;

        sxy0 = sxy1; sxy1 = sxy2; sxy2 = sxy3;
      }
      
      /** pml region IV */
      for ( _vyend = _vy + gxe_pml_IV-gxs_pml_IV+1; _vy < _vyend; )
      {
        sxy3 = *_sxy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
        dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vy_x_IV) = ((*_vy_x_IV) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
          (*_vy_y_IV) = (*_vy_y_IV) + dfdy*(*_mvy);

          *_vy = *_vy_x_IV + *_vy_y_IV;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vy++; _mvy++;
        _vy_x_IV++; _vy_y_IV++;
        
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
      _vy += vy_a; _mvy += mvy_a; _sxy3 += sxy_a; _py0 += py_a;
      _py1 += py_a; _py2 += py_a; _py3 += py_a; _epx -= nx;
      _vy_x_III += vy_pml_III_a; 
      _vy_y_III += vy_pml_III_a;
      _vy_x_IV  += vy_pml_IV_a;
      _vy_y_IV  += vy_pml_IV_a;
    
    }
    /*************************************************************************************/

    /** pml region II ********************************************************************/
    /** adjust pointers */
    _vy   += (gys_pml_II + tid - iy) * s[D_V1 ]._dims[0].n0;
    _mvy  += (gys_pml_II + tid - iy) * s[D_MV1]._dims[0].n0;
    _sxy3 += (gys_pml_II + tid - iy) * s[D_S0 ]._dims[0].n0;
    _py2  += (gys_pml_II + tid - iy) * s[D_P1 ]._dims[0].n0;
    
    _py3 = _py2 + s[D_P1]._dims[0].n0;
    _py1 = _py2 - s[D_P1]._dims[0].n0;
    _py0 = _py1 - s[D_P1]._dims[0].n0;
  
    _epx = s[D_EP[0]]._s + (gxs - s[D_EP[0]]._dims[0].gs);                        /* 1D */
    _epy = s[D_EV[1]]._s + (gys_pml_II - s[D_EV[1]]._dims[0].gs + tid);           /* 1D */
    
    s_pml = ld_pml[2]._s;
    _vy_x_II = s_pml[D_V1]._s + (gxs - s_pml[D_V1]._dims[0].gs) + (gys_pml_II - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
    
    s_pml = ld_pml[3]._s;
    _vy_y_II = s_pml[D_V1]._s + (gxs - s_pml[D_V1]._dims[0].gs) + (gys_pml_II - s_pml[D_V1]._dims[1].gs + tid) * s_pml[D_V1]._dims[0].n0;
 
    /* adjoint formula not derived */
    //_rmpx = rs[D_MP00]._s + (gxs - s[D_MP00]._dims[0].gs) + (gys - s[D_MP00]._dims[1].gs + tid) * s[D_MP00]._dims[0].n0;
    
    for ( iy = gys_pml_II+tid; iy < gye_pml_II+1; iy += tsz )
    {
      sxy2 = _sxy3[-1]; sxy1 = _sxy3[-2]; sxy0 = _sxy3[-3];
      etaydt = (*_epy) * dt2;
      
      /* swap pointers when _fwd = 0 (imaging accumulation), not derived yet */
      /*
        if(!_fwd) { 
        tmp = _px;
        _px=_mpx; 
        _mpx=tmp;
        }
      */
      for ( _vyend = _vy + nx; _vy < _vyend; )
      {
        sxy3 = *_sxy3++;
        etaxdt = (*_epx++) * dt2;
        
        dfdx = (sxy3 - sxy0 + (sxy1 - sxy2) * 27.0) * lax; 
        dfdy = ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay;
        if(_fwd) {
          (*_vy_x_II) = ((*_vy_x_II) * (1.0f - etaxdt) + dfdx*(*_mvy))/(1.0f + etaxdt);
          (*_vy_y_II) = ((*_vy_y_II) * (1.0f - etaydt) + dfdy*(*_mvy))/(1.0f + etaydt);

          *_vy = *_vy_x_II + *_vy_y_II;
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
        }
        else {
          /*
          (*_px) = (*_px) + delta/(*_rmpx++);
          */
          // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
        }
        _vy++; _mvy++;
        _vy_x_II++; _vy_y_II++;
	
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
      _vy += vy_a; _mvy += mvy_a; _sxy3 += sxy_a; _py0 += py_a;
      _py1 += py_a; _py2 += py_a; _py3 += py_a; _epx -= nx; _epy += tsz;
      _vy_x_II += vy_pml_II_a; 
      _vy_y_II += vy_pml_II_a;
    }
    /*************************************************************************************/

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/


/*---- END POINTER BRANCH ----------------------------------------------------*/
#endif


/*********************************** original sgn functions ****************************************/
