/* 
time-step functions
D.S. 02.08.12 
WWS 10.04.11 
After D.S. 08.11.09
After IST 2008-9
********************************************************************************
Implementations of 2-k staggered grid schemes in 1D, 2D, and 3D.
*/

#include "sgn.h"

/*============================================================================*/

#define NSTORE 20000
static ireal _delta[NSTORE];
static ireal _rdelta[NSTORE];
#define VEC

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/
/*
1/24 constant.
*/
#define C24 4.166666666666666666666667e-2
/*----------------------------------------------------------------------------*/

int asg_fts2d_24p01(RDOM *dom, void *pars);
int asg_ftsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars);

int asg_fts2d_24v0(RDOM *dom, void *pars);
int asg_ftsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars);

int asg_fts2d_24v1(RDOM *dom, void *pars);
int asg_ftsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars);
int asg_atsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars);

int asg_fts2d_24(RDOM *dom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return asg_fts2d_24p01(dom, pars);
  if ( iarr == D_V0 ) return asg_fts2d_24v0(dom, pars);
  if ( iarr == D_V1 ) return asg_fts2d_24v1(dom, pars);

  return 0;
}

int asg_ftsm2d_24(RDOM *dom, RDOM *refdom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return asg_ftsm2d_24p01(dom, refdom, pars);
  if ( iarr == D_V0 ) return asg_ftsm2d_24v0(dom, refdom, pars);
  if ( iarr == D_V1 ) return asg_ftsm2d_24v1(dom, refdom, pars);

  return 0;
}

int asg_atsm2d_24(RDOM *dom, RDOM *refdom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return asg_atsm2d_24p01(dom, refdom, pars);
  if ( iarr == D_V0 ) return asg_atsm2d_24v0(dom, refdom, pars);
  if ( iarr == D_V1 ) return asg_atsm2d_24v1(dom, refdom, pars);

  return 0;
}

/*----------------------------------------------------------------------------*/

#ifdef VEC

/*-------------- SIM FUNCS --------------------------------------*/
/*-------------- fts_p01sv ----------------------------------*/
int asg_fts2d_24p01(RDOM *dom, void *pars) {
  int nx, ny, gxs, gys, iy, px_nx0, py_nx0, mpx_nx0, vx_nx0, vy_nx0;
  register ireal * restrict _px, * restrict _py,
    * restrict _mpx, * restrict _epx, * restrict _epy, * restrict _vx3,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0, * restrict _py_ol, * restrict _py_or, * restrict _py_l, * restrict _py_r, * restrict _px_ol;

  register ireal lax, lay, dt2, etaydt, rfac0, fac0;
  register int i, lbcx,rbcx,lbcy,rbcy, k;    /* auxiliary variables */
  
  RARR *s;

  s = dom->_s;
 
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  px_nx0 = s[D_P0]._dims[0].n0;
  py_nx0 = s[D_P1]._dims[0].n0;
  mpx_nx0 = s[D_MP0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;
 

  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;

  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  
  k = ((SGN_TS_PARS*)pars)->k;

  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];
  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0;
  _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _px_ol = _px;

  _py_ol = _py;
  _py_l = _py_ol - py_nx0;
  _py_r = _py + ny * py_nx0;
  _py_or = _py_r - py_nx0;
  
  _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * mpx_nx0;

  _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0 + 1;
  _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  _vy3 = _vy2 + vy_nx0;
  _vy1 = _vy2 - vy_nx0;
  _vy0 = _vy1 - vy_nx0;
  _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);      /* 1D */
  _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);      /* 1D */
  
  for ( iy = 0; iy < ny; iy ++) {
    etaydt = (*_epy) * dt2;
    rfac0 = 1.0f / (1.0f + etaydt);
    fac0 = (1.0f - etaydt) * rfac0;

    for ( i = 0; i < nx; ++i )  
      _delta[i] = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
		   (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay);
    
    
    for ( i = 0; i < nx; ++i ) { 
      _px[i] = (_px[i] * (1.0f - _epx[i] * dt2) + _delta[i]*_mpx[i]) / (1.0f + _epx[i] * dt2);
      _py[i] = _py[i] * fac0 + _delta[i]*_mpx[i] * rfac0;	 
    }
    
    _px  += px_nx0;
    _py  += py_nx0; 

    _mpx += mpx_nx0;
    _vx3 += vx_nx0; 
    _vy0 += vy_nx0;
    _vy1 += vy_nx0; 
    _vy2 += vy_nx0; 
    _vy3 += vy_nx0; 
    _epy ++;
  }

  /* reflection boundary for _px  */
  if (lbcx){
    /* -- range [gxs0 : gxs-1] ) -- */
    _px = _px_ol;
    for ( iy = 0; iy < ny; iy ++) {
      for (i = -k; i < -1; i++)
	_px[i] = -_px[-2 - i];
      _px[-1] = REAL_ZERO;
      _px  += px_nx0;
    }
  }
  if (rbcx){
    /* -- range [gxe+1 : gxe0] -- */
    _px = _px_ol + nx;
    for ( iy = 0; iy < ny; iy ++) {
      _px[0] = REAL_ZERO;
      for (i = 1; i < k; i++)
	_px[i] = -_px[-i];
      _px  += px_nx0;
    }
  }

  /* reflection boundary for _py  */
  if (lbcy) {
    /* -- gys-1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_l[i] = REAL_ZERO; 
    _py_l  -= py_nx0; 
    /* -- range [gys0 : gys-2] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_py_l[i] = -_py_ol[i];
      _py_l  -= py_nx0; 
      _py_ol += py_nx0; 
    }
  }
  
  if (rbcy) {
    /* -- gye+1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_r[i] = REAL_ZERO; 
    _py_r  += py_nx0;     
    /* -- range [gye+2 : gye0] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_py_r[i] = -_py_or[i];
      _py_r  += py_nx0; 
      _py_or -= py_nx0; 
    }
  }
  
  return 0;
}

/*--------------- fts_v0 ----------------------------------*/
int asg_fts2d_24v0(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, iy, px_nx0, mvx_nx0, vx_nx0;
  register ireal * restrict _vx, * restrict _mvx,
    * restrict _evx, * restrict _px3, * restrict _vx_ol;
  register ireal lax, dt2;
  register int i, k, km1, lbcx, rbcx; /* auxiliary variables */ 
  RARR *s;
     
  s = dom->_s;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;// / 2;
  km1 = k - 1;

  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];

  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  /* px_nx0, mvx_nx0, vx_nx0 */
  px_nx0 = s[D_P0]._dims[0].n0;
  mvx_nx0 = s[D_MV0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;
  
  gxs = s[D_V0]._dims[0].gs;
  gys = s[D_V0]._dims[1].gs;

   _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0;
   _vx_ol = _vx;

  _mvx = s[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * mvx_nx0;
  _px3 = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0 + 2;
  _evx = s[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */

  for ( iy = 0; iy < ny; iy ++ ) {            
    for ( i = 0;  i < nx; ++i ) 
	_vx[i] = (_vx[i] * (1.0f - _evx[i] * dt2) + (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * 27.0f) * lax * _mvx[i]) / (1.0f + _evx[i] * dt2);

    _vx += vx_nx0;
    _mvx += mvx_nx0;
    _px3 += px_nx0;

  } 

  /* reflection boundary for _vx */ 
  // km1 = k - 1 
  if(lbcx) {
    /* -- range [gxs0 : gxs-1] -- */
    _vx = _vx_ol;
    for ( iy = 0; iy < ny; iy ++ ) {            
      for (i = -km1; i < 0; i++)
	_vx[i] = _vx[-1 - i];
      _vx += vx_nx0;
    }
  }
  if(rbcx) {
    /* -- range [gxe+1 : gxe0] -- */
    _vx = _vx_ol + nx;
    for ( iy = 0; iy < ny; iy ++ ) {            
       for (i = 0; i < km1; i++)
	 _vx[i] = _vx[-1 - i];
       _vx += vx_nx0;
    }
  }

  return 0;
}

/*------------- fts_v1 ----------------------------------*/

int asg_fts2d_24v1(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, iy, py_nx0, mvy_nx0, vy_nx0;
  register ireal * restrict _vy, * restrict _mvy, * restrict _evy,
    * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0, * restrict _vy_ol, * restrict _vy_or, * restrict _vy_l, * restrict _vy_r;
  register ireal lay, dt2, etaydt, rfac0, fac0;
  register int i, k, lbcy, rbcy;

  RARR *s;

  s = dom->_s;
    
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  /* py_nx0, mvy_nx0, vy_nx0 */
  py_nx0 = s[D_P1]._dims[0].n0;
  mvy_nx0 = s[D_MV1]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;
    
  gxs = s[D_V1]._dims[0].gs;
  gys = s[D_V1]._dims[1].gs;
    
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;// / 2;

  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  
  _vy_ol = _vy;
  _vy_l = _vy_ol - vy_nx0;
  _vy_r = _vy + ny * vy_nx0;
  _vy_or = _vy_r - vy_nx0;

  _mvy = s[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * mvy_nx0;
  _py1 = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _py2 = _py1 + py_nx0;
  _py3 = _py2 + py_nx0;
  _py0 = _py1 - py_nx0;
  _evy = s[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs);        /* 1D */
  
  for ( iy = 0; iy < ny; iy ++) {            
      
    etaydt = (*_evy) * dt2;
    rfac0 = 1.0f / (1.0f + etaydt);
    fac0 = (1.0f - etaydt) * rfac0;
         
    for ( i = 0; i < nx; ++i ) 
      _vy[i] = _vy[i] * fac0 + 
	(_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * 27.0f) * lay * _mvy[i] * rfac0;
    
    _vy += vy_nx0;
    _mvy += mvy_nx0;
    _py0 += py_nx0;
    _py1 += py_nx0;
    _py2 += py_nx0;
    _py3 += py_nx0;
    _evy ++;
    
  }

  /* reflection boundary for _vy  */
  if (lbcy) {
    /* -- range [gys0 : gys-1] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_vy_l[i] = _vy_ol[i];
      _vy_l  -= vy_nx0; 
      _vy_ol += vy_nx0; 
    }
  }
  if (rbcy) {
    /* -- range [gye+1 : gye0] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_vy_r[i] = _vy_or[i];
      _vy_r  += vy_nx0; 
      _vy_or -= vy_nx0; 
    }
  }
 
  return 0;
}

/*------------ Born-SIM FUNCS -----------------------------------*/
/*-------------- ftsm_p01sv -------------------------*/
int asg_ftsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars) {
  int nx, ny, gxs, gys, iy, px_nx0, py_nx0, mpx_nx0, vx_nx0, vy_nx0;
  register ireal * restrict _px, * restrict _py,
    * restrict _mpx, * restrict _rmpx, * restrict _epx, * restrict _epy, * restrict _vx3,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0, * restrict _rvx3,
    * restrict _rvy3, * restrict _rvy2, * restrict _rvy1, * restrict _rvy0, * restrict _py_ol, * restrict _py_or, * restrict _py_l, * restrict _py_r, * restrict _px_ol;

  register ireal lax, lay, dt2, etaydt, rfac0, fac0;
  register int i, lbcx,rbcx,lbcy,rbcy, k;    /* auxiliary variables */

  RARR *s, *rs;

  s = dom->_s;
  rs = rdom->_s;
 
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
 
  px_nx0 = s[D_P0]._dims[0].n0;
  py_nx0 = s[D_P1]._dims[0].n0;
  mpx_nx0 = s[D_MP0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;


  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;

  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  
  k = ((SGN_TS_PARS*)pars)->k;
 
  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];
  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0;
  _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _px_ol = _px;
  
  _py_ol = _py;
  _py_l = _py_ol - py_nx0;
  _py_r = _py + ny * py_nx0;
  _py_or = _py_r - py_nx0;
  
  _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * mpx_nx0;
  _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0 + 1;
  _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  _vy3 = _vy2 + vy_nx0;
  _vy1 = _vy2 - vy_nx0;
  _vy0 = _vy1 - vy_nx0;

  _rmpx = rs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * mpx_nx0;
  _rvx3 = rs[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0 + 1;
  _rvy2 = rs[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  _rvy3 = _rvy2 + vy_nx0;
  _rvy1 = _rvy2 - vy_nx0;
  _rvy0 = _rvy1 - vy_nx0;

  _epx = rs[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);      /* 1D */
  _epy = rs[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);      /* 1D */
  
  for ( iy = 0; iy < ny; iy ++) {
    etaydt = (*_epy) * dt2;
    rfac0 = 1.0f / (1.0f + etaydt);
    fac0 = (1.0f - etaydt) * rfac0;

    for ( i = 0; i < nx; ++i )  {
      _delta[i] = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
		   (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay);

      _rdelta[i] = ((_rvx3[i] - _rvx3[i-3] + (_rvx3[i-2] - _rvx3[i-1]) * (ireal)27.0) * lax + 
		   (_rvy3[i] - _rvy0[i  ] + (_rvy1[i  ] - _rvy2[i  ]) * (ireal)27.0) * lay);
    }
    
    for ( i = 0; i < nx; ++i ) { 
      _px[i] = (_px[i] * (1.0f - _epx[i] * dt2) + 
		_delta[i]*_rmpx[i] +
		_rdelta[i]*_mpx[i]
		) / (1.0f + _epx[i] * dt2);

      _py[i] = _py[i] * fac0 + (_delta[i]*_rmpx[i] + _rdelta[i]*_mpx[i]) * rfac0;
    }
    
    _px  += px_nx0;
    _py  += py_nx0; 

    _mpx += mpx_nx0;
    _vx3 += vx_nx0; 
    _vy0 += vy_nx0;
    _vy1 += vy_nx0; 
    _vy2 += vy_nx0; 
    _vy3 += vy_nx0; 
    _epy ++;

    _rmpx += mpx_nx0;
    _rvx3 += vx_nx0; 
    _rvy0 += vy_nx0;
    _rvy1 += vy_nx0; 
    _rvy2 += vy_nx0; 
    _rvy3 += vy_nx0; 
  }

  /* reflection boundary for _px  */
  if (lbcx){
    /* -- range [gxs0 : gxs-1] ) -- */
    _px = _px_ol;
    for ( iy = 0; iy < ny; iy ++) {
      for (i = -k; i < -1; i++)
	_px[i] = -_px[-2 - i];
      _px[-1] = REAL_ZERO;
      _px  += px_nx0;
    }
  }
  if (rbcx){
    /* -- range [gxe+1 : gxe0] -- */
    _px = _px_ol + nx;
    for ( iy = 0; iy < ny; iy ++) {
      _px[0] = REAL_ZERO;
      for (i = 1; i < k; i++)
	_px[i] = -_px[-i];
      _px  += px_nx0;
    }
  }

  /* reflection boundary for _py  */
  if (lbcy) {
    /* -- gys-1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_l[i] = REAL_ZERO; 
    _py_l  -= py_nx0; 
    /* -- range [gys0 : gys-2] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_py_l[i] = -_py_ol[i];
      _py_l  -= py_nx0; 
      _py_ol += py_nx0; 
    }
  }
  
  if (rbcy) {
    /* -- gye+1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_r[i] = REAL_ZERO; 
    _py_r  += py_nx0;     
    /* -- range [gye+2 : gye0] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_py_r[i] = -_py_or[i];
      _py_r  += py_nx0; 
      _py_or -= py_nx0; 
    }
  }


  return 0;
}

/*--------------- ftsm_v0 ----------------------------------*/
int asg_ftsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars)
{
  int nx, ny, gxs, gys, iy, px_nx0, mvx_nx0, vx_nx0;
  register ireal * restrict _vx, * restrict _mvx, * restrict _rmvx,
    * restrict _evx, * restrict _px3, * restrict _rpx3, * restrict _vx_ol;
  register ireal lax, dt2;
  register int i, k, km1, lbcx, rbcx; /* auxiliary variables */ 
  RARR *s, *rs;
     
  s = dom->_s;
  rs = rdom->_s;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;
  km1 = k - 1;

  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];

  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  /* px_nx0, mvx_nx0, vx_nx0 */
  px_nx0 = s[D_P0]._dims[0].n0;
  mvx_nx0 = s[D_MV0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;

  gxs = s[D_V0]._dims[0].gs;
  gys = s[D_V0]._dims[1].gs;
    
   _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0;
   _vx_ol = _vx;

  _mvx = s[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * mvx_nx0;
  _px3 = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0 + 2;
  _evx = rs[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */

  _rmvx = rs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * mvx_nx0;      
  _rpx3 = rs[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0 + 2;

  for ( iy = 0; iy < ny; iy ++ ) {            
    for ( i = 0;  i < nx; ++i ) 
	_vx[i] = (_vx[i] * (1.0f - _evx[i] * dt2) +
		  (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * 27.0f) * lax * _rmvx[i] +
		  (_rpx3[i] - _rpx3[i-3] + (_rpx3[i-2] - _rpx3[i-1]) * 27.0f) * lax * _mvx[i]
		  ) / (1.0f + _evx[i] * dt2);

    _vx += vx_nx0;
    _mvx += mvx_nx0;
    _px3 += px_nx0;

    _rmvx += mvx_nx0;
    _rpx3 += px_nx0;
  } 

  /* reflection boundary for _vx */ 
  // km1 = k - 1 
  if(lbcx) {
    /* -- range [gxs0 : gxs-1] -- */
    _vx = _vx_ol;
    for ( iy = 0; iy < ny; iy ++ ) {            
      for (i = -km1; i < 0; i++)
	_vx[i] = _vx[-1 - i];
      _vx += vx_nx0;
    }
  }
  if(rbcx) {
    /* -- range [gxe+1 : gxe0] -- */
    _vx = _vx_ol + nx;
    for ( iy = 0; iy < ny; iy ++ ) {            
       for (i = 0; i < km1; i++)
	 _vx[i] = _vx[-1 - i];
       _vx += vx_nx0;
    }
  }

  return 0;
}

/*------------- ftsm_v1 ----------------------------------*/
int asg_ftsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars)
{
  int nx, ny, gxs, gys, iy, py_nx0, mvy_nx0, vy_nx0;
  register ireal * restrict _vy, * restrict _mvy, 
    * restrict _rmvy, * restrict _evy,
    * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0, 
    * restrict _rpy3, * restrict _rpy2, * restrict _rpy1, * restrict _rpy0, 
    * restrict _vy_ol, * restrict _vy_or, * restrict _vy_l, * restrict _vy_r;

  register ireal lay, dt2, etaydt, rfac0, fac0;
  register int i, k, lbcy, rbcy;

  RARR *s, *rs;

  s = dom->_s;
  rs = rdom->_s;
    
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  /* py_nx0, mvy_nx0, vy_nx0 */
  py_nx0 = s[D_P1]._dims[0].n0;
  mvy_nx0 = s[D_MV1]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;
    
  gxs = s[D_V1]._dims[0].gs;
  gys = s[D_V1]._dims[1].gs;
    
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;

  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;

  _vy_ol = _vy;
  _vy_l = _vy_ol - vy_nx0;
  _vy_r = _vy + ny * vy_nx0;
  _vy_or = _vy_r - vy_nx0;

  _mvy = s[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * mvy_nx0;
  _py1 = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _py2 = _py1 + py_nx0;
  _py3 = _py2 + py_nx0;
  _py0 = _py1 - py_nx0;
  _evy = rs[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs);        /* 1D */
  
  _rmvy = rs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * mvy_nx0;
  _rpy1 = rs[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _rpy2 = _rpy1 + py_nx0;
  _rpy3 = _rpy2 + py_nx0;
  _rpy0 = _rpy1 - py_nx0;

  for ( iy = 0; iy < ny; iy ++) {            
      
    etaydt = (*_evy) * dt2;
    rfac0 = 1.0f / (1.0f + etaydt);
    fac0 = (1.0f - etaydt) * rfac0;
    for ( i = 0; i < nx; ++i ) 
      _vy[i] = _vy[i] * fac0 + 
	((_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * 27.0f) * lay * _rmvy[i] + 
	 (_rpy3[i] - _rpy0[i] + (_rpy1[i] - _rpy2[i]) * 27.0f) * lay * _mvy[i]) * rfac0;
    
    _vy += vy_nx0;
    _mvy += mvy_nx0;
    _py0 += py_nx0;
    _py1 += py_nx0;
    _py2 += py_nx0;
    _py3 += py_nx0;
    _evy ++;
    
    _rmvy += mvy_nx0;
    _rpy0 += py_nx0;
    _rpy1 += py_nx0;
    _rpy2 += py_nx0;
    _rpy3 += py_nx0;
  }

  /* reflection boundary for _vy  */
  if (lbcy) {
    /* -- range [gys0 : gys-1] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_vy_l[i] = _vy_ol[i];
      _vy_l  -= vy_nx0; 
      _vy_ol += vy_nx0; 
    }
  }
  if (rbcy) {
    /* -- range [gye+1 : gye0] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i )  
	_vy_r[i] = _vy_or[i];
      _vy_r  += vy_nx0; 
      _vy_or -= vy_nx0; 
    }
  }

  return 0;
}

/*-------------- ADJ FUNCS --------------------------------------*/
/*-------------- atsm_p01sv -------------------------*/
int asg_atsm2d_24p01(RDOM *dom, RDOM *rdom, void *pars) {
  int nx, ny, gxs, gys, iy, px_nx0, py_nx0, mpx_nx0, vx_nx0, vy_nx0;
  register ireal * restrict _px, * restrict _py,
    * restrict _mpx, * restrict _rmpx, * restrict _epx, * restrict _epy, * restrict _vx3,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0, * restrict _rvx3,
    * restrict _rvy3, * restrict _rvy2, * restrict _rvy1, * restrict _rvy0, * restrict _py_ol, * restrict _py_or, * restrict _py_l, * restrict _py_r, * restrict _px_ol;

  register ireal lax, lay, dt2, etaydt;
  register ireal fac0, fac0y, rfac0, rfac1, facx, facy, fac1, fac2;
  register int i, lbcx,rbcx,lbcy,rbcy, k;    /* auxiliary variables */

  RARR *s, *rs;

  s = dom->_s;
  rs = rdom->_s;
 
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  px_nx0 = s[D_P0]._dims[0].n0;
  py_nx0 = s[D_P1]._dims[0].n0;
  mpx_nx0 = s[D_MP0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;


  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
  
  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;

  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  
  k = ((SGN_TS_PARS*)pars)->k;

  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];
  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0;
  _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _px_ol = _px;
  
  _py_ol = _py;
  _py_l = _py_ol - py_nx0;
  _py_r = _py + ny * py_nx0;
  _py_or = _py_r - py_nx0;
  
  _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * mpx_nx0;
  _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0 + 1;
  _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  _vy3 = _vy2 + vy_nx0;
  _vy1 = _vy2 - vy_nx0;
  _vy0 = _vy1 - vy_nx0;

  _rmpx = rs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * mpx_nx0;
  _rvx3 = rs[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0 + 1;
  _rvy2 = rs[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;
  _rvy3 = _rvy2 + vy_nx0;
  _rvy1 = _rvy2 - vy_nx0;
  _rvy0 = _rvy1 - vy_nx0;

  _epx = rs[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);      /* 1D */
  _epy = rs[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);      /* 1D */
  

  /* adjoint reflection bcs, _py  */
  if (rbcy) {
    /* -- gye+1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_r[i] = REAL_ZERO;
    _py_r  += py_nx0; 
    /* -- range [gye+2 : gye0] */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i ) { 
	_py_or[i] += -_py_r[i];
	_py_r[i] = REAL_ZERO;
      }
      _py_r  += py_nx0; 
      _py_or -= py_nx0; 
    }
  }
  if (lbcy) {
    /* -- gys-1 -- */
    for ( i = 0; i < nx; ++i )  
      _py_l[i] = REAL_ZERO; 
    _py_l  -= py_nx0; 
    /* -- range [gys0 : gys-2] */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i ) { 
	_py_ol[i] += -_py_l[i];
	_py_l[i] = REAL_ZERO;
      }
      _py_l  -= py_nx0; 
      _py_ol += py_nx0; 
    }
  }

  /* adjoint reflection bcs, _px */
  if (rbcx) {
    /* -- range [gxe+1 : gxe0]  */
    _px = _px_ol + nx;
    for ( iy = 0; iy < ny; iy ++) {
      for (i = 1; i < k; i++){
	_px[-i] += -_px[i];
	_px[i] = REAL_ZERO;
      }
      _px[0] = REAL_ZERO;
      _px  += px_nx0;
    }
  }
  if (lbcx) {
    /* -- range [gxs0 : gxs-1] */
    _px = _px_ol;
    for ( iy = 0; iy < ny; iy ++) {
      for (i = -k; i < -1; i++){
	_px[-2 - i] += -_px[i];
	_px[i] = REAL_ZERO;
      }
      _px[-1] = REAL_ZERO;
      _px  += px_nx0;
    }
  }
  
  /* adjoint deriv, vx,vy->px,py  */ 
  _px = _px_ol;
  for ( iy = 0; iy < ny; iy ++) {
    etaydt = (*_epy) * dt2;
    rfac1= 1.0f/(1.0f + etaydt);
    fac0y = (1.0f - etaydt) * rfac1;

    // Note: in ftsm_p01:
    // _px[i] += -_px[i] + (_px[i] * (1.0f - _epx[i] * dt2) + 
    //  	      _delta[i]*_rmpx[i] +
    //	              _rdelta[i]*_mpx[i]
    //		   ) / (1.0f + _epx[i] * dt2);
    // += (-2*_epx[i]*dt2/(1.0 + _epx[i]*dt2)) _px[i] 
    //     + _delta[i] * _rmpx[i]/(1.0 + _epx[i]*dt2)
    //     + _rdelta[i] * _mpx[i]/(1.0 + _epx[i]*dt2)    
    //
    // _delta[i] = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay);
    // ..................
    
    for ( i = 0; i < nx; ++i )  {
      rfac0= 1/(1.0f + _epx[i]*dt2);
      fac0  = (1.0f - _epx[i] * dt2) * rfac0;
      fac1 = lax * _rmpx[i];
      fac2 = lay * _rmpx[i];
      _rdelta[i] = ((_rvx3[i] - _rvx3[i-3] + (_rvx3[i-2] - _rvx3[i-1]) * (ireal)27.0) * lax + (_rvy3[i] - _rvy0[i  ] + (_rvy1[i  ] - _rvy2[i  ]) * (ireal)27.0) * lay);

      // step 1: adjoint state deriv, vx,vy->py
      facy = fac1 * rfac1;
      _vx3[i] += facy * _py[i];
      _vx3[i-3] += -facy * _py[i];
      facy *= (ireal) 27.0;
      _vx3[i-2] += facy * _py[i];
      _vx3[i-1] += -facy * _py[i];

      facy = fac2 * rfac1;      
      _vy3[i]   += facy * _py[i];
      _vy0[i]   += -facy * _py[i];
      facy *= (ireal) 27.0;      
      _vy1[i]   += facy * _py[i];
      _vy2[i]   += -facy * _py[i];
      
      // step 2: adjoint control deriv, vz,vx->py
      _mpx[i] += _rdelta[i] * _py[i] * rfac1;
      
      _py[i] = fac0y * _py[i];

      // step 3: adjoint state deriv, vx,vy->px
      facx = fac1 * rfac0;
      _vx3[i] += facx * _px[i];
      _vx3[i-3] += -facx * _px[i] ;
      facx *= (ireal) 27.0;
      _vx3[i-2] += facx * _px[i];
      _vx3[i-1] += -facx * _px[i];
      
      facx = fac2 * rfac0;
      _vy3[i]   += facx * _px[i];
      _vy0[i]   += -facx * _px[i];
      facx *= (ireal) 27.0;
      _vy1[i]   += facx * _px[i];
      _vy2[i]   += -facx * _px[i];
      
      // step 4: adjoint control deriv, vz,vx->px
      _mpx[i] += _rdelta[i] * _px[i] * rfac0;
      _px[i] = fac0 * _px[i];
	
    }
     
    _px  += px_nx0;
    _py  += py_nx0; 

    _mpx += mpx_nx0;
    _vx3 += vx_nx0; 
    _vy0 += vy_nx0;
    _vy1 += vy_nx0; 
    _vy2 += vy_nx0; 
    _vy3 += vy_nx0; 
    _epy ++;

    _rmpx += mpx_nx0;
    _rvx3 += vx_nx0; 
    _rvy0 += vy_nx0;
    _rvy1 += vy_nx0; 
    _rvy2 += vy_nx0; 
    _rvy3 += vy_nx0; 
  }
    
  return 0;
}

/*--------------- atsm_v0 ----------------------------------*/
int asg_atsm2d_24v0(RDOM *dom, RDOM *rdom, void *pars)
{
  int nx, ny, gxs, gys, iy, px_nx0, mvx_nx0, vx_nx0;
  register ireal * restrict _vx, * restrict _mvx, * restrict _rmvx,
    * restrict _evx, * restrict _px3, * restrict _rpx3, * restrict _vx_ol;
  register ireal lax, dt2;
  register int i, k, km1, lbcx, rbcx; /* auxiliary variables */
  register ireal rfac0, fac0, fac1, fac2, fac3;

  RARR *s, *rs;
     
  s = dom->_s;
  rs = rdom->_s;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;
  km1 = k - 1;

  lbcx = ((SGN_TS_PARS*)pars)->lbc[0];
  rbcx = ((SGN_TS_PARS*)pars)->rbc[0];

  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;

  /* px_nx0, mvx_nx0, vx_nx0 */
  px_nx0 = s[D_P0]._dims[0].n0;
  mvx_nx0 = s[D_MV0]._dims[0].n0;
  vx_nx0 = s[D_V0]._dims[0].n0;

  gxs = s[D_V0]._dims[0].gs;
  gys = s[D_V0]._dims[1].gs;
    
   _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * vx_nx0;
   _vx_ol = _vx;

  _mvx = s[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * mvx_nx0;
  _px3 = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0 + 2;
  _evx = rs[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */

  _rmvx = rs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * mvx_nx0;      
  _rpx3 = rs[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * px_nx0 + 2;


  /* adjoint reflection bcs, _vx */
  if(rbcx) {
    /* -- range [gxe+1 : gxe0]  */
    _vx = _vx_ol + nx;
    for ( iy = 0; iy < ny; iy ++ ) {            
      for (i = 0; i < km1; i++){
	_vx[-1 - i] += _vx[i];
	_vx[i] = REAL_ZERO;
      }
      _vx += vx_nx0;
    }
  }
  if(lbcx) {
    /* -- range [gxs0 : gxs-1] -- */
    _vx = _vx_ol;
    for ( iy = 0; iy < ny; iy ++ ) {            
      for (i = -km1; i < 0; i++){
	_vx[-1 - i] += _vx[i];
	_vx[i] = REAL_ZERO;
      }
      _vx += vx_nx0;
    }
  }

  //  adjoint deriv, px->vx
  _vx = _vx_ol;
  for ( iy = 0; iy < ny; iy ++ ) {            

    for ( i = 0;  i < nx; ++i ) { 
      rfac0 = 1.0f / (1.0f + _evx[i] * dt2);
      fac0 = (1.0f - _evx[i] * dt2) * rfac0;
      fac1 = lax * rfac0;
      fac2 = fac1 * _rmvx[i];
      fac3 = (ireal) 27.0 * fac2;
      // step 1: adjoint state deriv, px->vx
      _px3[i] += fac2 * _vx[i];
      _px3[i-3] += -fac2 * _vx[i];
      _px3[i-2] += fac3 * _vx[i];
      _px3[i-1] += -fac3 * _vx[i];
      // step 2: adjoint control deriv, px->vx
      _mvx[i] += _vx[i] * (_rpx3[i] - _rpx3[i-3] + (_rpx3[i-2] - _rpx3[i-1]) * 27.0f) * fac1;
      _vx[i] = fac0 * _vx[i];
    }
    
    _vx += vx_nx0;
    _mvx += mvx_nx0;
    _px3 += px_nx0;

    _rmvx += mvx_nx0;
    _rpx3 += px_nx0;
  } 

  return 0;
}

/*------------- atsm_v1 ----------------------------------*/
int asg_atsm2d_24v1(RDOM *dom, RDOM *rdom, void *pars)
{
  int nx, ny, gxs, gys, iy, py_nx0, mvy_nx0, vy_nx0;
  register ireal * restrict _vy, * restrict _mvy, 
    * restrict _rmvy, * restrict _evy,
    * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0, 
    * restrict _rpy3, * restrict _rpy2, * restrict _rpy1, * restrict _rpy0, 
    * restrict _vy_ol, * restrict _vy_or, * restrict _vy_l, * restrict _vy_r;

  register ireal lay, dt2, etaydt;
  register int i, k, lbcy, rbcy;
  register ireal rfac0, fac0, fac1, fac2, fac3;

  RARR *s, *rs;

  s = dom->_s;
  rs = rdom->_s;
    
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  /* py_nx0, mvy_nx0, vy_nx0 */
  py_nx0 = s[D_P1]._dims[0].n0;
  mvy_nx0 = s[D_MV1]._dims[0].n0;
  vy_nx0 = s[D_V1]._dims[0].n0;
  
  gxs = s[D_V1]._dims[0].gs;
  gys = s[D_V1]._dims[1].gs;
    
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  k = ((SGN_TS_PARS*)pars)->k;

  lbcy = ((SGN_TS_PARS*)pars)->lbc[1];
  rbcy = ((SGN_TS_PARS*)pars)->rbc[1];

  _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * vy_nx0;

  _vy_ol = _vy;
  _vy_l = _vy_ol - vy_nx0;
  _vy_r = _vy + ny * vy_nx0;
  _vy_or = _vy_r - vy_nx0;

  _mvy = s[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * mvy_nx0;
  _py1 = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _py2 = _py1 + py_nx0;
  _py3 = _py2 + py_nx0;
  _py0 = _py1 - py_nx0;
  _evy = rs[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs);        /* 1D */
  
  _rmvy = rs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * mvy_nx0;
  _rpy1 = rs[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * py_nx0;
  _rpy2 = _rpy1 + py_nx0;
  _rpy3 = _rpy2 + py_nx0;
  _rpy0 = _rpy1 - py_nx0;

  /* adjoint reflection bcs, vy */
  if (rbcy) {
    /* -- range [gye+1 : gye0] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i ) {
	_vy_or[i] += _vy_r[i];
	_vy_r[i] = REAL_ZERO; 
      }
      _vy_r  += vy_nx0; 
      _vy_or -= vy_nx0; 
    }
  }
  if (lbcy) {
    /* -- range [gys0 : gys-1] -- */
    for ( iy = 0; iy < k-1; iy++) { 
      for ( i = 0; i < nx; ++i ){
	_vy_ol[i] += _vy_l[i]; 
	_vy_l[i] = REAL_ZERO;
      }
      _vy_l  -= vy_nx0; 
      _vy_ol += vy_nx0; 
    }
  }

  // adjoint deriv, py->vy
  for ( iy = 0; iy < ny; iy ++) {            
      
    etaydt = (*_evy) * dt2;
    rfac0 = 1.0f / (1.0f + etaydt);
    fac0 = (1.0f - etaydt) * rfac0;

    for ( i = 0; i < nx; ++i ) {
      fac1 = lay * rfac0;
      fac2 = fac1 * _rmvy[i];
      fac3 = (ireal) 27.0 * fac2;
      // step 1: adjoint state deriv, py->vy
      _py3[i] += fac2 * _vy[i];
      _py0[i] += -fac2 * _vy[i];
      _py1[i] += fac3 * _vy[i];
      _py2[i] += -fac3 * _vy[i];
      // step 2: adjoint control deriv, py->vy      
      _mvy[i] += _vy[i] * (_rpy3[i] - _rpy0[i] + (_rpy1[i] - _rpy2[i]) * 27.0f) * fac1;
      _vy[i] = fac0 * _vy[i];
    }
    
    _vy += vy_nx0;
    _mvy += mvy_nx0;
    _py0 += py_nx0;
    _py1 += py_nx0;
    _py2 += py_nx0;
    _py3 += py_nx0;
    _evy ++;
    
    _rmvy += mvy_nx0;
    _rpy0 += py_nx0;
    _rpy1 += py_nx0;
    _rpy2 += py_nx0;
    _rpy3 += py_nx0;
  }

  return 0;
}


/*---- END VECTOR BRANCH -----------------------------------------------------*/

#else 
/*---- BEGIN POINTER BRANCH --------------------------------------------------*/

/*---- END POINTER BRANCH ----------------------------------------------------*/
#endif
