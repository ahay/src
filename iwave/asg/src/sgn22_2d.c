/* 
   sgn22_2d.c
   Igor Terentyev.
   Prototype - Tanya Vdovina.
   ********************************************************************************
   Implementation of 2-2 scheme in 2D.
*/
/*============================================================================*/

#include "utils.h"
#include "sgn.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22p01(RDOM *, void *);
int sgn_ts2d_22p0(RDOM *, void *);
int sgn_ts2d_22p1(RDOM *, void *);
int sgn_ts2d_22p(RDOM *, void *);
int sgn_ts2d_22v0(RDOM *, void *);
int sgn_ts2d_22v1(RDOM *, void *);

int sgn_ts2d_22(RDOM *dom, int iarr, void *pars)
{
  int e0, e1;                   /* short-hands */

  if ( iarr == D_P0 )           /* covers D_P0, D_P1 in one function */
    {
      /*
      if ( ((SGN_TS_PARS*)pars)->psingle ) return sg_ts2d_22p(dom, &(((SGN_TS_PARS*)pars)->sgpars));
      else
      */
	{
	  e0 = ((SGN_TS_PARS*)pars)->eflags[0];
	  e1 = ((SGN_TS_PARS*)pars)->eflags[1];
	  if ( e0 && e1 ) return sgn_ts2d_22p01(dom, pars);
	  else if ( e0 )  return sgn_ts2d_22p0(dom, pars);
	  else if ( e1 )  return sgn_ts2d_22p1(dom, pars);
	  else            return sgn_ts2d_22p(dom, pars);
	}
    }
  if ( iarr == D_V0 )
    {
      /*
      if ( ((SGN_TS_PARS*)pars)->eflags[0] ) 
      */
      return sgn_ts2d_22v0(dom, pars);
      //      else return sg_ts2d_22v0(dom, &(((SGN_TS_PARS*)pars)->sgpars));
    }
  if ( iarr == D_V1 )
    {
      /*
      if ( ((SGN_TS_PARS*)pars)->eflags[1] )
      */
      return sgn_ts2d_22v1(dom, pars);
      //      else return sg_ts2d_22v1(dom, &(((SGN_TS_PARS*)pars)->sgpars), D_P1);
    }

        
  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22p01(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pxend,
    * restrict _mpx, * restrict _epx, * restrict _epy,
    * restrict _vx1, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, dt2, vx1, vx0, delta, etaxdt, etaydt;
  RARR *s;
    
  s = dom->_s;

  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
            
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
   
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

#pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epx,_epy,_vx1,_vy1,_vy0,vx1,vx0,delta,etaxdt,etaydt)
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
      px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
      py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
      mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
      vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
      vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
    }
#pragma omp barrier

    _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
    _vx1 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
    _vy1 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
    _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */

    for ( iy = tid; iy < ny; iy += tsz )
      {
	vx0 = _vx1[-1];
	etaydt = (*_epy) * dt2;
            
	for ( _pxend = _px + nx; _px < _pxend; )
	  {
	    vx1 = *_vx1++;
	    etaxdt = (*_epx++) * dt2;
                
	    delta = ((vx0 - vx1) * lax + ((*_vy0++) - (*_vy1++)) * lay) * (*_mpx++);
            
	    (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
	    _px++;
                
	    (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
	    _py++;
                
	    vx0 = vx1;
	  }
            
	_px += px_a; _py += py_a; _mpx += mpx_a; _vx1 += vx_a; _vy0 += vy_a; _vy1 += vy_a; _epx -= nx; _epy += tsz;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22p(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pxend,
    * restrict _mpx, * restrict _vx1, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, vx1, vx0, delta;
  RARR *s;

  s = dom->_s;
    
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
    
#pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_vx1,_vy1,_vy0,vx1,vx0,delta)
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
      px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
      py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
      mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
      vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
      vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
    }
#pragma omp barrier
        
    _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
    _vx1 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
    _vy1 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        
    for ( iy = tid; iy < ny; iy += tsz )
      {
	vx0 = _vx1[-1];
	    
	for ( _pxend = _px + nx; _px < _pxend; )
	  {
	    vx1 = *_vx1++;
		  
	    delta = ((vx0 - vx1) * lax + ((*_vy0++) - (*_vy1++)) * lay) * (*_mpx++);
	    (*_px++) += delta;
	    (*_py++) += delta;
		  
	    vx0 = vx1;
	  }

	_px += px_a; _py += py_a; _mpx += mpx_a; _vx1 += vx_a; _vy0 += vy_a; _vy1 += vy_a;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22p0(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pxend, * restrict _mpx,
    * restrict _epx, * restrict _vx1, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, dt2, vx1, vx0, delta, etaxdt;
  RARR *s;
    
  s = dom->_s;
    
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epx,_vx1,_vy1,_vy0,vx1,vx0,delta,etaxdt)
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
      px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
      py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
      mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
      vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
      vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
    }
#pragma omp barrier
        
    _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
    _vx1 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
    _vy1 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        
    for ( iy = tid; iy < ny; iy += tsz )
      {
	vx0 = _vx1[-1];
            
	for ( _pxend = _px + nx; _px < _pxend; )
	  {
	    vx1 = *_vx1++;
	    etaxdt = (*_epx++) * dt2;
		  
	    delta = ((vx0 - vx1) * lax + ((*_vy0++) - (*_vy1++)) * lay) * (*_mpx++);
	    (*_py++) += delta;            
		  
	    (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
	    _px++;
		  
	    vx0 = vx1;
	  }

	_px += px_a; _py += py_a; _mpx += mpx_a; _vx1 += vx_a; _vy0 += vy_a; _vy1 += vy_a; _epx -= nx;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22p1(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pxend, * restrict _mpx,
    * restrict _epy, * restrict _vx1, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, dt2, vx1, vx0, delta, etaydt;
  RARR *s;
    
  s = dom->_s;
    
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epy,_vx1,_vy1,_vy0,vx1,vx0,delta,etaydt)
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
      px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
      py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
      mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
      vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
      vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
    }
#pragma omp barrier
        
    _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
    _vx1 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
    _vy1 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
    _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */
            
    for ( iy = tid; iy < ny; iy += tsz )
      {
	vx0 = _vx1[-1];
	etaydt = (*_epy) * dt2;
            
	for ( _pxend = _px + nx; _px < _pxend; )
	  {
	    vx1 = *_vx1++;
                
	    delta = ((vx0 - vx1) * lax + ((*_vy0++) - (*_vy1++)) * lay) * (*_mpx++);
	    (*_px++) += delta;            
                
	    (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
	    _py++;
                
	    vx0 = vx1;
	  }
            
	_px += px_a; _py += py_a; _mpx += mpx_a; _vx1 += vx_a; _vy0 += vy_a; _vy1 += vy_a; _epy += tsz;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22v0(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, vx_a, mvx_a, px_a, iy, tsz, tid;
  register ireal * restrict _vx, * restrict _vxend, * restrict _mvx,
    * restrict _evx, * restrict _px1;
  register ireal lax, dt2, px1, px0, etaxdt;
  RARR *s;
    
  s = dom->_s;
    
  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_V0]._dims[0].gs;
  gys = s[D_V0]._dims[1].gs;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,_vx,_vxend,_mvx,_evx,_px1,px1,px0,etaxdt)
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
      vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
      mvx_a = tsz * s[D_MV0]._dims[0].n0 - nx;
      px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
    }
#pragma omp barrier
        
    _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
    _mvx = s[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs + tid) * s[D_MV0]._dims[0].n0;
    _px1 = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0 + 1;
    _evx = s[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
        
    for ( iy = tid; iy < ny; iy += tsz )
      {
	px0 = _px1[-1];
            
	for ( _vxend = _vx + nx; _vx < _vxend; )
	  {
	    px1 = *_px1++;
	    etaxdt = (*_evx++) * dt2;
                
	    (*_vx) = ((*_vx) * (1.0 - etaxdt) + (px0 - px1) * lax * (*_mvx++)) / (1.0 + etaxdt);
	    _vx++;
                
	    px0 = px1;
	  }
            
	_vx += vx_a; _mvx += mvx_a; _px1 += px_a; _evx -= nx;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_22v1(RDOM *dom, void *pars)
{
  int nx, ny, gxs, gys, vy_a, mvy_a, py_a, iy, tsz, tid;
  register ireal * restrict _vy, * restrict _vyend, * restrict _mvy,
    * restrict _evy, * restrict _py1, * restrict _py0;
  register ireal lay, dt2, etaydt;
  RARR *s;
    
  s = dom->_s;
    
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_V1]._dims[0].gs;
  gys = s[D_V1]._dims[1].gs;
    
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,_vy,_vyend,_mvy,_evy,_py1,_py0,etaydt)
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
      vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
      mvy_a = tsz * s[D_MV1]._dims[0].n0 - nx;
      py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
    }
#pragma omp barrier
        
    _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _mvy = s[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs + tid) * s[D_MV1]._dims[0].n0;
    _py0 = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _py1 = _py0 + s[D_P1]._dims[0].n0;
    _evy = s[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs + tid);        /* 1D */
        
    for ( iy = tid; iy < ny; iy += tsz )
      {            
	etaydt = (*_evy) * dt2;
            
	for ( _vyend = _vy + nx; _vy < _vyend; )
	  {
	    (*_vy) = ((*_vy) * (1.0 - etaydt) + ((*_py0++) - (*_py1++)) * lay * (*_mvy++)) / (1.0 + etaydt);
	    _vy++;
	  }
            
	_vy += vy_a; _mvy += mvy_a; _py0 += py_a; _py1 += py_a; _evy += tsz;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/
