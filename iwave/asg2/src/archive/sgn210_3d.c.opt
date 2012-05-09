/* 
   sgn210_3d.c
   Igor Terentyev.
   ********************************************************************************
   Implementation of 2-10 scheme in 3D.

   Optimized by Performancejones, summer 2011
*/
/*============================================================================*/

#include "utils.h"
#include "sgn.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/
/*
  Constants
*/
#define C1 ( -19845.0e0/16384.0e0  )
#define C2 (    735.0e0/8192.0e0   )
#define C3 (   -567.0e0/40960.0e0  )
#define C4 (    405.0e0/229376.0e0 )
#define C5 (    -35.0e0/294912.0e0 )
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p012(RDOM *, void *);
int sgn_ts3d_210p0j(RDOM *, void *, int);
int sgn_ts3d_210p12(RDOM *, void *);
int sgn_ts3d_210p0(RDOM *, void *);
int sgn_ts3d_210pj(RDOM *, void *, int);
int sgn_ts3d_210p(RDOM *, void *);
int sgn_ts3d_210v0(RDOM *, void *);
int sgn_ts3d_210vj(RDOM *, void *, int);

int sgn_ts3d_210(RDOM *dom, int iarr, void *pars)
{
  int e0, e1, e2;               /* short-hands */
  if ( iarr == D_P0 )           /* covers D_P0, D_P1 in one function */
    {
      /*
        if ( ((SGN_TS_PARS*)pars)->psingle ) return sg_ts3d_210p(dom, &(((SGN_TS_PARS*)pars)->sgpars));
        else
      */
      {
	e0 = ((SGN_TS_PARS*)pars)->eflags[0];
	e1 = ((SGN_TS_PARS*)pars)->eflags[1];
	e2 = ((SGN_TS_PARS*)pars)->eflags[2];
	if ( e0 && e1 && e2 ) return sgn_ts3d_210p012(dom, pars);
	else if ( e0 && e1 )  return sgn_ts3d_210p0j(dom, pars, 1);
	else if ( e0 && e2 )  return sgn_ts3d_210p0j(dom, pars, 2);
	else if ( e1 && e2 )  return sgn_ts3d_210p12(dom, pars);
	else if ( e0 )        return sgn_ts3d_210p0(dom, pars);
	else if ( e1 )        return sgn_ts3d_210pj(dom, pars, 1);
	else if ( e2 )        return sgn_ts3d_210pj(dom, pars, 2);
	else                  return sgn_ts3d_210p(dom, pars);
      }
    }
  if ( iarr == D_V0 )
    {
      /*
        if ( ((SGN_TS_PARS*)pars)->eflags[0] ) 
      */
      return sgn_ts3d_210v0(dom, pars);
      /*  else
	  return sg_ts3d_210v0(dom, &(((SGN_TS_PARS*)pars)->sgpars)); */
      
    }
  if ( iarr == D_V1 )
    {
      /*
        if ( ((SGN_TS_PARS*)pars)->eflags[1] ) 
      */
      return sgn_ts3d_210vj(dom, pars, 1);
      /*
        else
	return sg_ts3d_210vj(dom, &(((SGN_TS_PARS*)pars)->sgpars), 1, D_P1); */
    }
  if ( iarr == D_V2 )
    {
      /*
        if ( ((SGN_TS_PARS*)pars)->eflags[2] ) 
      */
      return sgn_ts3d_210vj(dom, pars, 2);
      /*
        else 
	return sg_ts3d_210vj(dom, &(((SGN_TS_PARS*)pars)->sgpars), 2, D_P2); */
    }
  // in all other cases, no-op
  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p012(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, ix, iy, iz, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _epx, * restrict _vx9,
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  ireal * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, dt2, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta, etaxdt, etaydt, etazdt;
  RARR *s;
  INFODIM *dims;
    
  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_P1 ]._dims[0].n0 - nx;
  pz_a  = s[D_P2 ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_epy,_epz,_vx9,	\
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,etaxdt,etaydt,etazdt,dims)
  {
    ireal C1x, C2x, C3x, C4x, C5x;
    ireal C1y, C2y, C3y, C4y, C5y;
    ireal C1z, C2z, C3z, C4z, C5z;

    ireal ratio_y, recip_y;
    ireal ratio_z, recip_z;

    ireal *ratio_x, *recip_x;

    ratio_x = (ireal *) malloc(nx*sizeof(ireal));
    recip_x = (ireal *) malloc(nx*sizeof(ireal));

#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif
        
#pragma omp single
    {            
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_P1 ]._dims[1].n0 - ny ) * s[D_P1 ]._dims[0].n0;
      pz_aa  = (tsz * s[D_P2 ]._dims[1].n0 - ny ) * s[D_P2 ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P1 ]._dims;
    _py  = s[D_P1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P2 ]._dims;
    _pz  = s[D_P2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;

    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
    _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);              /* 1D */
    _epz = s[D_EP2]._s + (tid - s[D_EP2]._dims[0].gs);              /* 1D */
        
    tid -= gzs; /* adjust back after calculation of pointers */

    C1x = C1 * lax;
    C2x = C2 * lax;
    C3x = C3 * lax;
    C4x = C4 * lax;
    C5x = C5 * lax;

    C1y = C1 * lay;
    C2y = C2 * lay;
    C3y = C3 * lay;
    C4y = C4 * lay;
    C5y = C5 * lay;

    C1z = C1 * laz;
    C2z = C2 * laz;
    C3z = C3 * laz;
    C4z = C4 * laz;
    C5z = C5 * laz;

    for ( ix = 0; ix < nx; ix++ )
      {
	etaxdt = _epx[ix] * dt2;

	ratio_x[ix] = (1.0 - etaxdt) / (1.0 + etaxdt);
	recip_x[ix] = 1.0 / (1.0 + etaxdt);
      }
        
    for ( iz = tid; iz < nz; iz += tsz )
      {
	etazdt = (*_epz) * dt2;

	ratio_z = (1.0 - etazdt) / (1.0 + etazdt);
	recip_z = 1.0 / (1.0 + etazdt);
            
	for ( iy = 0; iy < ny; ++iy )
	  {
	    etaydt = (*_epy++) * dt2;

	    ratio_y = (1.0 - etaydt) / (1.0 + etaydt);
	    recip_y = 1.0 / (1.0 + etaydt);
                
#pragma ivdep
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
		vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
		vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];

		vx9 = *_vx9++;

		delta = ((vx9 - vx0) * C5x + (vx8 - vx1) * C4x + (vx7 - vx2) * C3x + (vx6 - vx3) * C2x + (vx5 - vx4) * C1x +
			 ((*_vy9++) - (*_vy0++)) * C5y + ((*_vy8++) - (*_vy1++)) * C4y + ((*_vy7++) - (*_vy2++)) * C3y +
			 ((*_vy6++) - (*_vy3++)) * C2y + ((*_vy5++) - (*_vy4++)) * C1y +
			 ((*_vz9++) - (*_vz0++)) * C5z + ((*_vz8++) - (*_vz1++)) * C4z + ((*_vz7++) - (*_vz2++)) * C3z +
			 ((*_vz6++) - (*_vz3++)) * C2z + ((*_vz5++) - (*_vz4++)) * C1z) * (*_mpx++);
                    
		(*_px) = (*_px) * (*ratio_x++) + delta * (*recip_x++);
		_px++;
                    
		(*_py) = (*_py) * ratio_y + delta * recip_y;
		_py++;
                    
		(*_pz) = (*_pz) * ratio_z + delta * recip_z;
		_pz++;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a; 
	    ratio_x -= nx;
	    recip_x -= nx;
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
	_epy -= ny; _epz += tsz;
      }

    free(ratio_x);
    free(recip_x);

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _vx9,
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  register ireal lax, lay, laz, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta;
  RARR *s;
  INFODIM *dims;
    
  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_P1 ]._dims[0].n0 - nx;
  pz_a  = s[D_P2 ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_vx9, \
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,dims)
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
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_P1 ]._dims[1].n0 - ny ) * s[D_P1 ]._dims[0].n0;
      pz_aa  = (tsz * s[D_P2 ]._dims[1].n0 - ny ) * s[D_P2 ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P1 ]._dims;
    _py  = s[D_P1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P2 ]._dims;
    _pz  = s[D_P2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;
        
    tid -= gzs; /* adjust back after calculation of pointers */
        
    for ( iz = tid; iz < nz; iz += tsz )
      {
	for ( iy = 0; iy < ny; ++iy )
	  {
	    vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
	    vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
	    vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];
                
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx9 = *_vx9++;
                    
		delta = (( (vx9 - vx0) * C5 + (vx8 - vx1) * C4 + (vx7 - vx2) * C3 + (vx6 - vx3) * C2 + (vx5 - vx4) * C1) * lax + 
			 ( ((*_vy9++) - (*_vy0++)) * C5 + ((*_vy8++) - (*_vy1++)) * C4 + ((*_vy7++) - (*_vy2++)) * C3 + 
			   ((*_vy6++) - (*_vy3++)) * C2 + ((*_vy5++) - (*_vy4++)) * C1) * lay +
			 ( ((*_vz9++) - (*_vz0++)) * C5 + ((*_vz8++) - (*_vz1++)) * C4 + ((*_vz7++) - (*_vz2++)) * C3 + 
			   ((*_vz6++) - (*_vz3++)) * C2 + ((*_vz5++) - (*_vz4++)) * C1) * laz) * (*_mpx++);
                    
		(*_px++) += delta;
		(*_py++) += delta;
		(*_pz++) += delta;
                    
		vx0 = vx1; vx1 = vx2; vx2 = vx3; vx3 = vx4; vx4 = vx5; vx5 = vx6; vx6 = vx7; vx7 = vx8; vx8 = vx9;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a; 
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p0j(RDOM *dom, void *pars, int ind)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a, ep_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, ep_aa, iy, iz, tsz, tid, D_ep, D_py, D_pz;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _epx, * restrict _vx9,
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  ireal * restrict _ep;
  register ireal lax, lay, laz, dt2, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta, etaxdt, etadt;
  RARR *s;
  INFODIM *dims;
    
  D_ep = D_EP[ind];
  D_pz = D_P[ind];         /* pz always plays a role of pj */
  D_py = D_P[3 - ind];
    
  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_py ]._dims[0].n0 - nx;
  pz_a  = s[D_pz ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
  ep_a  = ((ind == 1) ? 1 : 0);
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_ep,_vx9, \
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,etaxdt,etadt,dims)
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
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_py ]._dims[1].n0 - ny ) * s[D_py ]._dims[0].n0;
      pz_aa  = (tsz * s[D_pz ]._dims[1].n0 - ny ) * s[D_pz ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
      ep_aa  = ((ind == 1) ? -ny : tsz);
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_py ]._dims;
    _py  = s[D_py ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_pz ]._dims;
    _pz  = s[D_pz ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;
    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
    _ep  = s[D_ep ]._s + ((ind == 1) ? gys : tid) - s[D_ep]._dims[0].gs ;
        
    tid -= gzs; /* adjust back after calculation of pointers */        
        
    for ( iz = tid; iz < nz; iz += tsz )
      {
	for ( iy = 0; iy < ny; ++iy )
	  {
	    vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
	    vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
	    vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];
	    etadt = (*_ep) * dt2;
                
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx9 = *_vx9++;
		etaxdt = (*_epx++) * dt2;
                    
		delta = (( (vx9 - vx0) * C5 + (vx8 - vx1) * C4 + (vx7 - vx2) * C3 + (vx6 - vx3) * C2 + (vx5 - vx4) * C1) * lax + 
			 ( ((*_vy9++) - (*_vy0++)) * C5 + ((*_vy8++) - (*_vy1++)) * C4 + ((*_vy7++) - (*_vy2++)) * C3 + 
			   ((*_vy6++) - (*_vy3++)) * C2 + ((*_vy5++) - (*_vy4++)) * C1) * lay +
			 ( ((*_vz9++) - (*_vz0++)) * C5 + ((*_vz8++) - (*_vz1++)) * C4 + ((*_vz7++) - (*_vz2++)) * C3 + 
			   ((*_vz6++) - (*_vz3++)) * C2 + ((*_vz5++) - (*_vz4++)) * C1) * laz) * (*_mpx++);
		(*_py++) += delta;
                    
		(*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
		_px++;
                    
		(*_pz) = ((*_pz) * (1.0 - etadt ) + delta) / (1.0 + etadt );
		_pz++;
                    
		vx0 = vx1; vx1 = vx2; vx2 = vx3; vx3 = vx4; vx4 = vx5; vx5 = vx6; vx6 = vx7; vx7 = vx8; vx8 = vx9;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a; 
	    _epx -= nx; _ep += ep_a;
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
	_ep += ep_aa;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p12(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _vx9,
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  ireal * restrict _epy, * restrict _epz;
  register ireal lax, lay, laz, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta, etaydt, etazdt;
  ireal dt2;
  RARR *s;
  INFODIM *dims;
    
  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_P1 ]._dims[0].n0 - nx;
  pz_a  = s[D_P2 ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epy,_epz,_vx9, \
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,etaydt,etazdt,dims)
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
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_P1 ]._dims[1].n0 - ny ) * s[D_P1 ]._dims[0].n0;
      pz_aa  = (tsz * s[D_P2 ]._dims[1].n0 - ny ) * s[D_P2 ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P1 ]._dims;
    _py  = s[D_P1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P2 ]._dims;
    _pz  = s[D_P2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;
    _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);              /* 1D */
    _epz = s[D_EP2]._s + (tid - s[D_EP2]._dims[0].gs);              /* 1D */
        
    tid -= gzs; /* adjust back after calculation of pointers */
        
    for ( iz = tid; iz < nz; iz += tsz )
      {
	etazdt = (*_epz) * dt2;
            
	for ( iy = 0; iy < ny; ++iy )
	  {
	    vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
	    vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
	    vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];
	    etaydt = (*_epy++) * dt2;
                
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx9 = *_vx9++;
                    
		delta = (( (vx9 - vx0) * C5 + (vx8 - vx1) * C4 + (vx7 - vx2) * C3 + (vx6 - vx3) * C2 + (vx5 - vx4) * C1) * lax + 
			 ( ((*_vy9++) - (*_vy0++)) * C5 + ((*_vy8++) - (*_vy1++)) * C4 + ((*_vy7++) - (*_vy2++)) * C3 + 
			   ((*_vy6++) - (*_vy3++)) * C2 + ((*_vy5++) - (*_vy4++)) * C1) * lay +
			 ( ((*_vz9++) - (*_vz0++)) * C5 + ((*_vz8++) - (*_vz1++)) * C4 + ((*_vz7++) - (*_vz2++)) * C3 + 
			   ((*_vz6++) - (*_vz3++)) * C2 + ((*_vz5++) - (*_vz4++)) * C1) * laz) * (*_mpx++);
		(*_px++) += delta;
                    
		(*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
		_py++;
                    
		(*_pz) = ((*_pz) * (1.0 - etazdt) + delta) / (1.0 + etazdt);
		_pz++;
                    
		vx0 = vx1; vx1 = vx2; vx2 = vx3; vx3 = vx4; vx4 = vx5; vx5 = vx6; vx6 = vx7; vx7 = vx8; vx8 = vx9;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a; 
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
	_epy -= ny; _epz += tsz;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210pj(RDOM *dom, void *pars, int ind)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a, ep_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, ep_aa, iy, iz, tsz, tid, D_ep, D_py, D_pz;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _vx9,
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  ireal * restrict _ep;
  register ireal lax, lay, laz, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta, etadt;
  ireal dt2;
  RARR *s;
  INFODIM *dims;
    
  D_ep = D_EP[ind];
  D_pz = D_P[ind];         /* pz always plays a role of pj */
  D_py = D_P[3 - ind];
    
  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_py ]._dims[0].n0 - nx;
  pz_a  = s[D_pz ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
  ep_a  = ((ind == 1) ? 1 : 0);
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_ep,_vx9, \
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,etadt,dims)
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
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_py ]._dims[1].n0 - ny ) * s[D_py ]._dims[0].n0;
      pz_aa  = (tsz * s[D_pz ]._dims[1].n0 - ny ) * s[D_pz ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
      ep_aa  = ((ind == 1) ? -ny : tsz);
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_py ]._dims;
    _py  = s[D_py ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_pz ]._dims;
    _pz  = s[D_pz ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;
    _ep  = s[D_ep ]._s + ((ind == 1) ? gys : tid) - s[D_ep]._dims[0].gs ;
        
    tid -= gzs; /* adjust back after calculation of pointers */
        
    for ( iz = tid; iz < nz; iz += tsz )
      {
	for ( iy = 0; iy < ny; ++iy )
	  {
	    vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
	    vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
	    vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];
	    etadt = (*_ep) * dt2;
                
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx9 = *_vx9++;
                    
		delta = (( (vx9 - vx0) * C5 + (vx8 - vx1) * C4 + (vx7 - vx2) * C3 + (vx6 - vx3) * C2 + (vx5 - vx4) * C1) * lax + 
			 ( ((*_vy9++) - (*_vy0++)) * C5 + ((*_vy8++) - (*_vy1++)) * C4 + ((*_vy7++) - (*_vy2++)) * C3 + 
			   ((*_vy6++) - (*_vy3++)) * C2 + ((*_vy5++) - (*_vy4++)) * C1) * lay +
			 ( ((*_vz9++) - (*_vz0++)) * C5 + ((*_vz8++) - (*_vz1++)) * C4 + ((*_vz7++) - (*_vz2++)) * C3 + 
			   ((*_vz6++) - (*_vz3++)) * C2 + ((*_vz5++) - (*_vz4++)) * C1) * laz) * (*_mpx++);
		(*_px++) += delta;
		(*_py++) += delta;
                    
		(*_pz) = ((*_pz) * (1.0 - etadt) + delta) / (1.0 + etadt);
		_pz++;
                    
		vx0 = vx1; vx1 = vx2; vx2 = vx3; vx3 = vx4; vx4 = vx5; vx5 = vx6; vx6 = vx7; vx7 = vx8; vx8 = vx9;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a;
	    _ep += ep_a;
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
	_ep += ep_aa;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210p0(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
    px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
  register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
    * restrict _mpx, * restrict _epx, * restrict _vx9, 
    * restrict _vy9, * restrict _vy8, * restrict _vy7, 
    * restrict _vy6, * restrict _vy5, * restrict _vy4,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0,
    * restrict _vz9, * restrict _vz8, * restrict _vz7, 
    * restrict _vz6, * restrict _vz5, * restrict _vz4,
    * restrict _vz3, * restrict _vz2, * restrict _vz1, * restrict _vz0;
  register ireal lax, lay, laz, dt2, vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0, delta, etaxdt;
  RARR *s;
  INFODIM *dims;

  s = dom->_s;
    
  dims = s[D_P0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
  py_a  = s[D_P1 ]._dims[0].n0 - nx;
  pz_a  = s[D_P2 ]._dims[0].n0 - nx;
  mpx_a = s[D_MP0]._dims[0].n0 - nx;
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  vy_a  = s[D_V1 ]._dims[0].n0 - nx;
  vz_a  = s[D_V2 ]._dims[0].n0 - nx;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  lay = ((SGN_TS_PARS*)pars)->lam[1];
  laz = ((SGN_TS_PARS*)pars)->lam[2];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_vx9, \
			     _vy9,_vy8,_vy7,_vy6,_vy5,_vy4,_vy3,_vy2,_vy1,_vy0,	\
			     _vz9,_vz8,_vz7,_vz6,_vz5,_vz4,_vz3,_vz2,_vz1,_vz0,	\
			     vx9, vx8, vx7, vx6, vx5, vx4, vx3, vx2, vx1, vx0,delta,etaxdt,dims)
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
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
      py_aa  = (tsz * s[D_P1 ]._dims[1].n0 - ny ) * s[D_P1 ]._dims[0].n0;
      pz_aa  = (tsz * s[D_P2 ]._dims[1].n0 - ny ) * s[D_P2 ]._dims[0].n0;
      mpx_aa = (tsz * s[D_MP0]._dims[1].n0 - ny ) * s[D_MP0]._dims[0].n0;
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      vy_aa  = (tsz * s[D_V1 ]._dims[1].n0 - ny ) * s[D_V1 ]._dims[0].n0;
      vz_aa  = (tsz * s[D_V2 ]._dims[1].n0 - ny ) * s[D_V2 ]._dims[0].n0;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_P0 ]._dims;
    _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P1 ]._dims;
    _py  = s[D_P1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_P2 ]._dims;
    _pz  = s[D_P2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MP0]._dims;
    _mpx = s[D_MP0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_V0 ]._dims;
    _vx9 = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 4;
    dims = s[D_V1 ]._dims;
    _vy5 = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vy6 = _vy5 + dims[0].n0;
    _vy7 = _vy6 + dims[0].n0;
    _vy8 = _vy7 + dims[0].n0;
    _vy9 = _vy8 + dims[0].n0;
    _vy4 = _vy5 - dims[0].n0;
    _vy3 = _vy4 - dims[0].n0;
    _vy2 = _vy3 - dims[0].n0;
    _vy1 = _vy2 - dims[0].n0;
    _vy0 = _vy1 - dims[0].n0;
    dims = s[D_V2]._dims;
    _vz5 = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _vz6 = _vz5 + dims[1].n0 * dims[0].n0;
    _vz7 = _vz6 + dims[1].n0 * dims[0].n0;
    _vz8 = _vz7 + dims[1].n0 * dims[0].n0;
    _vz9 = _vz8 + dims[1].n0 * dims[0].n0;
    _vz4 = _vz5 - dims[1].n0 * dims[0].n0;
    _vz3 = _vz4 - dims[1].n0 * dims[0].n0;
    _vz2 = _vz3 - dims[1].n0 * dims[0].n0;
    _vz1 = _vz2 - dims[1].n0 * dims[0].n0;
    _vz0 = _vz1 - dims[1].n0 * dims[0].n0;
    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        
    tid -= gzs; /* adjust back after calculation of pointers */

    for ( iz = tid; iz < nz; iz += tsz )
      {
	for ( iy = 0; iy < ny; ++iy )
	  {
	    vx8 = _vx9[-1]; vx7 = _vx9[-2]; vx6 = _vx9[-3];
	    vx5 = _vx9[-4]; vx4 = _vx9[-5]; vx3 = _vx9[-6];
	    vx2 = _vx9[-7]; vx1 = _vx9[-8]; vx0 = _vx9[-9];
                
	    for ( _pxend = _px + nx; _px < _pxend; )
	      {
		vx9 = *_vx9++;
		etaxdt = (*_epx++) * dt2;
                    
		delta = (( (vx9 - vx0) * C5 + (vx8 - vx1) * C4 + (vx7 - vx2) * C3 + (vx6 - vx3) * C2 + (vx5 - vx4) * C1) * lax + 
			 ( ((*_vy9++) - (*_vy0++)) * C5 + ((*_vy8++) - (*_vy1++)) * C4 + ((*_vy7++) - (*_vy2++)) * C3 + 
			   ((*_vy6++) - (*_vy3++)) * C2 + ((*_vy5++) - (*_vy4++)) * C1) * lay +
			 ( ((*_vz9++) - (*_vz0++)) * C5 + ((*_vz8++) - (*_vz1++)) * C4 + ((*_vz7++) - (*_vz2++)) * C3 + 
			   ((*_vz6++) - (*_vz3++)) * C2 + ((*_vz5++) - (*_vz4++)) * C1) * laz) * (*_mpx++);
		(*_py++) += delta;
		(*_pz++) += delta;
                    
		(*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
		_px++;
                    
		vx0 = vx1; vx1 = vx2; vx2 = vx3; vx3 = vx4; vx4 = vx5; vx5 = vx6; vx6 = vx7; vx7 = vx8; vx8 = vx9;
	      }
                
	    _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx9 += vx_a;
	    _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _vy4 += vy_a; 
	    _vy5 += vy_a; _vy6 += vy_a; _vy7 += vy_a; _vy8 += vy_a; _vy9 += vy_a; 
	    _vz0 += vz_a; _vz1 += vz_a; _vz2 += vz_a; _vz3 += vz_a; _vz4 += vz_a; 
	    _vz5 += vz_a; _vz6 += vz_a; _vz7 += vz_a; _vz8 += vz_a; _vz9 += vz_a; 
	    _epx -= nx;
	  }
            
	_px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx9 += vx_aa;
	_vy0 += vy_aa; _vy1 += vy_aa; _vy2 += vy_aa; _vy3 += vy_aa; _vy4 += vy_aa; 
	_vy5 += vy_aa; _vy6 += vy_aa; _vy7 += vy_aa; _vy8 += vy_aa; _vy9 += vy_aa; 
	_vz0 += vz_aa; _vz1 += vz_aa; _vz2 += vz_aa; _vz3 += vz_aa; _vz4 += vz_aa; 
	_vz5 += vz_aa; _vz6 += vz_aa; _vz7 += vz_aa; _vz8 += vz_aa; _vz9 += vz_aa; 
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

/* Pad buffer (Cx) rows for 32-byte alignment */
#define C_NCOL 8

int sgn_ts3d_210v0(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, vx_a, mvx_a, px_a, vx_aa, mvx_aa, px_aa,
    ix, iy, iz, tsz, tid;
  register ireal * restrict _vx, * restrict _mvx,
    * restrict _evx, * restrict _px9;
  register ireal lax, dt2, etaxdt;
  RARR *s;
  INFODIM *dims;
    
  s = dom->_s;
    
  dims = s[D_V0]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  vx_a  = s[D_V0 ]._dims[0].n0 - nx;
  mvx_a = s[D_MV0]._dims[0].n0 - nx;
  px_a  = s[D_P0 ]._dims[0].n0 - nx;
    
  lax = ((SGN_TS_PARS*)pars)->lam[0];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_vx,_mvx,_evx,_px9,etaxdt,dims)
  {
    ireal (*Cx)[C_NCOL];

    Cx = (ireal (*)[C_NCOL]) malloc(nx*(C_NCOL)*sizeof(ireal));

#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif
        
#pragma omp single
    {
      vx_aa  = (tsz * s[D_V0 ]._dims[1].n0 - ny ) * s[D_V0 ]._dims[0].n0;
      mvx_aa = (tsz * s[D_MV0]._dims[1].n0 - ny ) * s[D_MV0]._dims[0].n0;
      px_aa  = (tsz * s[D_P0 ]._dims[1].n0 - ny ) * s[D_P0 ]._dims[0].n0;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_V0 ]._dims;
    _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_MV0]._dims;
    _mvx = s[D_MV0]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;        
    dims = s[D_P0 ]._dims;
    _px9 = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 5;
    _evx = s[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
        
    tid -= gzs; /* adjust back after calculation of pointers */

#pragma ivdep
    for ( ix = 0; ix < nx; ix++ )
      {
	ireal ratio;

	etaxdt = _evx[ix] * dt2;

	ratio = (1.0 - etaxdt) / (1.0 + etaxdt);

	Cx[ix][0] = ratio;
	Cx[ix][1] = C1 * lax / (1.0 + etaxdt);
	Cx[ix][2] = C2 * lax / (1.0 + etaxdt);
	Cx[ix][3] = C3 * lax / (1.0 + etaxdt);
	Cx[ix][4] = C4 * lax / (1.0 + etaxdt);
	Cx[ix][5] = C5 * lax / (1.0 + etaxdt);
      }
        
    for ( iz = tid; iz < nz; iz += tsz )
      {

	for ( iy = 0; iy < ny; ++iy )
	  {

#pragma ivdep
	    for ( ix = 0; ix < nx; ix++ )
	      {
		ireal ratio = Cx[ix][0];

		_vx[ix] = (_vx[ix] * ratio + 
			   ( (_px9[ix  ] - _px9[ix-9]) * Cx[ix][5] +
			     (_px9[ix-1] - _px9[ix-8]) * Cx[ix][4] +
			     (_px9[ix-2] - _px9[ix-7]) * Cx[ix][3] +
			     (_px9[ix-3] - _px9[ix-6]) * Cx[ix][2] +
			     (_px9[ix-4] - _px9[ix-5]) * Cx[ix][1]   ) * _mvx[ix]);
	      }
                
	    _vx += nx + vx_a; _mvx += nx + mvx_a; _px9 += nx + px_a;
	  }
            
	_vx += vx_aa; _mvx += mvx_aa; _px9 += px_aa;
      }

    free(Cx);

  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_210vj_1(RDOM *dom, void *pars);
int sgn_ts3d_210vj_0(RDOM *dom, void *pars, int ind);

int sgn_ts3d_210vj(RDOM *dom, void *pars, int ind)
{
  if ( ind == 1 )
    {
      return sgn_ts3d_210vj_1(dom, pars);
    }
  else
    {
      return sgn_ts3d_210vj_0(dom, pars, ind);
    }
}

int sgn_ts3d_210vj_1(RDOM *dom, void *pars)
{
  int nx, ny, nz, gxs, gys, gzs, v_a, mv_a, p_a, ev_a, v_aa, mv_aa, p_aa,
    ev_aa, iy, iz, tsz, tid, D_v, D_mv, D_ev, D_p, p_shift;
  register ireal * restrict _v, * restrict _p4, * restrict _mv;
  ireal * restrict _ev;
  register ireal la, etadt;
  ireal dt2;
  RARR *s;
  INFODIM *dims;
    
  D_v  = D_V[1];
  D_mv = D_MV[1];
  D_ev = D_EV[1];
  D_p  = D_P[1];

  s = dom->_s;
    
  dims = s[D_v]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  v_a  = s[D_v ]._dims[0].n0 - nx;
  mv_a = s[D_mv]._dims[0].n0 - nx;
  p_a  = s[D_p ]._dims[0].n0 - nx;
  ev_a = 1;
    
  la = ((SGN_TS_PARS*)pars)->lam[1];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

  p_shift = s[D_p]._dims[0].n0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_v,_p4,_mv,_ev,etadt,dims)
  {
    ireal (*Cy)[C_NCOL];

    Cy = (ireal (*)[C_NCOL]) malloc(ny*(C_NCOL)*sizeof(ireal));

#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif
        
#pragma omp single
    {
      v_aa  = (tsz * s[D_v ]._dims[1].n0 - ny ) * s[D_v ]._dims[0].n0;
      mv_aa = (tsz * s[D_mv]._dims[1].n0 - ny ) * s[D_mv]._dims[0].n0;
      p_aa  = (tsz * s[D_p ]._dims[1].n0 - ny ) * s[D_p ]._dims[0].n0;
      ev_aa = -ny;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_v ]._dims;
    _v   = s[D_v ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_mv]._dims;
    _mv  = s[D_mv]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;        
    dims = s[D_p ]._dims;
    _p4  = s[D_p ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _ev  = s[D_ev]._s + gys - s[D_ev]._dims[0].gs ;

    tid -= gzs; /* adjust back after calculation of pointers */

    for ( iz = tid; iz < nz; iz += tsz )
      {
	int ix;

	for ( iy = 0; iy < ny; iy++ )
	  {
	    ireal ratio;
	    ireal C1x, C2x, C3x, C4x, C5x;

	    etadt = _ev[iy] * dt2;

	    ratio = (1.0 - etadt) / (1.0 + etadt);

	    C1x = C1 * la / (1.0 + etadt);
	    C2x = C2 * la / (1.0 + etadt);
	    C3x = C3 * la / (1.0 + etadt);
	    C4x = C4 * la / (1.0 + etadt);
	    C5x = C5 * la / (1.0 + etadt);

	    Cy[iy][0] = ratio;
	    Cy[iy][1] = C1x;
	    Cy[iy][2] = C2x;
	    Cy[iy][3] = C3x;
	    Cy[iy][4] = C4x;
	    Cy[iy][5] = C5x;
	  }

	for ( iy = 0; iy < ny; iy++ )
	  {
#pragma ivdep
	    for ( ix = 0; ix < nx; ix++ )
	      {
		_v[ix] = _v[ix] * Cy[iy][0] +
		  ((_p4[ix + 5*p_shift] - _p4[ix - 4*p_shift]) * Cy[iy][5] +
		   (_p4[ix + 4*p_shift] - _p4[ix - 3*p_shift]) * Cy[iy][4] +
		   (_p4[ix + 3*p_shift] - _p4[ix - 2*p_shift]) * Cy[iy][3] +
		   (_p4[ix + 2*p_shift] - _p4[ix -   p_shift]) * Cy[iy][2] +
		   (_p4[ix +   p_shift] - _p4[ix            ]) * Cy[iy][1]   ) * _mv[ix];
	      }

	    _v += nx + v_a; _mv += nx + mv_a;
	    _p4 += p_shift;
	  }

	_v += v_aa; _mv += mv_aa;
	_p4 += p_aa; _ev += ny + ev_aa;
      }

    free(Cy);

  } /* omp parallel */

  return 0;
}

int sgn_ts3d_210vj_0(RDOM *dom, void *pars, int ind)
{
  int nx, ny, nz, gxs, gys, gzs, v_a, mv_a, p_a, ev_a, v_aa, mv_aa, p_aa,
    ev_aa, iy, iz, tsz, tid, D_v, D_mv, D_ev, D_p, p_shift;
  register ireal * restrict _v, * restrict _p4, * restrict _mv;
  ireal * restrict _ev;
  register ireal la, etadt;
  ireal dt2;
  RARR *s;
  INFODIM *dims;
    
  D_v  = D_V[ind];
  D_mv = D_MV[ind];
  D_ev = D_EV[ind];
  D_p  = D_P[ind];

  s = dom->_s;
    
  dims = s[D_v]._dims;
    
  nx = dims[0].n;
  ny = dims[1].n;
  nz = dims[2].n;
  if ( nx * ny * nz == 0 ) return 0;
    
  gxs = dims[0].gs;
  gys = dims[1].gs;
  gzs = dims[2].gs;
    
  v_a  = s[D_v ]._dims[0].n0 - nx;
  mv_a = s[D_mv]._dims[0].n0 - nx;
  p_a  = s[D_p ]._dims[0].n0 - nx;
  ev_a = 0;
    
  la = ((SGN_TS_PARS*)pars)->lam[ind];
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

  p_shift = s[D_p]._dims[0].n0 * s[D_p]._dims[1].n0;
    
#pragma omp parallel private(tsz,tid,iy,iz,_v,_p4,_mv,_ev,etadt,dims)
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
      v_aa  = (tsz * s[D_v ]._dims[1].n0 - ny ) * s[D_v ]._dims[0].n0;
      mv_aa = (tsz * s[D_mv]._dims[1].n0 - ny ) * s[D_mv]._dims[0].n0;
      p_aa  = (tsz * s[D_p ]._dims[1].n0 - ny ) * s[D_p ]._dims[0].n0;
      ev_aa = tsz;
    }
#pragma omp barrier
        
    tid += gzs; /* adjust for shorter calculation of pointers */
        
    dims = s[D_v ]._dims;
    _v   = s[D_v ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    dims = s[D_mv]._dims;
    _mv  = s[D_mv]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;        
    dims = s[D_p ]._dims;
    _p4  = s[D_p ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
    _ev  = s[D_ev]._s + tid - s[D_ev]._dims[0].gs ;

    tid -= gzs; /* adjust back after calculation of pointers */

    for ( iz = tid; iz < nz; iz += tsz )
      {
	ireal ratio;
	ireal C1x, C2x, C3x, C4x, C5x;

	int ix;

	etadt = (*_ev) * dt2;

	ratio = (1.0 - etadt) / (1.0 + etadt);

	C1x = C1 * la / (1.0 + etadt);
	C2x = C2 * la / (1.0 + etadt);
	C3x = C3 * la / (1.0 + etadt);
	C4x = C4 * la / (1.0 + etadt);
	C5x = C5 * la / (1.0 + etadt);

	for ( iy = 0; iy < ny; ++iy )
	  {
#pragma ivdep
	    for ( ix = 0; ix < nx; ix++ )
	      {
		_v[ix] = _v[ix] * ratio +
		  ((_p4[ix + 5*p_shift] - _p4[ix - 4*p_shift]) * C5x +
		   (_p4[ix + 4*p_shift] - _p4[ix - 3*p_shift]) * C4x +
		   (_p4[ix + 3*p_shift] - _p4[ix - 2*p_shift]) * C3x +
		   (_p4[ix + 2*p_shift] - _p4[ix -   p_shift]) * C2x +
		   (_p4[ix +   p_shift] - _p4[ix            ]) * C1x   ) * _mv[ix];
	      }

	    _v += nx + v_a; _mv += nx + mv_a;
	    _p4 += nx + p_a;
	  }

	_v += v_aa; _mv += mv_aa;
	_p4 += p_aa; _ev += ev_aa;
      }
  } /* omp parallel */

  return 0;
}
/*----------------------------------------------------------------------------*/
