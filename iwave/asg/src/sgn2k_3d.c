/* 
sgn2k_3d.c
Igor Terentyev.
********************************************************************************
Implementation of 2-2k scheme in 3D.
*/
/*============================================================================*/

#include "utils.h"
#include "sgn.h"
#include "sgcoeffs.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp012(RDOM *, void *);
int sgn_ts3d_2Kp0j(RDOM *, void *, int);
int sgn_ts3d_2Kp12(RDOM *, void *);
int sgn_ts3d_2Kp0(RDOM *, void *);
int sgn_ts3d_2Kp012(RDOM *, void *);
int sgn_ts3d_2Kpj(RDOM *, void *, int);
int sgn_ts3d_2Kp(RDOM *, void *);
int sgn_ts3d_2Kv0(RDOM *, void *);
int sgn_ts3d_2Kvj(RDOM *, void *, int);

int sgn_ts3d_2K(RDOM *dom, int iarr, void *pars)
{
    int e0, e1, e2;               /* short-hands */
    if ( iarr == D_P0 )           /* covers D_P0, D_P1 in one function */
    {
      /*
        if ( ((SGN_TS_PARS*)pars)->psingle ) return sg_ts3d_2Kp(dom, &(((SGN_TS_PARS*)pars)->sgpars));
        else
      */
        {
            e0 = ((SGN_TS_PARS*)pars)->eflags[0];
            e1 = ((SGN_TS_PARS*)pars)->eflags[1];
            e2 = ((SGN_TS_PARS*)pars)->eflags[2];
            if ( e0 && e1 && e2 ) return sgn_ts3d_2Kp012(dom, pars);
            else if ( e0 && e1 )  return sgn_ts3d_2Kp0j(dom, pars, 1);
            else if ( e0 && e2 )  return sgn_ts3d_2Kp0j(dom, pars, 2);
            else if ( e1 && e2 )  return sgn_ts3d_2Kp12(dom, pars);
            else if ( e0 )        return sgn_ts3d_2Kp0(dom, pars);
            else if ( e1 )        return sgn_ts3d_2Kpj(dom, pars, 1);
            else if ( e2 )        return sgn_ts3d_2Kpj(dom, pars, 2);
            else                  return sgn_ts3d_2Kp(dom, pars);
        }
    }
    if ( iarr == D_V0 )
    {
      //if ( ((SGN_TS_PARS*)pars)->eflags[0] ) 
      return sgn_ts3d_2Kv0(dom, pars);
      //  else return sg_ts3d_2Kv0(dom, &(((SGN_TS_PARS*)pars)->sgpars));
    }
    if ( iarr == D_V1 )
    {
      //if ( ((SGN_TS_PARS*)pars)->eflags[1] ) 
      return sgn_ts3d_2Kvj(dom, pars, 1);
      //else return sg_ts3d_2Kvj(dom, &(((SGN_TS_PARS*)pars)->sgpars), 1, D_P1);
    }
    if ( iarr == D_V2 )
    {
      //if ( ((SGN_TS_PARS*)pars)->eflags[2] ) 
      return sgn_ts3d_2Kvj(dom, pars, 2);
      //  else return sg_ts3d_2Kvj(dom, &(((SGN_TS_PARS*)pars)->sgpars), 2, D_P2);
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp012(RDOM *dom, void *pars)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _epx, * restrict _vx, 
                  * restrict _vy,  * restrict _vz;
    ireal * restrict _epy, * restrict _epz;
    register ireal lax, lay, laz, dt2, delta, etaxdt, etaydt, etazdt, px, py, pz;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];
    
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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

   #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_epy,_epz,_vx,_vy,_vz,delta,etaxdt,etaydt,etazdt,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2]._dims;
        _vz  = s[D_V2]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);              /* 1D */
        _epz = s[D_EP2]._s + (tid - s[D_EP2]._dims[0].gs);              /* 1D */
        
        tid -= gzs; /* adjust back after calculation of pointers */
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            etazdt = (*_epz) * dt2;
            
            for ( iy = 0; iy < ny; ++iy )
            {
                etaydt = (*_epy++) * dt2;
                
                for ( _pxend = _px + nx; _px < _pxend; )
                {
                    etaxdt = (*_epx++) * dt2;
                    
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    
                    (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
                    ++_px;
                    
                    (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
                    ++_py;
                    
                    (*_pz) = ((*_pz) * (1.0 - etazdt) + delta) / (1.0 + etazdt);
                    ++_pz;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx += vx_a;
                _vy += vy_a; _vz += vz_a; _epx -= nx;
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx += vx_aa;
            _vy += vy_aa; _vz += vz_aa; _epy -= ny; _epz += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp(RDOM *dom, void *pars)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _vx, * restrict _vy, * restrict _vz;
    register ireal lax, lay, laz, delta, px, py, pz;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];

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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

    #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_vx,_vy,_vz,delta,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2 ]._dims;
        _vz  = s[D_V2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        
        tid -= gzs; /* adjust back after calculation of pointers */
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                for ( _pxend = _px + nx; _px < _pxend; )
                {
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    (*_px++) += delta;
                    (*_py++) += delta;
                    (*_pz++) += delta;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a;
                _vx += vx_a; _vy += vy_a; _vz += vz_a; 
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa;
            _vx += vx_aa; _vy += vy_aa; _vz += vz_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp0j(RDOM *dom, void *pars, int ind)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, py_a, px_a, pz_a, mpx_a, vx_a, vy_a, vz_a, ep_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, ep_aa, iy, iz, tsz, tid, D_ep, D_py, D_pz;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _epx, * restrict _vx, 
                  * restrict _vy, * restrict _vz;
    ireal * restrict _ep;
    register ireal lax, lay, laz, dt2, delta, etaxdt, etadt, px, py, pz;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];

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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

    #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_ep,_vx,_vy,_vz,delta,etaxdt,etadt,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2 ]._dims;
        _vz  = s[D_V2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        _ep  = s[D_ep ]._s + ((ind == 1) ? gys : tid) - s[D_ep]._dims[0].gs ;
        
        tid -= gzs; /* adjust back after calculation of pointers */        
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                etadt = (*_ep) * dt2;
                for ( _pxend = _px + nx; _px < _pxend; )
                {
                    etaxdt = (*_epx++) * dt2;
                    
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    (*_py++) += delta;
                    
                    (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
                    ++_px;
                    
                    (*_pz) = ((*_pz) * (1.0 - etadt ) + delta) / (1.0 + etadt );
                    ++_pz;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx += vx_a;
                _vy += vy_a; _vz += vz_a; _epx -= nx; _ep += ep_a;
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx += vx_aa;
            _vy += vy_aa; _vz += vz_aa; _ep += ep_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp12(RDOM *dom, void *pars)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _vx,
                  * restrict _vy, * restrict _vz;
    ireal * restrict _epy, * restrict _epz;
    register ireal lax, lay, laz, delta, etaydt, etazdt, px, py, pz;
    ireal dt2;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];
    
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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

    #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epy,_epz,_vx,_vy,_vz,delta,etaydt,etazdt,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2 ]._dims;
        _vz  = s[D_V2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);              /* 1D */
        _epz = s[D_EP2]._s + (tid - s[D_EP2]._dims[0].gs);              /* 1D */
        
        tid -= gzs; /* adjust back after calculation of pointers */
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            etazdt = (*_epz) * dt2;
            
            for ( iy = 0; iy < ny; ++iy )
            {
                etaydt = (*_epy++) * dt2;
                
                for ( _pxend = _px + nx; _px < _pxend; )
                {
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    (*_px++) += delta;
                    
                    (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
                    ++_py;
                    
                    (*_pz) = ((*_pz) * (1.0 - etazdt) + delta) / (1.0 + etazdt);
                    ++_pz;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a;
                _vx += vx_a; _vy += vy_a; _vz += vz_a;
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx += vx_aa;
            _vy += vy_aa; _vz += vz_aa; _epy -= ny; _epz += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kpj(RDOM *dom, void *pars, int ind)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a, ep_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, ep_aa, iy, iz, tsz, tid, D_ep, D_py, D_pz;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _vx, * restrict _vy, * restrict _vz;
    ireal * restrict _ep;
    register ireal lax, lay, laz, delta, etadt, px, py, pz;
    ireal dt2;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];
    
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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

    #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_ep,_vx,_vy,_vz,delta,etadt,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2 ]._dims;
        _vz  = s[D_V2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        _ep  = s[D_ep ]._s + ((ind == 1) ? gys : tid) - s[D_ep]._dims[0].gs ;
        
        tid -= gzs; /* adjust back after calculation of pointers */
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                etadt = (*_ep) * dt2;
                
                for ( _pxend = _px + nx; _px < _pxend; )
                {
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    (*_px++) += delta;
                    (*_py++) += delta;
                    
                    (*_pz) = ((*_pz) * (1.0 - etadt) + delta) / (1.0 + etadt);
                    ++_pz;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx += vx_a;
                _vy += vy_a; _vz += vz_a; _ep += ep_a;
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa; _vx += vx_aa;
            _vy += vy_aa; _vz += vz_aa; _ep += ep_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kp0(RDOM *dom, void *pars)
{
	int i, k, vys, vzs;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, px_a, py_a, pz_a, mpx_a, vx_a, vy_a, vz_a,
        px_aa, py_aa, pz_aa, mpx_aa, vx_aa, vy_aa, vz_aa, iy, iz, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pz, * restrict _pxend,
                  * restrict _mpx, * restrict _epx, * restrict _vx, 
                  * restrict _vy, * restrict _vz;
    register ireal lax, lay, laz, dt2, delta, etaxdt, px, py, pz;
    RARR *s;
    INFODIM *dims;

	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];
    
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
    
 	vys = s[D_V1]._dims[0].n0;
 	vzs = s[D_V2]._dims[1].n0 * s[D_V2]._dims[0].n0;

    #pragma omp parallel private(tsz,tid,iy,iz,_px,_py,_pz,_pxend,_mpx,_epx,_vx,_vy,_vz,delta,etaxdt,px,py,pz,i,dims)
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
        _vx  = s[D_V0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V1 ]._dims;
        _vy  = s[D_V1 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_V2 ]._dims;
        _vz  = s[D_V2 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        
        tid -= gzs; /* adjust back after calculation of pointers */

        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                for ( _pxend = _px + nx; _px < _pxend; )
                {
                    etaxdt = (*_epx++) * dt2;
                    
					px = py = pz = 0.0;
					for ( i = 0; i < k; ++i )
					{
						px += (_vx[i      ] - _vx[(-1-i)      ]) * coeffs[i];
						py += (_vy[i * vys] - _vy[(-1-i) * vys]) * coeffs[i];
						pz += (_vz[i * vzs] - _vz[(-1-i) * vzs]) * coeffs[i];
					}
					++_vx;
					++_vy;
					++_vz;
                    delta = (px * lax + py * lay + pz * laz) * (*_mpx++);
                    (*_py++) += delta;
                    (*_pz++) += delta;
                    
                    (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
                    ++_px;
                }
                
                _px += px_a; _py += py_a; _pz += pz_a; _mpx += mpx_a; _vx += vx_a;
                _vy += vy_a; _vz += vz_a; _epx -= nx;
            }
            
            _px += px_aa; _py += py_aa; _pz += pz_aa; _mpx += mpx_aa;
            _vx += vx_aa; _vy += vy_aa; _vz += vz_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kv0(RDOM *dom, void *pars)
{
	int i, k;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, vx_a, mvx_a, px_a, vx_aa, mvx_aa, px_aa,
        iy, iz, tsz, tid;
    register ireal * restrict _vx, * restrict _vxend, * restrict _mvx,
                  * restrict _evx, * restrict _px;
    register ireal lax, dt2, etaxdt, vx;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];

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
    
    #pragma omp parallel private(tsz,tid,iy,iz,_vx,_vxend,_mvx,_evx,_px,vx,etaxdt,dims,i)
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
        _px  = s[D_P0 ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + 1;
        _evx = s[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
        
        tid -= gzs; /* adjust back after calculation of pointers */
        
        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                for ( _vxend = _vx + nx; _vx < _vxend; )
                {
                    etaxdt = (*_evx++) * dt2;
                
					vx = 0.0;
					for ( i = 0; i < k; ++i ) vx += (_px[i] - _px[-1-i]) * coeffs[i];
					++_px;
                    
                    (*_vx) = ((*_vx) * (1.0 - etaxdt) + vx * lax * (*_mvx++)) / (1.0 + etaxdt);
                    ++_vx;
                }
                
                _vx += vx_a; _mvx += mvx_a; _px += px_a; _evx -= nx;
            }
            
            _vx += vx_aa; _mvx += mvx_aa; _px += px_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts3d_2Kvj(RDOM *dom, void *pars, int ind)
{
	int i, k, ps;
	const ireal *coeffs;
    int nx, ny, nz, gxs, gys, gzs, v_a, mv_a, p_a, ev_a, v_aa, mv_aa, p_aa,
        ev_aa, iy, iz, tsz, tid, D_v, D_mv, D_ev, D_p/*, p_shift - replaced with ps */;
    register ireal * restrict _v, * restrict _vend, * restrict _mv,
                  * restrict _p;
    ireal * restrict _ev;
    register ireal la, etadt, v;
    ireal dt2;
    RARR *s;
    INFODIM *dims;
    
	k = ((SGN_TS_PARS*)pars)->k;
	if ( (k < 1) || (k > MAXK) ) return E_BADINPUT;
	coeffs = SCHEME_COEFFS[k - 1];

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
    ev_a = ((ind == 1) ? 1 : 0);
    
    la = ((SGN_TS_PARS*)pars)->lam[ind];
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

    ps = s[D_p]._dims[0].n0 * ((ind == 1) ?  1 : s[D_p]._dims[1].n0);
    
    #pragma omp parallel private(tsz,tid,iy,iz,_v,_vend,_mv,_ev,_p,v,etadt,dims,i)
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
            ev_aa = ((ind == 1) ? -ny : tsz);
        }
        #pragma omp barrier
        
        tid += gzs; /* adjust for shorter calculation of pointers */
        
        dims = s[D_v ]._dims;
        _v   = s[D_v ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;
        dims = s[D_mv]._dims;
        _mv  = s[D_mv]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0;        
        dims = s[D_p ]._dims;
        _p   = s[D_p ]._s + (gxs - dims[0].gs) + (gys - dims[1].gs + (tid - dims[2].gs) * dims[1].n0) * dims[0].n0 + ps;
        _ev  = s[D_ev]._s + ((ind == 1) ? gys : tid) - s[D_ev]._dims[0].gs ;

        tid -= gzs; /* adjust back after calculation of pointers */

        for ( iz = tid; iz < nz; iz += tsz )
        {
            for ( iy = 0; iy < ny; ++iy )
            {
                etadt = (*_ev) * dt2;

                for ( _vend = _v + nx; _v < _vend; )
                {
					v = 0.0;
					for ( i = 0; i < k; ++i ) v += (_p[i * ps] - _p[(-1-i) * ps]) * coeffs[i];
					++_p;
                    (*_v) = ((*_v) * (1.0 - etadt) + v * la * (*_mv++)) / (1.0 + etadt);
                    ++_v;
                }

                _v += v_a; _mv += mv_a; _p += p_a; _ev += ev_a;
            }

            _v += v_aa; _mv += mv_aa; _p += p_aa; _ev += ev_aa;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/
