/* 
sgn210_1d.c
Igor Terentyev.
********************************************************************************
Implementation of 2-10 scheme in 1D.
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

int sgn_ts1d_210(RDOM *dom, int iarr, void *pars)
{
    /* rx - recomputed array, qx - participating array */
    int nx, gxs, D_rx, D_mrx, D_qx, D_erx, qx_shift;
    long ltsz, ltid, lthr_shift;     /* thread related */
    register ireal * restrict _rx, * restrict _rxend, * restrict _mrx,
                  * restrict _erx, * restrict _qx9;
    register ireal lax, dt2, qx9, qx8, qx7, qx6, qx5, qx4, qx3, qx2, qx1, qx0, etaxdt;
    RARR *s;
    
    /*
    if ( ((SGN_TS_PARS*)pars)->eflags[0] == 0 )  
        return sg_ts1d_210(dom, iarr,&(((SGN_TS_PARS*)pars)->sgpars));
    */
    if ( iarr == D_P0 )
    {
        D_rx  = D_P0;
        D_mrx = D_MP0;
        D_erx = D_EP0;
        D_qx  = D_V0;
        qx_shift = 0;
    }
    else if ( iarr == D_V0 )
    {
        D_rx  = D_V0;
        D_mrx = D_MV0;
        D_erx = D_EV0;
        D_qx  = D_P0;
        qx_shift = 1;
    }
    else return 0;

    s = dom->_s;
    
    nx = s[D_rx]._dims[0].n;
    if ( nx == 0 ) return 0;
            
    gxs = s[D_rx]._dims[0].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0];
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

    #pragma omp parallel private(ltsz,ltid,lthr_shift,_rx,_rxend,_mrx,_qx9,_erx,qx9,qx8,qx7,qx6,qx5,qx4,qx3,qx2,qx1,qx0,etaxdt)
    {
        #ifdef _OPENMP
        ltsz = (long)omp_get_num_threads();
        ltid = (long)omp_get_thread_num();
        #else
        ltsz = 1L;
        ltid = 0L;
        #endif
        
        lthr_shift = (ltid * (long)nx) / ltsz;
        
        _rx  = s[D_rx ]._s + (gxs - s[D_rx ]._dims[0].gs) + lthr_shift;
        _mrx = s[D_mrx]._s + (gxs - s[D_mrx]._dims[0].gs) + lthr_shift;
        _erx = s[D_erx]._s + (gxs - s[D_erx]._dims[0].gs) + lthr_shift;
        _qx9 = s[D_qx ]._s + (gxs - s[D_qx ]._dims[0].gs) + lthr_shift + qx_shift + 4;
        
        lthr_shift =  ((ltid + 1L) * (long)nx) / ltsz - lthr_shift;
        _rxend = _rx + lthr_shift;
        
        qx8 = _qx9[-1]; qx7 = _qx9[-2]; qx6 = _qx9[-3];
        qx5 = _qx9[-4]; qx4 = _qx9[-5]; qx3 = _qx9[-6];
        qx2 = _qx9[-7]; qx1 = _qx9[-8]; qx0 = _qx9[-9];
        
        while ( _rx < _rxend )
        {
            qx9 = *_qx9++;
            etaxdt = (*_erx++) * dt2;
            
            (*_rx) = ((*_rx) * (1.0 - etaxdt) + 
                      ( (qx9 - qx0) * C5 + (qx8 - qx1) * C4 + (qx7 - qx2) * C3 + (qx6 - qx3) * C2 + (qx5 - qx4) * C1 ) *
			          lax * (*_mrx++)) / (1.0 + etaxdt);
            _rx++;
            
            qx0 = qx1; qx1 = qx2; qx2 = qx3; qx3 = qx4; qx4 = qx5; qx5 = qx6; qx6 = qx7; qx7 = qx8; qx8 = qx9;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/
