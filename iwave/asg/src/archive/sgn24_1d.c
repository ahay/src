/* 
sgn24_1d.c
Igor Terentyev.
********************************************************************************
Implementation of 2-4 scheme in 1D.
*/
/*============================================================================*/

#include "utils.h"
#include "sgn.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/
/*
1/24 constant.
*/
#define C24 4.166666666666666666666667e-2
/*----------------------------------------------------------------------------*/

int sgn_ts1d_24(RDOM *dom, int iarr, void *pars)
{
    /* rx - recomputed array, qx - participating array */
    int nx, gxs, D_rx, D_mrx, D_qx, D_erx, qx_shift;
    long ltsz, ltid, lthr_shift;     /* thread related */
    register ireal * restrict _rx, * restrict _rxend, * restrict _mrx,
                  * restrict _erx, * restrict _qx3;
    register ireal lax, dt2, qx3, qx2, qx1, qx0, etaxdt;
    RARR *s;
    
    /*
    if ( ((SGN_TS_PARS*)pars)->eflags[0] == 0 )
        return sg_ts1d_24(dom, iarr, it, &(((SGN_TS_PARS*)pars)->sgpars));
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
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

    #pragma omp parallel private(ltsz,ltid,lthr_shift,_rx,_rxend,_mrx,_qx3,_erx,qx3,qx2,qx1,qx0,etaxdt)
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
        _qx3 = s[D_qx ]._s + (gxs - s[D_qx ]._dims[0].gs) + lthr_shift + qx_shift + 1;
        
        lthr_shift =  ((ltid + 1L) * (long)nx) / ltsz - lthr_shift;
        _rxend = _rx + lthr_shift;
        
        qx2 = _qx3[-1]; qx1 = _qx3[-2]; qx0 = _qx3[-3];
        
        while ( _rx < _rxend )
        {
            qx3 = *_qx3++;
            etaxdt = (*_erx++) * dt2;
            
            (*_rx) = ((*_rx) * (1.0 - etaxdt) + (qx3 - qx0 + (qx1 - qx2) * 27.0) *
			                                    lax * (*_mrx++)) / (1.0 + etaxdt);
            _rx++;
            
            qx0 = qx1; qx1 = qx2; qx2 = qx3;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/
