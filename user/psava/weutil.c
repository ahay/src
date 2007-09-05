#include <rsf.h>

#include "weutil.h"
/*^*/

#ifndef _weutil_h

/*------------------------------------------------------------*/
typedef struct fft *fft2d;
/*^*/

struct fft{
    int            n1,n2;
    kiss_fft_cfg  *forw1; /* FFT on axis 1 */
    kiss_fft_cfg  *invs1;
    sf_complex     *shf1;
    kiss_fft_cfg  *forw2; /* FFT on axis 2 */
    kiss_fft_cfg  *invs2;
    sf_complex     *shf2;
    kiss_fft_cpx**ctrace;
    float       fftscale;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct ssr *ssr3d;
/*^*/

struct ssr{
    sf_axa  az,axx,ayy;
    sf_axa     bxx,byy;
    sf_axa     alx,aly;
    float         **kk;
    int            *lx;
    int            *ly;
    sf_complex   ***pk;
    sf_complex   ***wk;
    float        ***wt;
    int         ompnth;
    float       dsmax2;
    fft2d          fft;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct tap *tap3d;
/*^*/

struct tap{
    int nt1, nt2, nt3;
    int  n1,  n2,  n3;
    bool b1,  b2,  b3;
    float *tap1;
    float *tap2;
    float *tap3;
};
/*^*/


/*------------------------------------------------------------*/
typedef struct sroperator *sroperator3d;
/*^*/

struct sroperator{
    ssr3d ssr;
    tap3d tap;
};
/*^*/

#endif
