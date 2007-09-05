#include <rsf.h>

#include "slice.h"
/*^*/

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
typedef struct slo *slo3d;
/*^*/

struct slo{
    fslice    slice; /* slowness slice */
    int         *nr; /* number of references */
    float      **sm; /* ref slo squared */
    float     ***ss; /* slowness */
    float     ***so; /* slowness */
    sf_axa      alx;
    sf_axa      aly;
    sf_axa      amz;
    int       nrmax;
    float     dsmax;
    int      ompnth;
};
/*^*/


/*------------------------------------------------------------*/
typedef struct sroperator *sroperator3d;
/*^*/

struct sroperator{
    bool verb;
    float eps;
    int ompnth;
    ssr3d ssr; // SSR operator
    tap3d tap; // taper parameters
    slo3d s_s; //   source slowness
    slo3d s_r; // receiver slowness
    sf_axa aw;
    sf_axa ae;
    sf_axa amz;
    sf_axa amx;
    sf_axa amy;
    sf_axa alx;
    sf_axa aly;
    sf_complex ***ww_s;
    sf_complex ***ww_r;
};
/*^*/

#endif
