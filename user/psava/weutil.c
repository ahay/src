#include <rsf.h>

#include "slice.h"
/*^*/

#include "weutil.h"
/*^*/

#ifndef _weutil_h


/*------------------------------------------------------------*/
typedef struct cub *cub3d;
/*^*/

struct cub{
    bool    verb;
    sf_axa  amx,amy,amz;
    sf_axa  alx,aly;
    sf_axa  aw;
    sf_axa  ae;
    float   eps;
    int     ompnth;
};
/*^*/
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
    sf_axa     bxx,byy;
    float         **kk;
    int            *lx;
    int            *ly;
    sf_complex   ***pk;
    sf_complex   ***wk;
    float        ***wt;
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
    int       nrmax;
    float     dsmax;
};
/*^*/


/*------------------------------------------------------------*/
typedef struct sroperator *sroperator3d;
/*^*/

struct sroperator{
    sf_complex ***ww_s;
    sf_complex ***ww_r;
};
/*^*/

#endif
