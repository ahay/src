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
    int  ompnth;
    int  ompchunk;
};
/*^*/
/*------------------------------------------------------------*/
typedef struct fft *fft2d;
/*^*/

struct fft{
    int            n1,n2;
    kiss_fft_cfg  *forw1; /*   FFT on axis 1 */
    kiss_fft_cfg  *invs1;
    kiss_fft_cfg  *forw2; /*   FFT on axis 2 */
    kiss_fft_cfg  *invs2;
    sf_complex     *shf1; /* shift on axis 1 */
    sf_complex     *shf2; /* shift on axis 2 */
    kiss_fft_cpx**ctrace; /* temp array */
    float       fftscale; /* FFT scale */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct ssr *ssr3d;
/*^*/

struct ssr{
    sf_axa     bxx,byy;
    float         **kk;
    int        *lx,*ly;
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
    float    twoway; /* two-way traveltime */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct img *img3d;
/*^*/

struct img{
    fslice imag;
    fslice cigs;
    sf_complex ****qs;      /* source wavefield */
    sf_complex ****qr;      /* receiver wavefield */
    float       ***qi;      /* image */
    float       ***qc;      /* cigs */
    sf_axa  acx,acy,acz;
    int     jcx,jcy,jcz;
    sf_axa  ahx,ahy,ahz,aht;
    sf_axa  ahh,aha,ahb;
    sf_complex   **tt;      /* time-lag I.C. weight */
    int LOx,HIx;
    int LOy,HIy;
    int LOz,HIz;            /* space-lags I.C. */
    float vpvs;             /* abs-lag I.C. */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weoperator *weoperator3d;
/*^*/

struct weoperator{
    sf_complex ***ww_s; /*   source wavefield */
    sf_complex ***ww_r; /* receiver wavefield */
    sf_complex ***ww  ; /*          wavefield */     
    fslice        wtmp; /* tmp wavefield (for SR modeling) */
    float      ***rr;   /*  reflectivity (for SR modeling) */
    float       **qq;   /* image (for ZO modeling/migration) */
};
/*^*/

#endif
