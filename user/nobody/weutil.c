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
    sf_axa  ahx,ahy;
    sf_axa  alx,aly;
    sf_axa  aw;
    sf_axa  ae;
    float   eps;
    int  ompnth;
    int  ompchunk;
};
/*^*/
/*------------------------------------------------------------*/
typedef struct f2d *fft2d;
/*^*/

struct f2d{
    int            n1,n2;
    kiss_fft_cfg  *forw1; /*   FFT on axis 1 */
    kiss_fft_cfg  *invs1;
    kiss_fft_cfg  *forw2; /*   FFT on axis 2 */
    kiss_fft_cfg  *invs2;
    sf_complex     *shf1; /* shift on axis 1 */
    sf_complex     *shf2; /* shift on axis 2 */
    kiss_fft_cpx**ctrace; /* temp array */
    float       fftscale; /* FFT scale  */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct f3d *fft3d;
/*^*/

struct f3d{
    int            n1,n2,n3;
    kiss_fft_cfg  *forw1; /*   FFT on axis 1 */
    kiss_fft_cfg  *invs1;
    kiss_fft_cfg  *forw2; /*   FFT on axis 2 */
    kiss_fft_cfg  *invs2;
    kiss_fft_cfg  *forw3; /*   FFT on axis 3 */
    kiss_fft_cfg  *invs3;
    sf_complex     *shf1; /* shift on axis 1 */
    sf_complex     *shf2; /* shift on axis 2 */
    sf_complex     *shf3; /* shift on axis 3 */
    kiss_fft_cpx**ctrace2; /* temp array */
    kiss_fft_cpx**ctrace3; /* temp array */
    float       fftscale; /* FFT scale  */
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
    fft2d          f2d;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct lsr *lsr3d;
/*^*/

struct lsr{
    sf_axa     bxx,byy;
    float         **kk;
    float         **kw;
    int        *lx,*ly;
    sf_complex   ***wk;
    sf_complex   ***wt;
    float       dsmax2;
    fft2d          f2d;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct cam *cam3d;
/*^*/

struct cam{
    sf_axa     bmx,bmy,bhx;
    float     **ksx,**krx;
    int        **is,**ir;
    int         *jx,*jy;
    sf_complex  ****pk;
    sf_complex  ****wk;
    float       ****wt;
    float       dsmax2;
    fft3d          f3d;
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
    float      ****qi;      /* image */
    float      ****qc;      /* cigs */
    float      ****qt;
    sf_axa  acx,acy,acz;
    int     jcx,jcy,jcz;
    sf_axa  ahx,ahy,ahz,aht;
    sf_axa  ahh,aha,ahb;
    sf_complex   **tt;      /*  time-lag I.C. weight */
    int LOx,HIx;
    int LOy,HIy;
    int LOz,HIz;            /* space-lags I.C. */
    float vpvs;             /*    abs-lag I.C. */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct ssroperator *ssroperator3d;
/*^*/

struct ssroperator{
    sf_complex ***ww_s; /*   source wavefield */
    sf_complex ***ww_r; /* receiver wavefield */
    sf_complex ***ww  ; /*          wavefield */   
    fslice        wtmp; /* tmp wavefield (for SR modeling) */
    float      ***rr;   /*  reflectivity (for SR modeling)  */
    sf_complex  **qq;   /* image (for ZO modeling/migration) */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct ssrmvaoperator *ssrmvaoperator3d;
/*^*/

struct ssrmvaoperator{
    sf_complex ***bw; /* wavefield */
    sf_complex ***dw; /* wavefield */
    sf_complex ***pw;
    sf_complex ***pwsum;
    sf_complex ***ds; /* slowness */
    sf_complex ***ps;
    sf_complex ***pssum;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct camoperator *camoperator3d;
/*^*/

struct camoperator{
    sf_complex ****ww; /* wavefield */   
    sf_complex  ***qq;  /* image */   
};
/*^*/

#endif
