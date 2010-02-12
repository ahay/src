#include <rsf.h>

#include "wexutl.h"
/*^*/

#ifndef _wexutl_h

/*------------------------------------------------------------*/
typedef struct wexcub *wexcub3d;
/*^*/

struct wexcub{
    bool    verb;
    sf_axa  amx,amy;
    sf_axa  alx,aly;
    sf_axa  az;
    sf_axa  aw;
    sf_axa  ae;
    float   eps;
    int  ompnth;
};
/*^*/
/*------------------------------------------------------------*/
typedef struct wexf2d *wexfft2d;
/*^*/

struct wexf2d{
    int            n1,n2;
    kiss_fft_cfg   forw1; /*   FFT on axis 1 */
    kiss_fft_cfg   invs1;
    kiss_fft_cfg   forw2; /*   FFT on axis 2 */
    kiss_fft_cfg   invs2;
    kiss_fft_cpx  *ctmp1; /* temp array */
    kiss_fft_cpx  *ctmp2; /* temp array */
    float       fftscale; /* FFT scale  */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct wexssr *wexssr3d;
/*^*/

struct wexssr{
    sf_axa     bxx,byy;
    float         **kk;
    int        *lx,*ly;
    sf_complex   ***pk;
    sf_complex   ***wk;
    float        ***wt;
    float       dsmax2;
    wexfft2d      *f2d;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct wextap *wextap3d;
/*^*/

struct wextap{
    int     nt1,   nt2;
    int      n1,    n2;
    bool     b1,    b2;
    float *tap1, *tap2;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct wexslo *wexslo3d;
/*^*/

struct wexslo{
    sf_fslice slice; /* slowness slice */
    int         *nr; /* number of references */
    float      **sm; /* ref slo squared */
    int       nrmax;
    float     dsmax;
    float      ***s;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct wexop *wexop3d;
/*^*/

struct wexop{
    sf_complex ***ww  ; /*          wavefield */   
    sf_complex ****w;   /*          wavefield */     
};
/*^*/

#endif
