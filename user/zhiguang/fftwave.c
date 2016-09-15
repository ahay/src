/* Simple 2-D wave propagation with multi-threaded fftw3 */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft2w.h"
#include "absorb.h"

#include "fftwave.h"

static int nx,nz,nk,nz2,nx2,nzx,nzx2,m2,pad1;
static bool abc,cmplx;
static float **wave, *curr, *prev;
static sf_complex *cwave, *cwavem;
static float **lt, **rt;
static bool taper;
static float *ktp;

int lrinit(int _nx, int _nz, int _m2, bool _cmplx, int _pad1, int _nb, float _c, bool _abc, float **_lt, float **_rt, bool _taper, float thresh, float dz, float dx)
/*<initialization>*/
{
  float kx,kz,ktmp;
  float dkz,dkx,kz0,kx0,kx_trs,kz_trs;
  int nkz;
  int iz,ix;

  nx = _nx; nz = _nz;
  m2 = _m2;
  cmplx = _cmplx;
  pad1 = _pad1;
  abc = _abc;
  lt = _lt; rt = _rt;
  taper = _taper;

  nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

  nzx = nz*nx;
  nzx2 = nz2*nx2;

  curr = sf_floatalloc(nzx2);
  prev = sf_floatalloc(nzx2);

  cwave  = sf_complexalloc(nk);
  cwavem = sf_complexalloc(nk);
  wave = sf_floatalloc2(nzx2,m2);

  if (taper) {
    ktp = sf_floatalloc(nk);
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    kx_trs = thresh*fabs(0.5/dx);
    kz_trs = thresh*fabs(0.5/dz);
    /* constructing the tapering op */
    for (ix=0; ix < nx2; ix++) {
      kx = kx0+ix*dkx;
      for (iz=0; iz < nkz; iz++) {
	kz = kz0+iz*dkz;
	ktmp = 1.;
	if (fabs(kx) > kx_trs)
	  ktmp *= powf((2*kx_trs - fabs(kx))/(kx_trs),2);
	if (fabs(kz) > kz_trs)
	  ktmp *= powf((2*kz_trs - fabs(kz))/(kz_trs),2);
	ktp[iz+ix*nkz] = ktmp;
      }
    }
  } else ktp = NULL;

  /* initialize abc */
  if (abc)
    abc_init(nz,nx,nz2,nx2,_nb,_nb,_nb,_nb,_c,_c,_c,_c);

  return 0;
}

int lrupdate(float *curr, float *prev)
/*<lowrank two-step update>*/
{
    int iz,ix,i,j,im,ik;     /* index variables */
    float c, old;
    
    /* matrix multiplication */
    fft2(curr,cwave);

    for (im = 0; im < m2; im++) {
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwave[ik]*rt[ik][im];
#else
	cwavem[ik] = sf_crmul(cwave[ik],rt[ik][im]);
#endif
      }
      ifft2(wave[im],cwavem);
    }

    for (ix = 0; ix < nx; ix++) {
      for (iz=0; iz < nz; iz++) {
	i = iz+ix*nz;  /* original grid */
	j = iz+ix*nz2; /* padded grid */
		
	old = c = curr[j];
	c += c - prev[j];
	prev[j] = old;

	for (im = 0; im < m2; im++) {
	  c += lt[im][i]*wave[im][j];
	}

	curr[j] = c;
      }	    
    }
    if (abc) {
      abc_apply(prev);
      abc_apply(curr);
    }

    if (taper) {
      fft2(curr,cwave);
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwave[ik]*ktp[ik];
#else
	cwavem[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
      }
      ifft2(curr,cwavem);
      fft2(prev,cwave);
      for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	cwavem[ik] = cwave[ik]*ktp[ik];
#else
	cwavem[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
      }
      ifft2(prev,cwavem);
    }
  
    return 0;
}

int lrclose()
/*<finalize>*/
{
    if (abc) abc_close();
    fft2_finalize();
    return 0;
}
