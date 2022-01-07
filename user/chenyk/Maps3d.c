/* 3D acoustic wavefield modeling using the pseudo-spectral method  
DEMO:
https://github.com/chenyk1990/tutorials/blob/main/demo/aps3d/SConstruct
*/
/*
  written by Chen et al., 2020
  
  Copyright: The University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
  Reference:
  Chen, Y., O.M. Saad, M. Bai, X. Liu, and S. Fomel, 2021, A compact program for 3D passive seismic source-location imaging, Seismological Research Letters, doi: 10.1785/0220210050.
  
  The fully reproducible package (including scripts for all data examples in the paper) is at https://github.com/chenyk1990/reproducible_research/tree/master/passive. 
*/

#include <rsf.h>

/** Part I: Ricker wavelet ********/
float Ricker(float t, float f0, float t0, float A) 
/*< ricker wavelet:
 * f0: peak frequency
 * t0: time lag
 * A: amplitude
 * ************************>*/
{
        float x=pow(SF_PI*f0*(t-t0),2);
        return -A*exp(-x)*(1-2*x);
}

static kiss_fft_cpx *shape;

void ricker_init(int nfft   /* time samples */, 
		 float freq /* frequency */,
		 int order  /* derivative order */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w;
    kiss_fft_cpx cw;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./(nfft*freq);
 
    shape = (kiss_fft_cpx*) sf_complexalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w *= w;

	switch (order) {
	    case 2: /* half-order derivative */
		cw.r = 2*SF_PI/nfft;
		cw.i = iw*2*SF_PI/nfft;
		cw = sf_csqrtf(cw);
		shape[iw].r = cw.r*w*expf(1-w)/nfft;
		shape[iw].i = cw.i*w*expf(1-w)/nfft;
		break;
	    case 0:
	    default:
		shape[iw].r = w*expf(1-w)/nfft;
		shape[iw].i = 0.;
		break;
	}
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_cset(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/** Part II: Absorbing boundary condition ********/
/*Note: more powerful and efficient ABC can be incorporated*/
static int nx, ny, nz, nx2, ny2, nz2, nbt, nbb, nblx, nbrx, nbly, nbry;
static float ct, cb, clx, crx, cly, cry;
static float *wt, *wb, *wlx, *wrx, *wly, *wry;

void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    /*const float pi=SF_PI;*/
    if(!nb) return;
    switch(abc) {
    default:
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
        for(ib=0; ib<nb; ib++){
	    w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	}
    }
}

void abc_init(int n1,  int n2, int n3,    /*model size*/
	      int n12, int n22, int n32,   /*padded model size*/
	      int nb1, int nb2,    /*top, bottom*/
	      int nb3, int nb4,   /*left x, right x*/
	      int nb5, int nb6,   /*left y, right y*/
	      float c1, float c2, /*top, bottom*/
	      float c3, float c4, /*left x, right x*/
	      float c5, float c6 /*left y, right y*/)
/*< initialization >*/
{
    int c;
    nz = n1;
    nx = n2;
    ny = n3;
    nz2= n12;
    nx2= n22;
    ny2= n32;
    nbt = nb1;
    nbb = nb2;
    nblx = nb3;
    nbrx = nb4;
    nbly = nb5;
    nbry = nb6;
    ct = c1;
    cb = c2;
    clx = c3;
    crx = c4;
    cly = c5;
    cry = c6;
    if(nbt) wt =  sf_floatalloc(nbt);
    if(nbb) wb =  sf_floatalloc(nbb);
    if(nblx) wlx =  sf_floatalloc(nblx);
    if(nbrx) wrx =  sf_floatalloc(nbrx);
    if(nbly) wly =  sf_floatalloc(nbly);
    if(nbry) wry =  sf_floatalloc(nbry);
    c=0;
    abc_cal(c,nbt,ct,wt);
    abc_cal(c,nbb,cb,wb);
    abc_cal(c,nblx,clx,wlx);
    abc_cal(c,nbrx,crx,wrx);
    abc_cal(c,nbly,cly,wly);
    abc_cal(c,nbry,cry,wry);      
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(nbt) free(wt);
    if(nbb) free(wb);
    if(nblx) free(wlx);
    if(nbrx) free(wrx);
    if(nbly) free(wly);
    if(nbry) free(wry);
}

void abc_apply(float *a /*2-D matrix*/) 
/*< boundary decay>*/
{
    int i;
    int iz, iy, ix;
	
    /* top */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix,iy,i)
{
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nx2; ix++) {
        for (iy=0; iy < ny2; iy++) {
	  i = nz2*nx2*iy + nz2*ix + iz;
	  a[i] *= wt[iz];
        }
        }
    }
    
    /* bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nx2; ix++) {
        for (iy=0; iy < ny2; iy++) {
	  i = nz2*nx2*iy + nz2*ix + nz2-1-iz;
	  a[i] *= wb[iz];
        }
    }
    }
      
    /* left x*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nz2; iz++) {  
        for (ix=0; ix < nblx; ix++) {
        for (iy=0; iy < ny2; iy++) { 
	  i = nz2*nx2*iy+nz2*ix + iz;
	  a[i] *= wlx[ix];
        }
        }
    }
    
    /* right x*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nz2; iz++) {  
        for (ix=0; ix < nbrx; ix++) {
        for (iy=0; iy < ny2; iy++) {     
	  i = nz2*nx2*iy + nz2*(nx2-1-ix) + iz;
          a[i] *= wrx[ix];
        }
        }
    }
        
    /* left y*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nz2; iz++) {  
       for (ix=0; ix < nx2; ix++) {
        for (iy=0; iy < nbly; iy++) { 
	  i = nz2*nx2*iy+nz2*ix + iz;
	  a[i] *= wly[iy];
        }
        }
    }
        
    /* right y*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nz2; iz++) {  
       for (ix=0; ix < nx2; ix++) {
        for (iy=0; iy < nbry; iy++) {    
	  i = nz2*nx2*(ny2-1-iy) + nz2*ix + iz;
          a[i] *= wry[iy];
        }
        }
    }
#ifdef _OPENMP
}
#endif
}


/** Part III: Fourier transform ********/
static bool cmplx;
static int n1, n2, n3, nk;
static float wwt;

static float ***ff=NULL;
static sf_complex ***cc=NULL;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fftr_cfg cfg, icfg;
static kiss_fft_cfg cfg1, icfg1, cfg2, icfg2, cfg3, icfg3;
static kiss_fft_cpx ***tmp, *ctrace2, *ctrace3;
static sf_complex *trace2, *trace3;
#endif

int fft3_init(bool cmplx1        /* if complex transform */,
	      int pad1           /* padding on the first axis */,
	      int nx,   int ny,   int nz   /* axis 1,2,3; input data size */, 
	      int *nx2, int *ny2, int *nz2 /* axis 1,2,3; padded data size */)
/*< initialize >*/
{
#ifndef SF_HAS_FFTW
    int i2, i3;
#endif

    cmplx = cmplx1;

    /* axis 1 */

    if (cmplx) {
	nk = n1 = kiss_fft_next_fast_size(nx*pad1);

#ifndef SF_HAS_FFTW
	cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
	icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
#endif
    } else {
	nk = kiss_fft_next_fast_size(pad1*(nx+1)/2)+1;
	n1 = 2*(nk-1);

#ifndef SF_HAS_FFTW
	cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
	icfg = kiss_fftr_alloc(n1,1,NULL,NULL);
#endif
    }

    /* axis 2 */

    n2 = kiss_fft_next_fast_size(ny);

#ifndef SF_HAS_FFTW
    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);

    trace2 = sf_complexalloc(n2);
    ctrace2 = (kiss_fft_cpx *) trace2;
#endif

    /* axis 3 */

    n3 = kiss_fft_next_fast_size(nz);

#ifndef SF_HAS_FFTW
    cfg3  = kiss_fft_alloc(n3,0,NULL,NULL);
    icfg3 = kiss_fft_alloc(n3,1,NULL,NULL);

    trace3 = sf_complexalloc(n3);
    ctrace3 = (kiss_fft_cpx *) trace3;

    /* --- */

    tmp = (kiss_fft_cpx***) sf_alloc (n3,sizeof(kiss_fft_cpx**));
    tmp[0] = (kiss_fft_cpx**) sf_alloc (n2*n3,sizeof(kiss_fft_cpx*));
    tmp[0][0] = (kiss_fft_cpx*) sf_alloc (nk*n2*n3,sizeof(kiss_fft_cpx));

    for (i2=1; i2 < n2*n3; i2++) {
	tmp[0][i2] = tmp[0][0]+i2*nk;
    }

    for (i3=1; i3 < n3; i3++) {
	tmp[i3] = tmp[0]+i3*n2;
    }
#endif

    if (cmplx) {
	cc = sf_complexalloc3(n1,n2,n3);
    } else {
	ff = sf_floatalloc3(n1,n2,n3);
    }

    *nx2 = n1;
    *ny2 = n2;
    *nz2 = n3;

    wwt =  1.0/(n3*n2*n1);

    return (nk*n2*n3);
}

void fft3(float *inp      /* [n1*n2*n3] */, 
	  sf_complex *out /* [nk*n2*n3] */)
/*< 3-D FFT >*/
{
    int i1, i2, i3;
    float f;

  #ifdef SF_HAS_FFTW
    if (NULL==cfg) {
	cfg = cmplx? 
	    fftwf_plan_dft_3d(n3,n2,n1,
			      (fftwf_complex *) cc[0][0], 
			      (fftwf_complex *) out,
			      FFTW_FORWARD, FFTW_MEASURE):
	    fftwf_plan_dft_r2c_3d(n3,n2,n1,
				  ff[0][0], (fftwf_complex *) out,
				  FFTW_MEASURE);
	if (NULL == cfg) sf_error("FFTW failure.");
    }
#endif  
    
    /* FFT centering */    
    for (i3=0; i3<n3; i3++) {
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1<n1; i1++) {
		f = inp[(i3*n2+i2)*n1+i1];
		if (cmplx) {
		    cc[i3][i2][i1] = sf_cmplx((((i3%2==0)==(i2%2==0))==(i1%2==0))? f:-f,0.);
		} else {
		    ff[i3][i2][i1] = ((i3%2==0)==(i2%2==0))? f:-f;
		}
	    }
	}
    }

#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);
#else

    /* FFT over first axis */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if (cmplx) {
		kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i3][i2],tmp[i3][i2],1);
	    } else {
		kiss_fftr (cfg,ff[i3][i2],tmp[i3][i2]);
	    }
	}
    }

    /* FFT over second axis */
    for (i3=0; i3 < n3; i3++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(cfg2,tmp[i3][0]+i1,ctrace2,nk);
	    for (i2=0; i2 < n2; i2++) {
		tmp[i3][i2][i1]=ctrace2[i2];
	    }
	}
    }

    /* FFT over third axis */
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(cfg3,tmp[0][0]+i2*nk+i1,ctrace3,nk*n2);
	    for (i3=0; i3<n3; i3++) {
		out[(i3*n2+i2)*nk+i1] = trace3[i3];
	    }
	}
    } 
   
#endif

}

void ifft3_allocate(sf_complex *inp /* [nk*n2*n3] */)
/*< allocate inverse transform >*/
{
#ifdef SF_HAS_FFTW
    icfg = cmplx? 
	fftwf_plan_dft_3d(n3,n2,n1,
			  (fftwf_complex *) inp, 
			  (fftwf_complex *) cc[0][0],
			  FFTW_BACKWARD, FFTW_MEASURE):
	fftwf_plan_dft_c2r_3d(n3,n2,n1,
			      (fftwf_complex *) inp, ff[0][0],
			      FFTW_MEASURE);
    if (NULL == icfg) sf_error("FFTW failure.");
 #endif
}

void ifft3(float *out      /* [n1*n2*n3] */, 
	   sf_complex *inp /* [nk*n2*n3] */)
/*< 3-D inverse FFT >*/
{
    int i1, i2, i3;

#ifdef SF_HAS_FFTW
    fftwf_execute(icfg);
#else

    /* IFFT over third axis */
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(icfg3,(kiss_fft_cpx *) (inp+i2*nk+i1),ctrace3,nk*n2);
	    for (i3=0; i3<n3; i3++) {
		tmp[i3][i2][i1] = ctrace3[i3];
	    }
	}
    }
    
    /* IFFT over second axis */
    for (i3=0; i3 < n3; i3++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(icfg2,tmp[i3][0]+i1,ctrace2,nk);		
	    for (i2=0; i2<n2; i2++) {
		tmp[i3][i2][i1] = ctrace2[i2];
	    }
	}
    }

    /* IFFT over first axis */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if (cmplx) {
		kiss_fft_stride(icfg1,tmp[i3][i2],(kiss_fft_cpx *) cc[i3][i2],1);		
	    } else {
		kiss_fftri(icfg,tmp[i3][i2],ff[i3][i2]);
	    }
	}
    }

#endif

    /* FFT centering and normalization */
    for (i3=0; i3<n3; i3++) {
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1<n1; i1++) {
		if (cmplx) {
		    out[(i3*n2+i2)*n1+i1] = ((((i3%2==0)==(i2%2==0))==(i1%2==0))? wwt:-wwt)*crealf(cc[i3][i2][i1]);
		} else {
		    out[(i3*n2+i2)*n1+i1] = (((i3%2==0)==(i2%2==0))? wwt: - wwt)*ff[i3][i2][i1];
		}
	    }
	}
    }
}


/** Part IV: pseudo-spectral wave extrapolation ********/
typedef struct Psmpar {
  /*survey parameters*/
  int   nx, ny, nz;
  float dx, dy, dz;
  int   ns;
  int   *spx, *spy, *spz;
  int   gpz, gpx, gpy, gplx, gply;
  int   gpz_v, gpx_v, gpl_v;
  int   jsnap;
  /*fft related*/
  bool  cmplx;
  int   pad1;
  /*absorbing boundary*/
  bool abc;
  int nbt, nbb, nblx, nbrx, nbly, nbry;
  float ct,cb,clx,crx,cly,cry;
  /*source parameters*/
  int src; /*source type*/
  int nt;
  float dt,*f0,*t0,*A;
  /*misc*/
  bool verb, ps;
  float vref;
} * psmpar; /*psm parameters*/
/*^*/


int psm(float **wvfld, float ***dat, float **dat_v, float *img, float *vel, psmpar par, bool tri)
/*< pseudo-spectral method >*/
{
    /*survey parameters*/
    int   nx, ny, nz;
    float dx, dy, dz;
    int   ns;
    int   *spx, *spy, *spz;
    int   gpz, gpy, gpx, gplx, gply;
    int   gpz_v, gpx_v, gpl_v;
    int   jsnap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nblx, nbrx, nbly, nbry;
    float ct,cb,clx,crx,cly,cry;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, ny1, nz1; /*domain of interest*/
    int it,iz,ik,ix,iy,i,j;     /* index variables */
    int nk,nzxy,nz2,nx2,ny2,nzxy2,nkz,nkx,nth;
    int it1, it2, its;
    float dkx,dky,dkz,kx0,ky0,kz0,vref2,kx,ky,kz,k,t;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*prev,*lapl;

    /*source*/
    float **rick;
    float freq;
    int fft_size;

    /*passing the parameters*/
    nx    = par->nx;
    ny    = par->ny;
    nz    = par->nz;
    dx    = par->dx;
    dy    = par->dy;
    dz    = par->dz;
    ns= par->ns;
    spx   = par->spx;
    spy   = par->spy;
    spz   = par->spz;
    gpz   = par->gpz;
    gpy   = par->gpy;
    gpx   = par->gpx;
    gplx   = par->gplx;
    gply   = par->gply;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    jsnap  = par->jsnap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nblx   = par->nblx;
    nbrx   = par->nbrx;
    nbly   = par->nbly;
    nbry   = par->nbry;
    ct    = par->ct;
    cb    = par->cb;
    clx    = par->clx;
    crx    = par->crx;
    cly    = par->cly;
    cry    = par->cry;
    src   = par->src;
    nt    = par->nt;
    dt    = par->dt;
    f0    = par->f0;
    t0    = par->t0;
    A     = par->A;
    verb  = par->verb;
    ps    = par->ps;
    vref  = par->vref;
    

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
#else
    nth = 1;
#endif
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nblx-nbrx;
    ny1 = ny-nbly-nbry;
    
    nk = fft3_init(cmplx,pad1,nz,nx,ny,&nz2,&nx2,&ny2);
    nzxy = nz*nx*ny;
    nzxy2 = nz2*nx2*ny2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    dky = 1./(ny2*dy); ky0 = -0.5/dy;
    nkz = (cmplx)? nz2:(nz2/2+1);
    nkx = (cmplx)? nx2:(nx2/2+1);
    
    if(nk!=ny2*nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,dky=%f,kz0=%f,kx0=%f,ky0=%f",dkz,dkx,dky,kz0,kx0,ky0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d,ny2=%d",nk,nkz,nz2,nx2,ny2);

    if(abc)
      abc_init(nz,nx,ny,nz2,nx2,ny2,nbt,nbb,nblx,nbrx,nbly,nbry,ct,cb,clx,crx,cly,cry);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzxy); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzxy2);
    curr   = sf_floatalloc(nzxy2);
    prev   = sf_floatalloc(nzxy2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    if (!tri && src==0) {

      rick = sf_floatalloc2(nt,ns);
      for (i=0; i<ns; i++) {
	for (it=0; it<nt; it++) {
	  rick[i][it] = 0.f;
	}
	rick[i][(int)(t0[i]/dt)] = A[i]; /*time delay*/
	freq = f0[i]*dt;           /*peak frequency*/
	fft_size = 2*kiss_fft_next_fast_size((nt+1)/2);
	ricker_init(fft_size, freq, 0);
	sf_freqfilt(nt,rick[i]);
	ricker_close();
      }
    } else{
    	 rick = NULL;}

    for (iz=0; iz < nzxy; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzxy2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (iy=0; iy < ny2; iy++) {
	ky = ky0+iy*dky;
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(ky,hypot(kx,kz));
	    if (ps) lapl[iz+ix*nkz+iy*nkz*nx2] = -k*k;
	    else lapl[iz+ix*nkz+iy*nkz*nx2] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }
    }

    if (tri) { /* time-reversal imaging */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* modeling */
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);

	/* matrix multiplication */
	fft3(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft3(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iy,ix,iz,i,j,old,c)
#endif
	for (iy = 0; iy < ny; iy++) {
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz+iy*nz*nx;  /* original grid */
		j = iz+ix*nz2+iy*nz2*nx2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}
	}

	if (tri) {
	  /* inject data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iy,ix)
#endif
	    for (iy = 0; iy < gply; iy++) {
	    for (ix = 0; ix < gplx; ix++) {
	      curr[gpz+(ix+gpx)*nz2+(iy+gpy)*nz2*nx2] += vv[gpz+(ix+gpx)*nz+(iy+gpy)*nz*nx]*dat[iy][ix][it];
	    }
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      curr[gpz_v+iz+(gpx_v)*nz2] += vv[gpz_v+iz+(gpx_v)*nz]*dat_v[iz][it];
	    }
	  }
	} else {
	  t = it*dt;
	  for (i=0; i<ns; i++) {
	    for(iy=-1;iy<=1;iy++) {
	    for(ix=-1;ix<=1;ix++) {
	      for(iz=-1;iz<=1;iz++) {
		ik = spz[i]+iz+nz*(spx[i]+ix)+nz*nx*(spy[i]+iy);
		j = spz[i]+iz+nz2*(spx[i]+ix)+nz2*nx2*(spy[i]+iy);
		if (src==0) {
		  curr[j] += vv[ik]*rick[i][it]/(abs(ix)+abs(iy)+abs(iz)+1);
		} else {
		  curr[j] += vv[ik]*Ricker(t, f0[i], t0[i], A[i])/(abs(ix)+abs(iy)+abs(iz)+1);
		}
	      }
	    }
	    }
	  }
	}
	
	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}

	if (!tri) {
	  /* record data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (iy = 0; iy < gply; iy++) {
	    for (ix = 0; ix < gplx; ix++) {
	      dat[iy][ix][it] = curr[gpz+(ix+gpx)*nz2+(iy+gpy)*nz2*nx2];
	    }
	    }
	  }
	  if (NULL!=dat_v) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iz)
#endif
	    for (iz = 0; iz < gpl_v; iz++) {
	      dat_v[iz][it] = curr[gpz_v+iz+(gpx_v)*nz2];
	    }
	  }
	}
		
	/* save wavefield */
	if (jsnap > 0 && it%jsnap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iy,ix,iz,i,j)
#endif
	  for (iy=0; iy<ny1; iy++) {
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix + nz1*nx1*iy;
	      j = iz+nbt + (ix+nblx)*nz2 + (iy+nbly)*nz2*nx2; /* padded grid */
	      wvfld[it/jsnap][i] = curr[j];
	    }
	  }
	}
	}
    }
    if(verb) sf_warning(".");
    if (tri) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(iy,ix,iz)
#endif
    for (iy = 0; iy < ny1; iy++) {
    for (ix = 0; ix < nx1; ix++) {
	for (iz = 0; iz < nz1; iz++) {
	  img[iz + nz1*ix + nz1*nx1*iy] = curr[iz+nbt + (ix+nblx)*nz2 + (iy+nbly)*nz2*nx2];
	}
    }
    }
    }

    /*free up memory*/
//     fft3_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}

/** Part V: main program ********/
int main(int argc, char* argv[])
{

    /*survey parameters*/
    int   nx, ny, nz;
    float dx, dy, dz;
    int   ns;
    int   *spx, *spy, *spz;
    int   gpz, gpx, gpy, gplx, gply; /*geophone positions (z,x,y) and geophone length (z,x,y)*/
    int   gpz_v, gpx_v, gpy_v, gpl_v;
    int   jsnap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nblx, nbrx, nbly, nbry; /*boundaries for top/bottom, left/right x, left/right y*/
    float ct,cb,clx,crx,cly,cry; 		  /*decaying parameter*/
    /*source parameters*/
    int src; /*source type*/
    int nt,ntsnap;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps, tri; /*tri: time-reversal imaging*/
    float vref;

    psmpar par;
    int nx1, ny1, nz1; /*domain of interest*/
    int it;
    float *vel,***dat,**dat_v,**wvfld,*img; /*velocity profile*/
    sf_file Fi,Fo,Fd,Fd_v,snaps; /* I/O files */
    sf_axis az,ax,ay; /* cube axes */

    sf_init(argc,argv);

    if (!sf_getint("jsnap",&jsnap)) jsnap=0; /* interval for snapshots */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex fft */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if(!sf_getbool("abc",&abc)) abc=false; /* absorbing flag */
    if (abc) {
      if(!sf_getint("nbt",&nbt)) sf_error("Need nbt!");
      if(!sf_getint("nbb",&nbb)) nbb = nbt;
      if(!sf_getint("nblx",&nblx)) nblx = nbt;
      if(!sf_getint("nbrx",&nbrx)) nbrx = nbt;
      if(!sf_getint("nbly",&nbly)) nbly = nbt;
      if(!sf_getint("nbry",&nbry)) nbry = nbt;
      if(!sf_getfloat("ct",&ct)) sf_error("Need ct!");
      if(!sf_getfloat("cb",&cb)) cb = ct;
      if(!sf_getfloat("clx",&clx)) clx = ct;
      if(!sf_getfloat("crx",&crx)) crx = ct;
      if(!sf_getfloat("cly",&cly)) cly = ct;
      if(!sf_getfloat("cry",&cry)) cry = ct;
    } else {
      nbt = 0; nbb = 0; nblx = 0; nbrx = 0; nbly = 0; nbry = 0;
      ct = 0; cb = 0; clx = 0; crx = 0; cly = 0; cry = 0;
    }
    if (!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getbool("ps",&ps)) ps=false; /* use pseudo-spectral */
    if (ps) sf_warning("Using pseudo-spectral...");
    else sf_warning("Using pseudo-analytical...");
    if (!sf_getbool("tri",&tri)) tri=false; /* if choose time reversal imaging */
    if (tri) sf_warning("Time-reversal imaging");
    else sf_warning("Forward modeling");
    if (!sf_getfloat("vref",&vref)) vref=1500; /* reference velocity (default using water) */

    /* setup I/O files */
    Fi = sf_input ("in");
    Fo = sf_output("out");
    if (tri) {
      gplx = -1;
      gply = -1;
      gpl_v = -1;
      if (NULL==sf_getstring("dat") && NULL==sf_getstring("dat_v"))
	sf_error("Need Data!");
      if (NULL!=sf_getstring("dat")) {
	Fd = sf_input("dat");
	sf_histint(Fd,"n1",&nt);
	sf_histfloat(Fd,"d1",&dt);
	sf_histint(Fd,"n2",&gplx);
	sf_histint(Fd,"n3",&gply);
      } else Fd = NULL;
      if (NULL!=sf_getstring("dat_v")) {
	Fd_v = sf_input("dat_v");
	sf_histint(Fd_v,"n1",&nt);
	sf_histfloat(Fd_v,"d1",&dt);
	sf_histint(Fd_v,"n2",&gpl_v);
      } else Fd_v = NULL;
      src = -1; ns = -1;
      spx = NULL; spy = NULL; spz = NULL;
      f0 = NULL; t0 = NULL; A = NULL;
    } else {
      Fd = NULL;
      if (!sf_getint("nt",&nt)) sf_error("Need nt!");
      if (!sf_getfloat("dt",&dt)) sf_error("Need dt!");
      if (!sf_getint("gplx",&gplx)) gplx = -1; /* geophone length X*/
      if (!sf_getint("gply",&gply)) gply = -1; /* geophone length Y*/
      if (!sf_getint("gpl_v",&gpl_v)) gpl_v = -1; /* geophone height */
      if (!sf_getint("src",&src)) src=0; /* source type */
      if (!sf_getint("ns",&ns)) ns=1; /* source type */
      spx = sf_intalloc(ns);
      spy = sf_intalloc(ns);
      spz = sf_intalloc(ns);
      f0  = sf_floatalloc(ns);
      t0  = sf_floatalloc(ns);
      A   = sf_floatalloc(ns);
      if (!sf_getints("spx",spx,ns)) sf_error("Need spx!"); /* shot position x */
      if (!sf_getints("spy",spy,ns)) sf_error("Need spy!"); /* shot position y */
      if (!sf_getints("spz",spz,ns)) sf_error("Need spz!"); /* shot position z */
      if (!sf_getfloats("f0",f0,ns)) sf_error("Need f0! (e.g. 30Hz)");   /*  wavelet peak freq */
      if (!sf_getfloats("t0",t0,ns)) sf_error("Need t0! (e.g. 0.04s)");  /*  wavelet time lag */
      if (!sf_getfloats("A",A,ns)) sf_error("Need A! (e.g. 1)");     /*  wavelet amplitude */
    }
    if (!sf_getint("gpx",&gpx)) gpx = -1; /* geophone position x */
    if (!sf_getint("gpy",&gpy)) gpy = -1; /* geophone position y */
    if (!sf_getint("gpz",&gpz)) gpz = -1; /* geophone position z */
    if (!sf_getint("gpx_v",&gpx_v)) gpx_v = -1; /* geophone position x */
    if (!sf_getint("gpy_v",&gpy_v)) gpy_v = -1; /* geophone position y */
    if (!sf_getint("gpz_v",&gpz_v)) gpz_v = -1; /* geophone position z */

    if (SF_FLOAT != sf_gettype(Fi)) sf_error("Need float input");

    /* Read/Write axes */
    az = sf_iaxa(Fi,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fi,2); nx = sf_n(ax); dx = sf_d(ax);
    ay = sf_iaxa(Fi,3); ny = sf_n(ay); dy = sf_d(ay);

    nz1 = nz-nbt-nbb;
    nx1 = nx-nblx-nbrx;
    ny1 = ny-nbly-nbry;
    if(verb)sf_warning("ny=%d,nbly=%d,nbry=%d",ny,nbly,nbry);
    if(verb)sf_warning("nz1=%d,nx1=%d,ny1=%d",nz1,nx1,ny1);
    if (gpx==-1) gpx = nblx;
    if (gpy==-1) gpy = nbly;
    if (gpz==-1) gpz = nbt;
    if (gplx==-1) gplx = nx1;
    if (gply==-1) gply = ny1;
    if (gpx_v==-1) gpx_v = nblx;
    if (gpy_v==-1) gpy_v = nbly;
    if (gpz_v==-1) gpz_v = nbt;
    if (gpl_v==-1) gpl_v = nz1;
    ntsnap=0;
    if (jsnap)
        for (it=0;it<nt;it++)
            if (it%jsnap==0) ntsnap++;
    if (tri) { /*output final wavefield*/
      sf_setn(az,nz1);
      sf_setn(ax,nx1);
      sf_setn(ay,ny1);
      sf_oaxa(Fo,az,1);
      sf_oaxa(Fo,ax,2);
      sf_oaxa(Fo,ay,3);   
      sf_settype(Fo,SF_FLOAT);
    } else { /*output data*/
      sf_setn(ax,gplx);
      sf_setn(ay,gply);
      sf_putint(Fo,"n3",gply);
      sf_warning("ny=%d,nbly=%d,nbry=%d",ny,nblx,nbly);
      sf_warning("gplx=%d,gply=%d",gplx,gply);
      /*output horizontal data is mandatory*/
      sf_putint(Fo,"n1",nt);
      sf_putfloat(Fo,"d1",dt);
      sf_putfloat(Fo,"o1",0.);
      sf_putstring(Fo,"label1","Time");
      sf_putstring(Fo,"unit1","s");
      sf_oaxa(Fo,ax,2);
      sf_settype(Fo,SF_FLOAT);
      /*output vertical data is optional*/
      if (NULL!=sf_getstring("dat_v")) {
	Fd_v = sf_output("dat_v");
	sf_setn(az,gpl_v);
	sf_putint(Fd_v,"n1",nt);
	sf_putfloat(Fd_v,"d1",dt);
	sf_putfloat(Fd_v,"o1",0.);
	sf_putstring(Fd_v,"label1","Time");
	sf_putstring(Fd_v,"unit1","s");
	sf_oaxa(Fd_v,az,2);
	sf_putint(Fd_v,"n3",1);
	sf_settype(Fd_v,SF_FLOAT);	
      } else Fd_v = NULL;
    }

    if (jsnap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	sf_setn(az,nz1);
	sf_setn(ax,nx1);
	sf_setn(ay,ny1);
	sf_oaxa(snaps,az,1);
	sf_oaxa(snaps,ax,2);
	sf_oaxa(snaps,ay,3);
	sf_putint(snaps,"n4",ntsnap);
	sf_putfloat(snaps,"d4",dt*jsnap);
	sf_putfloat(snaps,"o4",0.);
	sf_putstring(snaps,"label4","Time");
	sf_putstring(snaps,"unit4","s");
    } else snaps = NULL;

    par = (psmpar) sf_alloc(1,sizeof(*par));
    vel = sf_floatalloc(nz*ny*nx);
    
    if (tri && NULL==Fd) {dat = NULL;  }
    else { dat = sf_floatalloc3(nt,gplx,gply);}

    if (NULL!=Fd_v) dat_v = sf_floatalloc2(nt,gpl_v);
    else dat_v = NULL;

    if (tri) img = sf_floatalloc(nz1*ny1*nx1);
    else img = NULL;

    if (jsnap>0) wvfld = sf_floatalloc2(nx1*ny1*nz1,ntsnap);
    else wvfld = NULL;
    

    sf_floatread(vel,nz*ny*nx,Fi);

    if (tri) {
      if (NULL!=Fd)   sf_floatread(dat[0][0],gplx*gply*nt,Fd);
      if (NULL!=Fd_v) sf_floatread(dat_v[0],gpl_v*nt,Fd_v);
    }
    /*passing the parameters*/
    par->nx    = nx;  
    par->ny    = ny;
    par->nz    = nz;
    par->dx    = dx;
    par->dy    = dy;
    par->dz    = dz;
    par->ns	   = ns;
    par->spx   = spx;
    par->spy   = spy;
    par->spz   = spz;
    par->gpx   = gpx;
    par->gpy   = gpy;
    par->gpz   = gpz;
    par->gplx   = gplx;
    par->gply   = gply;
    par->gpz_v = gpz_v;
    par->gpx_v = gpx_v;
    par->gpl_v = gpl_v;
    par->jsnap  = jsnap;
    par->cmplx = cmplx;
    par->pad1  = pad1;
    par->abc   = abc;
    par->nbt   = nbt;
    par->nbb   = nbb;
    par->nblx   = nblx;
    par->nbrx   = nbrx;
    par->nbly   = nbly;
    par->nbry   = nbry;
    par->ct    = ct;
    par->cb    = cb;
    par->clx    = clx;
    par->crx    = crx;
    par->cly    = cly;
    par->cry    = cry;
    par->src   = src;
    par->nt    = nt;
    par->dt    = dt;
    par->f0    = f0;
    par->t0    = t0;
    par->A     = A;
    par->verb  = verb;
    par->ps    = ps;
    par->vref  = vref;

    /*do the work*/
    psm(wvfld, dat, dat_v, img, vel, par, tri);

    if (tri) {
      sf_floatwrite(img,nz1*ny1*nx1,Fo);
    } else {
      sf_floatwrite(dat[0][0],gplx*gply*nt,Fo);
      if (NULL!=Fd_v)
	sf_floatwrite(dat_v[0],gpl_v*nt,Fd_v);
    }

    if (jsnap>0)
      sf_floatwrite(wvfld[0],nz1*nx1*ny1*ntsnap,snaps);
    
    exit (0);
}
