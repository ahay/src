/* Pseudo-spectral wave extrapolation for second-order two-way wave equation */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft2w.h"
#include "absorb.h"
#include "ricker.h"
#include "ricker1.h"
#include "psp.h"

#ifndef _psp_h

typedef struct Pspar {
  /*survey parameters*/
  int   nx, nz;
  float dx, dz;
  int   n_srcs;
  int   *spx, *spz;
  int   gpz, gpx, gpl;
  int   gpz_v, gpx_v, gpl_v;
  int   snap;
  /*fft related*/
  bool  cmplx;
  int   pad1;
  /*absorbing boundary*/
  bool abc;
  int nbt, nbb, nbl, nbr;
  float ct,cb,cl,cr;
  /*source parameters*/
  int src; /*source type*/
  int nt;
  float dt,*f0,*t0,*A;
  /*misc*/
  bool verb, ps;
  float vref;
} * pspar; /*psp parameters*/
/*^*/

#endif

int psp(float **wvfld, float **dat, float **dat_v, float *img, float *vel, pspar par, bool mig)
/*< pseudo-spectral wave extrapolation >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k,t;
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
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
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
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    if (!mig && src==0) {
      rick = sf_floatalloc2(nt,n_srcs);
      for (i=0; i<n_srcs; i++) {
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
    } else rick = NULL;

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    if (mig) { /* time-reversal propagation */
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
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	if (mig) {
	  /* inject data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
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
	  for (i=0; i<n_srcs; i++) {
	    for(ix=-1;ix<=1;ix++) {
	      for(iz=-1;iz<=1;iz++) {
		ik = spz[i]+iz+nz*(spx[i]+ix);
		j = spz[i]+iz+nz2*(spx[i]+ix);
		if (src==0) {
		  curr[j] += vv[ik]*rick[i][it]/(abs(ix)+abs(iz)+1);
		} else {
		  curr[j] += vv[ik]*Ricker(t, f0[i], t0[i], A[i])/(abs(ix)+abs(iz)+1);
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

	if (!mig) {
	  /* record data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
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
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    if (mig) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz)
#endif
      for (ix = 0; ix < nx1; ix++) {
	for (iz = 0; iz < nz1; iz++) {
	  img[iz + nz1*ix] = curr[iz+nbt + (ix+nbl)*nz2];
	}
      }
    }

    /*free up memory*/
    fft2_finalize();
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

int psp2(float **wvfld1, float **wvfld, float **dat, float **dat_v, float *img, float *vel, pspar par, bool adj)
/*< pseudo-spectral de-migration/migration >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*prev,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
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
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }
    if (adj) {
      for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
      }
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    if (adj) { /* RTM */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* de-migration */
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);

	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	if (adj) {
	  /* inject data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
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
	  /*adj of cross-correlation*/
	  if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	      for (iz = 0; iz < nz1; iz++) {
		i = iz + nz1*ix;
		j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
		curr[j] += wvfld1[it/snap][i]*img[i]; // !!!IMPORTANT!!! vv[ik]*
	      }
	    }
	  }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}
	
	if (adj) {
	  /* cross-correlation imaging condition */
	  if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	      for (iz = 0; iz < nz1; iz++) {
		i = iz + nz1*ix;
		j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
		img[i] += wvfld1[it/snap][i]*curr[j]; // !!!IMPORTANT!!! vv[ik]*
	      }
	    }
	  }
	} else {
	  /* record data */
	  if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	      dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
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
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
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

int psp3(float **wvfld, float **wvfld1, float **dat, float **dat1, float *img, float *vel, pspar par)
/*< pseudo-spectral back propagation of two receiver wavefields >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,*curr1,*prev,*prev1,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
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
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    curr1  = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    prev1  = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	curr1[iz] = 0.;
	prev[iz] = 0.;
	prev1[iz] = 0.;
    }
    for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* step backward in time */
    it1 = nt-1;
    it2 = -1;
    its = -1;	

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);
	
        /* first receiver wavefield */
	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

	/* inject data */
	if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	        curr[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat[ix][it];
	    }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr);
	  abc_apply(prev);
	}
	
	/* second receiver wavefield */
	/* matrix multiplication */
	fft2(curr1,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr1[j];
		c += c - prev1[j];
		prev1[j] = old;
		c += wave[j]*vv[i];
		curr1[j] = c;
	    }
	}

	/* inject data */
	if (NULL!=dat1) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
	    for (ix = 0; ix < gpl; ix++) {
	        curr1[gpz+(ix+gpx)*nz2] += vv[gpz+(ix+gpx)*nz]*dat1[ix][it];
	    }
	}

	/*apply abc*/
	if (abc) {
	  abc_apply(curr1);
	  abc_apply(prev1);
	}

	/* cross-correlation imaging condition */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	    for (ix = 0; ix < nx1; ix++) {
	        for (iz = 0; iz < nz1; iz++) {
		    i = iz + nz1*ix;
		    j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
		    img[i] += curr[j]*curr1[j]; //!!!IMPORTANT!!!
		}
	    }
	}
	
	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
              wvfld1[it/snap][i] = curr1[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(curr1);
    free(prev);
    free(prev1);
    free(cwave);
    free(cwavem);
    
    return 0;
}

int dt2v2(float **wvfld, float *vel, pspar par)
/*< inplace adjusting the background wavefield w to d^2w/dt^2/v^2 >*/
{
  int nx1, nz1, wfnt;
  int ix, iz, it, i, j;
  float wfdt;
  float *vel2;
  
  nz1 = par->nz - par->nbt - par->nbb;
  nx1 = par->nx - par->nbl - par->nbr;
  wfnt = par->nt/par->snap;
  wfdt = par->snap*par->dt;

  vel2 = sf_floatalloc(nz1*nx1);

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i,j)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = ix*nz1+iz;
      j = (ix+par->nbl)*par->nz + iz + par->nbt;
      vel2[i] = vel[j]*vel[j]*wfdt*wfdt;
    }
  }
  
  /*2nd order forward FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = iz + nz1*ix;
      wvfld[0][i] = (wvfld[2][i] - 2.f*wvfld[1][i] + wvfld[0][i])/vel2[i];
    }
  }

  for (it=1; it<wfnt-1; it++) { /*2nd order central FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
    for (ix=0; ix<nx1; ix++) {
      for (iz=0; iz<nz1; iz++) {
	i = iz + nz1*ix;
	wvfld[it][i] = (wvfld[it+1][i] - 2.f*wvfld[it][i] + wvfld[it-1][i])/vel2[i];
      }
    }
  } /*it iteration*/

  /*2nd order backward FD*/
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,i)
#endif
  for (ix=0; ix<nx1; ix++) {
    for (iz=0; iz<nz1; iz++) {
      i = iz + nz1*ix;
      wvfld[wfnt-1][i] = (wvfld[wfnt-1][i] - 2.f*wvfld[wfnt-2][i] + wvfld[wfnt-3][i])/vel2[i];
    }
  }

  return 0;
}

int psp4(float **wvfld, float **wvfld0, float **dat, float **dat_v, float *vel, pspar par)
/*< pseudo-spectral wave extrapolation using wavefield injection >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j;     /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k,t;
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
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
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
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    if (src==0) {
      rick = sf_floatalloc2(nt,n_srcs);
      for (i=0; i<n_srcs; i++) {
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
    } else rick = NULL;

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = 0.;
	prev[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* modeling */
    /* step forward in time */
    it1 = 0;
    it2 = nt;
    its = +1;

    /* MAIN LOOP */
    for (it=it1; it!=it2; it+=its) {
      
        if(verb) sf_warning("it=%d/%d;",it,nt);
	
	/* matrix multiplication */
	fft2(curr,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*lapl[ik];
#else
	  cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
	}
	
	ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		old = c = curr[j];
		c += c - prev[j];
		prev[j] = old;
		c += wave[j]*vv[i];
		curr[j] = c;
	    }
	}

        if (NULL!=wvfld0) { /* wavefield injection */
          if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,ik)
#endif
            for (ix=0; ix<nx1; ix++) {
              for (iz=0; iz<nz1; iz++) {
                i = iz + nz1*ix;
                j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                ik = iz+nbt + (ix+nbl)*nz; /* padded grid */
                curr[j] += vv[ik]*wvfld0[it/snap][i];
              }
            }
          }
        } else { /* source injection */
          t = it*dt;
          for (i=0; i<n_srcs; i++) {
            for(ix=-1;ix<=1;ix++) {
              for(iz=-1;iz<=1;iz++) {
                ik = spz[i]+iz+nz*(spx[i]+ix);
                j = spz[i]+iz+nz2*(spx[i]+ix);
                if (src==0) {
                  curr[j] += vv[ik]*rick[i][it]/(abs(ix)+abs(iz)+1);
                } else {
                  curr[j] += vv[ik]*Ricker(t, f0[i], t0[i], A[i])/(abs(ix)+abs(iz)+1);
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

        /* record data */
        if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
            for (ix = 0; ix < gpl; ix++) {
                dat[ix][it] = curr[gpz+(ix+gpx)*nz2];
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

	/* save wavefield */
	if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
	  for (ix=0; ix<nx1; ix++) {
	    for (iz=0; iz<nz1; iz++) {
	      i = iz + nz1*ix;
	      j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
	      wvfld[it/snap][i] = curr[j];
	    }
	  }
	}
    }
    if(verb) sf_warning(".");

    /*free up memory*/
    fft2_finalize();
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

int psp5(int split, float **wvfld1, float **wvfld, float **dat, float *img, float *vel, pspar par)
/*< pseudo-spectral back propagation of multiple receiver wavefields and cross-correlation >*/
{
    /*survey parameters*/
    int   nx, nz;
    float dx, dz;
    int   n_srcs;
    int   *spx, *spz;
    int   gpz, gpx, gpl;
    int   gpz_v, gpx_v, gpl_v;
    int   snap;
    /*fft related*/
    bool  cmplx;
    int   pad1;
    /*absorbing boundary*/
    bool abc;
    int nbt, nbb, nbl, nbr;
    float ct,cb,cl,cr;
    /*source parameters*/
    int src; /*source type*/
    int nt;
    float dt,*f0,*t0,*A;
    /*misc*/
    bool verb, ps;
    float vref;
    
    int nx1, nz1; /*domain of interest*/
    int it,iz,ik,ix,i,j,ir;  /* index variables */
    int nk,nzx,nz2,nx2,nzx2,nkz,nth;
    int it1, it2, its;
    float dkx,dkz,kx0,kz0,vref2,kx,kz,k;
    float c, old;
    int d_gpl, o_gpl, n_gpl;

    /*wave prop arrays*/
    float *vv;
    sf_complex *cwave,*cwavem;
    float *wave,*curr,**currs,*prev,*lapl;

    /*passing the parameters*/
    nx    = par->nx;
    nz    = par->nz;
    dx    = par->dx;
    dz    = par->dz;
    n_srcs= par->n_srcs;
    spx   = par->spx;
    spz   = par->spz;
    gpz   = par->gpz;
    gpx   = par->gpx;
    gpl   = par->gpl;
    gpz_v = par->gpz_v;
    gpx_v = par->gpx_v;
    gpl_v = par->gpl_v;
    snap  = par->snap;
    cmplx = par->cmplx;
    pad1  = par->pad1;
    abc   = par->abc;
    nbt   = par->nbt;
    nbb   = par->nbb;
    nbl   = par->nbl;
    nbr   = par->nbr;
    ct    = par->ct;
    cb    = par->cb;
    cl    = par->cl;
    cr    = par->cr;
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
    nx1 = nx-nbl-nbr;

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    
    dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
    dkx = 1./(nx2*dx); kx0 = -0.5/dx;
    nkz = (cmplx)? nz2:(nz2/2+1);
    if(nk!=nx2*nkz) sf_error("wavenumber dimension mismatch!");
    sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
    sf_warning("nk=%d,nkz=%d,nz2=%d,nx2=%d",nk,nkz,nz2,nx2);

    if(abc)
      abc_init(nz,nx,nz2,nx2,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    /* allocate and read/initialize arrays */
    vv     = sf_floatalloc(nzx); 
    lapl   = sf_floatalloc(nk);
    wave   = sf_floatalloc(nzx2);
    curr   = sf_floatalloc(nzx2);
    currs  = sf_floatalloc2(nzx2,nt/snap);
    prev   = sf_floatalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    for (iz=0; iz < nzx; iz++) {
        vv[iz] = vel[iz]*vel[iz]*dt*dt;
    }
    vref *= dt;
    vref2 = vref*vref;
    for (it=0; it<nt/snap; it++) {
        for (iz=0; iz < nzx2; iz++) {
            currs[it][iz] = 1.;
        }
    }
    for (iz=0; iz < nz1*nx1; iz++) {
        img[iz] = 0.;
    }

    /* constructing the pseudo-analytical op */
    for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nkz; iz++) {
	    kz = kz0+iz*dkz;
	    k = 2*SF_PI*hypot(kx,kz);
	    if (ps) lapl[iz+ix*nkz] = -k*k;
	    else lapl[iz+ix*nkz] = 2.*(cos(vref*k)-1.)/vref2;
	}
    }

    /* step backward in time */
    it1 = nt-1;
    it2 = -1;
    its = -1;	

    d_gpl = (int) ceil(gpl*1.0/split);
    /* Loop over groups of receivers */
    for (ir=0; ir<split; ir++) {
        o_gpl = ir*d_gpl;
        n_gpl = (ir==split-1) ? gpl-ir*d_gpl : d_gpl;

        sf_warning("ir=%d/%d",ir,split);
        sf_warning("gpl=%d,split=%d,o_gpl=%d,n_gpl=%d",gpl,split,o_gpl,n_gpl);

        for (iz=0; iz < nzx2; iz++) {
            curr[iz] = 0.;
            prev[iz] = 0.;
        }
        /* MAIN LOOP */
        for (it=it1; it!=it2; it+=its) {

            if(verb) sf_warning("it=%d/%d;",it,nt);

            /* first receiver wavefield */
            /* matrix multiplication */
            fft2(curr,cwave);

            for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                cwavem[ik] = cwave[ik]*lapl[ik];
#else
                cwavem[ik] = sf_cmul(cwave[ik],lapl[ik]);
#endif
            }

            ifft2(wave,cwavem);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j,old,c)
#endif
            for (ix = 0; ix < nx; ix++) {
                for (iz=0; iz < nz; iz++) {
                    i = iz+ix*nz;  /* original grid */
                    j = iz+ix*nz2; /* padded grid */

                    old = c = curr[j];
                    c += c - prev[j];
                    prev[j] = old;
                    c += wave[j]*vv[i];
                    curr[j] = c;
                    if (snap > 0 && it%snap==0) {
                        currs[it/snap][j] *= c;
                    }
                }
            }

            /* inject data */
            if (NULL!=dat) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix)
#endif
                for (ix = 0; ix < n_gpl; ix++) {
                    curr[gpz+(ix+o_gpl+gpx)*nz2] += vv[gpz+(ix+o_gpl+gpx)*nz]*dat[ix+o_gpl][it];
                }
            }

            /*apply abc*/
            if (abc) {
                abc_apply(curr);
                abc_apply(prev);
            }

            if (ir == split-1) {
                /* cross-correlation imaging condition */
                if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
                    for (ix = 0; ix < nx1; ix++) {
                        for (iz = 0; iz < nz1; iz++) {
                            i = iz + nz1*ix;
                            j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                            img[i] += wvfld1[it/snap][i]*currs[it/snap][j];
                        }
                    }
                }

                /* save wavefield */
                if (snap > 0 && it%snap==0) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
                    for (ix=0; ix<nx1; ix++) {
                        for (iz=0; iz<nz1; iz++) {
                            i = iz + nz1*ix;
                            j = iz+nbt + (ix+nbl)*nz2; /* padded grid */
                            wvfld[it/snap][i] = currs[it/snap][j];
                        }
                    }
                }
            }
        }
        if(verb) sf_warning(".");
    }

    /*free up memory*/
    fft2_finalize();
    if (abc) abc_close();
    free(vv);
    free(lapl);   
    free(wave);
    free(curr);
    free(*currs); free(currs);
    free(prev);
    free(cwave);
    free(cwavem);
    
    return 0;
}

