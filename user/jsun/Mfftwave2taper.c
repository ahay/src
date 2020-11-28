/* Simple 2-D wave propagation with multi-threaded fftw3 */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
*/
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft2w.h"

int main(int argc, char* argv[])
{
    bool verb, cmplx;        
    int it,iz,im,ik,ix,i,j,snap,taper;     /* index variables */
    int nt,nz,nx,m2,nk,nzx,nz2,nx2,nzx2,n2,pad1,nkz;
    float c, old;
    float dt,dx,dz,dkx=0,dkz=0,kx0=0,kz0=0,kx,kz,ktmp,kx_trs,kz_trs,thresh;

    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem;
    float **wave, *curr, *prev;
    float *ktp=NULL;

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    float **lt, **rt;
    sf_file left, right, snaps;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if (!sf_getint("taper",&taper)) taper=0; /* tapering in the frequency domain */
    if (!sf_getfloat("thresh",&thresh)) thresh=0.92; /* tapering threshold */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fr,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    
    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */

    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	
	sf_oaxa(snaps,az,1); 
	sf_oaxa(snaps,ax,2);
	sf_oaxa(snaps,at,3);

	sf_putint(snaps,"n3",nt/snap);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.);
    } else {
	snaps = NULL;
    }

#ifdef _OPENMP
    int nth;
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

    if (taper!=0) {
      dkz = 1./(nz2*dz); kz0 = (cmplx)? -0.5/dz:0.;
      dkx = 1./(nx2*dx); kx0 = -0.5/dx;
      sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
      sf_warning("nk=%d,nkz=%d,nkx=%d",nk,nz2/2+1,nx2);
    }

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_floatalloc2(nzx,m2);
    rt = sf_floatalloc2(m2,nk);

    sf_floatread(lt[0],nzx*m2,left);
    sf_floatread(rt[0],m2*nk,right);

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx2);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    if (taper!=0) {
      nkz = (cmplx)? nz2:(nz2/2+1);
      kx_trs = thresh*fabs(0.5/dx);
      kz_trs = thresh*fabs(0.5/dz);
      sf_warning("Applying kz tapering below %f",kz_trs);
      sf_warning("Applying kx tapering below %f",kx_trs);
      ktp = sf_floatalloc(nk);
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
    }

    for (iz=0; iz < nzx2; iz++) {
      prev[iz]=0.;
      curr[iz]=0.;
    }

    ifft2_allocate(cwavem);

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

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
		c += c + ww[it] * rr[i] - prev[j];
		prev[j] = old;

		for (im = 0; im < m2; im++) {
		    c += lt[im][i]*wave[im][j];
		}

		curr[j] = c;
	    }
	    	
	    if (NULL != snaps && 0 == it%snap) {
		/* write wavefield snapshots */
		sf_floatwrite(curr+ix*nz2,nz,snaps);
	    }
	}
	
	if (taper!=0) {
	  if (it%taper == 0) {
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
	}
    }
    if(verb) sf_warning(".");   

    /* write final wavefield to output */
    for (ix = 0; ix < nx; ix++) {
	sf_floatwrite(curr+ix*nz2,nz,Fo); 
    }

    fft2_finalize();    
    exit (0);
}
