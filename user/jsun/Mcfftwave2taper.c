/* Complex 2-D wave propagation (with multi-threaded FFTW3)*/
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
#include "cfft2w.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb,complx,sub,os;
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
#ifdef _OPENMP
    int nth;
#endif
    sf_complex c,old;
    int snap;

    /* I/O arrays*/
    sf_complex *ww,*curr,*prev,*cwave,*cwavem,**wave,**lt, **rt;
    float *rcurr,*rr;
    
    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */
    sf_file left, right;
    sf_file snaps;

    /*for tapering*/
    float dt,dx,dz,dkx,dkz,kx0,kz0,kx,kz,ktmp,kx_trs,kz_trs,thresh;
    float *ktp;
    int taper;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("cmplx",&complx)) complx=true; /* outputs complex wavefield */
    if(!sf_getbool("os",&os)) os=true; /* one-step flag */
    if (os) {
      sf_warning("One-step wave extrapolation");
      if(!sf_getbool("sub",&sub)) sub=false; /* subtraction flag */
    } else {
      sf_warning("Two-step wave extrapolation");
      if(!sf_getbool("sub",&sub)) sub=true; /* subtraction flag */
    }
    if (!sf_getint("taper",&taper)) taper=0; /* tapering in the frequency domain */
    if (!sf_getfloat("thresh",&thresh)) thresh=0.92; /* tapering threshold */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");
    if (SF_FLOAT != sf_gettype(Fr)) sf_error("Need float ref");

    if(complx)
	sf_settype(Fo,SF_COMPLEX);
    else
	sf_settype(Fo,SF_FLOAT);

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

        if(complx)
          sf_settype(snaps,SF_COMPLEX);
        else
          sf_settype(snaps,SF_FLOAT);
    } else {
	snaps = NULL;
    }
    
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
  
    lt = sf_complexalloc2(nzx,m2);
    rt = sf_complexalloc2(m2,nk);

    sf_complexread(lt[0],nzx*m2,left);
    sf_complexread(rt[0],m2*nk,right);

    sf_fileclose(left);
    sf_fileclose(right);

    /* read wavelet & reflectivity */
    ww=sf_complexalloc(nt);  
    sf_complexread(ww,nt ,Fw);

    rr=sf_floatalloc(nzx); 
    sf_floatread(rr,nzx,Fr);

    curr   = sf_complexalloc(nzx2);
    if (!os) prev = sf_complexalloc(nzx2);
    else prev = NULL;
    if(!complx) rcurr  = sf_floatalloc(nzx2);
    else rcurr=NULL;

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nzx2,m2);

    /*icfft2_allocate(cwavem);*/

    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
	if (!os) prev[iz] = sf_cmplx(0.,0.);
	if(!complx) rcurr[iz]= 0.;
    }

    if (taper!=0) {
      dkz = 1./(nz2*dz); kz0 = -0.5/dz;
      dkx = 1./(nx2*dx); kx0 = -0.5/dx;
      kx_trs = thresh*fabs(0.5/dx);
      kz_trs = thresh*fabs(0.5/dz);
      sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
      sf_warning("nk=%d,nkz=%d,nkx=%d",nk,nz2,nx2);
      sf_warning("Applying kz tapering below %f",kz_trs);
      sf_warning("Applying kx tapering below %f",kx_trs);
      ktp = sf_floatalloc(nk);
      /* constructing the tapering op */
      for (ix=0; ix < nx2; ix++) {
	kx = kx0+ix*dkx;
	for (iz=0; iz < nz2; iz++) {
	  kz = kz0+iz*dkz;
	  ktmp = 1.;
	  if (fabs(kx) > kx_trs)
	    ktmp *= (fabs(kx)>kx_trs)? powf((fabs(kx0)-fabs(kx)+kx_trs)/kx0,2) : 1.;
	  if (fabs(kz) > kz_trs)
	    ktmp *= (fabs(kz)>kz_trs)? powf((fabs(kz0)-fabs(kz)+kz_trs)/kz0,2) : 1.;
     	  ktp[iz+ix*nz2] = ktmp;
	}
      }
    }

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* matrix multiplication */
	cfft2(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]); //complex multiplies complex
#endif
	    }
	    icfft2(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
	        i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		c = ww[it] * rr[i]; // source term
#else
		c = sf_crmul(ww[it], rr[i]); // source term
#endif
		if (sub) c += curr[j];
		if (!os) {
		  old = curr[j];
#ifdef SF_HAS_COMPLEX_H
		  c += sub? (old-prev[j]) : -prev[j];
#else
		  c = sf_cadd(c,sub? sf_csub(old,prev[j]) : sf_cneg(prev[j]));
#endif
		  prev[j] = old;
		}
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c += sf_cmul(lt[im][i], wave[im][j]);
#endif
		}

		curr[j] = c;
		if (!complx) rcurr[j] = crealf(c);

	    }

	    if (NULL != snaps && 0 == it%snap) {
              /* write wavefield snapshots */
	    if (complx)
              sf_complexwrite(curr+ix*nz2,nz,snaps);
            else
              sf_floatwrite(rcurr+ix*nz2,nz,snaps);
	    }
	}

	if (taper!=0) {
	  if (it%taper == 0) {
	    cfft2(curr,cwave);
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	      cwavem[ik] = cwave[ik]*ktp[ik];
#else
	      cwavem[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
	    }
	    icfft2(curr,cwavem);
            if (!os) {
              cfft2(prev,cwave);
              for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                cwavem[ik] = cwave[ik]*ktp[ik];
#else
                cwavem[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
              }
              icfft2(prev,cwavem);
            }
	  }
	}

    }
    if(verb) sf_warning("."); 

    /* write final wavefield to output */
    for (ix = 0; ix < nx; ix++) {
      if (complx)
        sf_complexwrite(curr+ix*nz2,nz,Fo);
      else
        sf_floatwrite(rcurr+ix*nz2,nz,Fo);
    }

    cfft2_finalize();
    exit (0);
}
