/* Simple 2-D wave propagation */
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

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

#include "fft2.h"

int main(int argc, char* argv[])
{
    bool verb, cmplx;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
    float c, old;

    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem;
    float **wave, *curr, *prev;

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    float **lt, **rt;
    sf_file left, right;

#ifdef SF_HAS_FFTW
    sf_complex *cc, **cw;
    fftwf_plan fft, *ifft;
#endif

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); 
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);
    
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

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
    prev = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    for (iz=0; iz < nzx; iz++) {
	prev[iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
    }

#ifdef SF_HAS_FFTW
    if (cmplx) {
	cc = sf_complexalloc(nzx2);
	fft = fftwf_plan_dft_2d(nx2,nz2,
				(fftwf_complex *) cc, 
				(fftwf_complex *) cwave,
				FFTW_FORWARD, FFTW_MEASURE);
	if (NULL == fft) sf_error("FFTW failure.");
	cw = sf_complexalloc2(nzx2,m2);
	for (im = 0; im < m2; im++) {
	    ifft[im] = fftwf_plan_dft_2d(nx2,nz2,
					 (fftwf_complex *) cwavem, 
					 (fftwf_complex *) cw[im],
					 FFTW_BACKWARD, FFTW_MEASURE);
	    if (NULL == ifft[im]) sf_error("FFTW failure.");
	} 
    } else {
	fft = fftwf_plan_dft_r2c_2d(nx2,nz2,
				    curr, (fftwf_complex *) cwave,
				    FFTW_MEASURE);
	if (NULL == fft) sf_error("FFTW failure.");
	ifft = (fftwf_plan *) sf_alloc(m2,sizeof(fftwf_plan));
	for (im = 0; im < m2; im++) {
	    ifft[im] = fftwf_plan_dft_c2r_2d(nx2,nz2,
					     (fftwf_complex *) cwavem,
					     wave[im],
					     FFTW_MEASURE);
	    if (NULL == ifft[im]) sf_error("FFTW failure.");
	}
    }
#endif

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* matrix multiplication */
	fft2_shift(curr);

#ifdef SF_HAS_FFTW
	if (cmplx) {
	    for (ix=0; ix < nzx2; ix++) {
		cc[ix] = sf_cmplx(curr[ix],0.0f);
	    }
	}
	fftwf_execute(fft);
#else
	fft2(curr,cwave);
#endif

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_crmul(cwave[ik],rt[ik][im]);
#endif
	    }
#ifdef SF_HAS_FFTW
	    fftwf_execute(ifft[im]);
	    if (cmplx) {
		for (ix=0; ix < nzx2; ix++) {
		    wave[im][ix] = crealf(cw[im][ix]);
		}
	    }	    
#else
	    ifft2(wave[im],cwavem);	    
#endif
	    fft2_unshift(wave[im]);
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		
		old = c = curr[j];
		c += c + ww[it] * rr[i] - prev[i];
		prev[i] = old;

		for (im = 0; im < m2; im++) {
		    c += lt[im][i]*wave[im][j];
		}

		curr[j] = c;
	    }
	    	
	    /* write wavefield to output */
	    sf_floatwrite(curr+ix*nz2,nz,Fo);
	}
    }
    if(verb) sf_warning(".");    
    
    exit (0);
}
