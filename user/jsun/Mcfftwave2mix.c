/* Complex 2-D wave propagation (with kiss-fft)*/
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

#include "cfft2.h"

int main(int argc, char* argv[])
{
    bool verb,complx;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
    sf_complex c;

    float  *rr;      /* I/O arrays*/
    sf_complex *ww, *cwave, *cwavem;

    sf_complex **wave,**wave2, *curr, *currm;
    float *rcurr;
    
    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    sf_complex **lt, **rt;
    sf_file left, right;


    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("cmplx",&complx)) complx=true; /* outputs complex wavefield */

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
    at = sf_iaxa(Fw,1); nt = sf_n(at); 
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);
    
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

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
    currm  = sf_complexalloc(nzx2);
    
    if(!complx) rcurr  = sf_floatalloc(nzx2);

    cwave  = sf_complexalloc(nk);
    cwavem  = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
	if(!complx) rcurr[iz]= 0.;
    }

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* matrix multiplication */
	for (im = 0; im < m2; im++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		    currm[j] = lt[im][i]*curr[j];
#else
		    currm[j] = sf_cmul(lt[im][i], curr[j]);
#endif
		}
	    }
	    cfft2(currm,wave[im]);
	}

/* approach no.2 */
	for (ik = 0; ik < nk; ik++) {
	    c = sf_cmplx(0.,0.);
	    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		c += wave[im][ik]*rt[ik][im];
#else
		c += sf_cmul(wave[im][ik],rt[ik][im]); //complex multiplies complex
#endif
	    }
	    cwave[ik] = c;
	}

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]); //complex multiplies complex
#endif
	    }
	    icfft2(wave2[im],cwavem);
	}

/*
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwave[ik] = wave[im][ik]*rt[ik][im]*rt[ik][im];
#else
		cwave[ik] = sf_cmul(sf_cmul(wave[im][ik],rt[ik][im]),rt[ik][im]);
#endif
	    }
	    icfft2(wave2[im],cwave);
	}
*/
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		c = ww[it] * rr[i]; // source term
#else
		c = sf_crmul(ww[it], rr[i]); // source term
#endif
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave2[im][j];
#else
		    c += sf_cmul(lt[im][i], wave2[im][j]);
#endif
		}

		curr[j] = c;
		if (!complx) rcurr[j] = crealf(c);

	    }

	    /* write wavefield to output */
	    if (complx)
		sf_complexwrite(curr+ix*nz2,nz,Fo);
	    else
		sf_floatwrite(rcurr+ix*nz2,nz,Fo);
	}
    }
    if(verb) sf_warning("."); 
    exit (0);
}
