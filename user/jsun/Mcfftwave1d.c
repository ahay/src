/* 1-D complex lowrank FFT wave extrapolation using complex to complex fft using initial condition*/
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "cfft1.h"
#include "fft1.h"

int main(int argc, char* argv[]) 
{
    bool sub,os,cft;
    int nx, nx2, nk, nt, m, ix, ik, it, im, n2;
    sf_complex *curr, **wave, *cwave, *cwavem, *prev, c, old;
    float dt;
    sf_complex **lft, **rht, **mat;
    sf_file Fw, Fo, right, left, prop;
    float *rcurr, **rwave, *rprev, f, rold, **rlft, **rrht;

    sf_init(argc,argv);
    Fw = sf_input("in");
    Fo = sf_output("out");

    if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");

    sf_settype(Fo,SF_FLOAT);

    /* Read/Write axes */

    if (!sf_histint(Fw,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");
    if (!sf_getbool("sub",&sub)) sub=false;
    /* if -1 is included in the matrix */
    if (!sf_getbool("os",&os)) os=true;
    if (!sf_getbool("cft",&cft)) cft=true;
    if (os) sf_warning("One-step wave extrapolation");
    else {
      sf_warning("Two-step wave extrapolation");
      if (!cft)
	sf_warning("... and real-complex FFT!");
    }

    sf_putint(Fo,"n2",nt);
    sf_putfloat(Fo,"d2",dt);
    sf_putfloat(Fo,"o2",0.);
    sf_putstring(Fo,"label2","Time");
    sf_putstring(Fo,"unit2","s");

    if (!os && !cft)
      nk = fft1_init(nx,&nx2);
    else {
      nk = cfft1_init(nx,&nx2);
      sf_warning("nk=%d\n", nk);
    }

    curr = sf_complexalloc(nx2);
    rcurr= sf_floatalloc(nx2);
    cwave = sf_complexalloc(nk);
    if (!os) {
      if (cft) {
	prev   = sf_complexalloc(nx2);
	rprev = NULL;
      } else {
	prev = NULL;
	rprev   = sf_floatalloc(nx2);
      }
    }

    if (NULL != sf_getstring("right")) {

      right = sf_input("right"); /* Right matrix */
      left  = sf_input("left");  /* Left matrix */

      if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nx);
      if (!sf_histint(left,"n2",&m)) sf_error("Need n2= in left");

      if (!sf_histint(right,"n1",&n2) || n2 != m) sf_error("Need n1=%d in right",m);
      if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right, now n2=%d",nk,n2);

      if (!os && !cft) {
	if (SF_FLOAT != sf_gettype(left) || SF_FLOAT != sf_gettype(right)) sf_error("Need float left and right");
	lft = NULL;
	rht = NULL;
	rlft = sf_floatalloc2(nx,m);
	rrht = sf_floatalloc2(m,nk);

	sf_floatread(rlft[0],nx*m,left);
	sf_floatread(rrht[0],nk*m,right);
      }	else {
	if (SF_COMPLEX != sf_gettype(left) || SF_COMPLEX != sf_gettype(right)) sf_error("Need complex left and right");
	lft = sf_complexalloc2(nx,m);
	rht = sf_complexalloc2(m,nk);
	rlft = NULL;
	rrht = NULL;

	sf_complexread(lft[0],nx*m,left);
	sf_complexread(rht[0],nk*m,right);
      }

      sf_fileclose(left);
      sf_fileclose(right);

      mat = NULL;

      cwavem = sf_complexalloc(nk);
      
      if (!os && !cft) {
	ifft1_allocate(cwavem);
	wave = NULL;
	rwave = sf_floatalloc2(nx2,m);
      } else {
	wave = sf_complexalloc2(nx2,m);
	rwave = NULL;
      }

    } else {

      prop  = sf_input("prop");  /* Propagator matrix */
      
      if (SF_COMPLEX != sf_gettype(prop)) sf_error("Need complex prop");

      lft = NULL;
      rht = NULL;
      rlft = NULL;
      rrht = NULL;

      if (!sf_histint(prop,"n1",&n2) || n2 != nk) sf_error("Need n1=%d in prop",nk);
      if (!sf_histint(prop,"n2",&n2) || n2 != nx) sf_error("Need n2=%d in prop, now n2=%d",nx,n2);

      mat = sf_complexalloc2(nk,nx);
      sf_complexread(mat[0],nx*nk,prop);
      sf_fileclose(prop);

      cwavem = NULL;
      wave = NULL;
      rwave = NULL;

    }

    sf_complexread (curr,nx,Fw);

    for (ix = nx; ix < nx2; ix++) {
	curr[ix] = sf_cmplx(0.,0.);
    }

    for (ix = 0; ix < nx2; ix++) {
      rcurr[ix]=crealf(curr[ix]);
      if (!os) {
	if (cft) {
	  prev[ix] = sf_cmplx(0.,0.);
	} else {
	  rprev[ix] = 0.0f;
	}
      }
    }

    if (!os && !cft) {
      for (it=0; it < nt; it++) {

	fft1(rcurr,cwave);

	if (NULL == mat) {
	  for (im = 0; im < m; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	      cwavem[ik] = cwave[ik]*rrht[ik][im];
#else
	      cwavem[ik] = sf_crmul(cwave[ik],rrht[ik][im]);
#endif
	    }
	    ifft1(rwave[im],cwavem);
	  }
	}
	

	for (ix = 0; ix < nx; ix++) {
	  rold = rcurr[ix];
	  f = sub? 2*rold: 0.0f;
	  f -= rprev[ix];
	  rprev[ix] = rold;

	  if (NULL == mat) {
	    for (im = 0; im < m; im++) {
	      f += rlft[im][ix]*rwave[im][ix];
	    }
	  } else {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	      f += crealf(mat[ix][ik] * cwave[ik]);
#else
	      f += crealf(sf_cmul(mat[ix][ik],cwave[ik]));
#endif
	    }
	  }

	  rcurr[ix] = f;
	} 

	sf_floatwrite(rcurr,nx,Fo);
      }
      
    } else {
      /* propagation in time */
      for (it=0; it < nt; it++) {
	sf_warning("it=%d;",it);

        /* FFT: curr -> cwave */
	cfft1(curr,cwave);

	if (NULL == mat) {

	  for (im = 0; im < m; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	      cwavem[ik] = cwave[ik]*rht[ik][im];
#else
	      cwavem[ik] = sf_cmul(cwave[ik],rht[ik][im]);
#endif
	    }
	    /* Inverse FFT: cwavem -> wave[im] */
	    icfft1(wave[im],cwavem);
	  }

	}

	for (ix = 0; ix < nx; ix++) {
	  
	  if (os) {
	    c = sub? curr[ix]:sf_cmplx(0.,0.);
	  } else {
	    old = curr[ix];
#ifdef SF_HAS_COMPLEX_H
	    c = sub? 2*old : sf_cmplx(0.,0.);
	    c -= prev[ix];
#else
	    c = sub? sf_crmul(old,2) : sf_cmplx(0.,0.);
	    c = sf_csub(c,prev[ix]);
#endif
	    prev[ix] = old;
	  }

	  if (NULL == mat) {
	    for (im = 0; im < m; im++) {
#ifdef SF_HAS_COMPLEX_H
	      c += lft[im][ix]*wave[im][ix];
#else
	      c = sf_cadd(c,sf_cmul(lft[im][ix], wave[im][ix]));
#endif
	    }
	  } else {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	      c += mat[ix][ik]*cwave[ik];
#else
	      c = sf_cadd(c,sf_cmul(mat[ix][ik],cwave[ik]));
#endif
	    }
	  }
	  curr[ix] = c;
	  rcurr[ix] = crealf(c);
	}

	sf_floatwrite(rcurr,nx,Fo);
      }
    }
    exit(0);
}
