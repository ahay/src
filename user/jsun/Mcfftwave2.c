/* 1-D complex lowrank FFT wave extrapolation */
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

#include "fft2.h"

int main(int argc, char* argv[]) 
{
  int it,iz,im,ik,ix,i,j; //index variables
  int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
  float dt, f, fscale;
  float *curr, **wave;
  //  float *ww;
  //  float *rr;
  sf_complex **lft, **rht, *cwave, *cwavem;
  sf_file Fr;
  sf_file inp, out, left, right;
  sf_axis at,az,ax;    /* cube axes */

  sf_init(argc,argv);
  inp = sf_input("in");
  out = sf_output("out");
  Fr  = sf_input("ref");
    
  /* Read/Write axes */
  at = sf_iaxa(inp,1); nt = sf_n(at); 
  az = sf_iaxa(Fr,1); nz = sf_n(az); 
  ax = sf_iaxa(Fr,2); nx = sf_n(ax); 
  
  sf_oaxa(out,az,1); 
  sf_oaxa(out,ax,2); 
  sf_oaxa(out,at,3);
    

  //  if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
  //  if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");
  //  if (!sf_getint("nt",&nt)) sf_error("No nt= in input");

  //  if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");
  if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

  //  sf_putint(out,"n2",nt);
  //  sf_putfloat(out,"d2",dt);
  //  sf_putfloat(out,"o2",0.);
  //  sf_putstring(out,"label2","Time");
  //  sf_putstring(out,"unit2","s");


    nk = fft2_init(true,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;

    left = sf_input("left");   /* Left matrix */
    right = sf_input("right"); /* Right matrix */

    if (SF_COMPLEX != sf_gettype(left) ||
	SF_COMPLEX != sf_gettype(right)) sf_error("Need complex left and right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2)) sf_error("Need n2=%d in left",m2);

    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

    lft = sf_complexalloc2(nzx,m2);
    rht = sf_complexalloc2(m2,nk);

    sf_complexread(lft[0],nzx*m2,left);
    sf_complexread(rht[0],nk*m2,right);

    sf_fileclose(left);
    sf_fileclose(right);

    /* read wavelet & reflectivity */
    //    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    //    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);  

    curr = sf_floatalloc(nzx2);
    cwave = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    ifft2_allocate(cwavem);
 
    /* read the initial data */
    sf_floatread (curr,nx,inp);

    for (iz = 0; iz < nzx2; iz++) {
	curr[iz] = 0.;
    }

    fscale = 1.0/nx2;

    /* propagation in time */
    for (it=0; it < nt; it++) {
	/* FFT: curr -> cwave */
	fft2(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rht[ik][im]*fscale;
#else
		cwavem[ik] = sf_crmul(sf_cmul(cwave[ik],rht[ik][im]),fscale);
#endif
	    }
	    /* Inverse FFT: cwavem -> wave[im] */
	    ifft2(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	  for (iz=0; iz < nz; iz++) {
	    i = iz+ix*nz;  /* original grid */
	    j = iz+ix*nz2; /* padded grid */
	    
	    f = 0.0f;
	    for (im = 0; im < m2; im++) {
		f += lft[im][i]*wave[im][j];
	    }
	    curr[ix] = f;
	} 
	  /* write wavefield to output */
	  sf_floatwrite(curr+ix*nz2,nz,out);
	}
    }
    
    exit(0);
}
