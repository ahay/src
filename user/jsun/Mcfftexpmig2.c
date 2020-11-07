/* Complex 2-D exploding reflector migration (read in initial complex wavefield in depth) */
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
    bool verb;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
    float z0,x0,dx,dz,dt;
    sf_complex c;

    /* I/O arrays*/
    sf_complex *cwave, *cwavem;

    sf_complex **wave, *curr;
    
    sf_file wvfld,image;    /* I/O files */

    sf_complex **lt, **rt;
    sf_file left, right;

    sf_init(argc,argv);

    /* setup I/O files */

    wvfld = sf_input("in");
    image = sf_output("out");
    sf_settype(image,SF_COMPLEX);

    /* Read/Write axes */

    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");

    if (!sf_histint(wvfld,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(wvfld,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(wvfld,"o1",&z0)) z0=0.; 

    if (!sf_histint(wvfld,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(wvfld,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(wvfld,"o2",&x0)) x0=0.; 

    sf_putint(image,"n1",nz);
    sf_putfloat(image,"d1",dz);
    sf_putfloat(image,"o1",z0);
    //sf_putfloat(image,"o1",0.);
    sf_putstring(image,"label1","Depth");
    
    sf_putint(image,"n2",nx);
    sf_putfloat(image,"d2",dx);
    sf_putfloat(image,"o2",x0);
    sf_putstring(image,"label2","Distance");

    sf_putint(image,"n3",nt);
    sf_putfloat(image,"d3",dt);
    sf_putfloat(image,"o3",0.);
    
    /* fft prep */
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

    curr   = sf_complexalloc(nzx2);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nzx2,m2);

    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }
   
    for (ix = 0; ix < nx; ix++) {
	sf_complexread(curr+ix*nz2,nz,wvfld);
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
//		sf_warning("realcwave=%g, imagcwave=%g", crealf(cwavem[ik]),cimagf(cwavem[ik]));
	    }
	    icfft2(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		c = 0.; // initialize
#else
		c = sf_cmplx(0.,0.); // initialize
#endif
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c += sf_cmul(lt[im][i], wave[im][j]);
#endif
		}

		curr[j] = c;
	    }

	    /* write wavefield to output */
	    sf_complexwrite(curr+ix*nz2,nz,image);
	}
    }
    if(verb) sf_warning("."); 
    
    exit (0);
}
