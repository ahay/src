/* Development stage-2-D complex FFT-based zero-offset exploding reflector modeling/migration  */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "cfft2w.h"

int main(int argc, char* argv[])
{
    bool adj;
    int it, nt, ix, nx, iz, nz, nx2, nz2, nzx, ntx, nzx2, pad1,n0;
    int im, i, j, m2, ik, n2, nk;
    float dt, dx, dz,x0,t0;
    sf_complex *curr, **img, **dat, **lft, **rht, **wave, *cwave, *cwavem, c;
    sf_complex *currm;
    sf_complex **wave2;
    sf_file data, image, left, right;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if n, modeling; if y, migration */

    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if(!sf_getint("n0",&n0)) n0=0;/* geophone surface */
    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");
	sf_settype(image,SF_COMPLEX);

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(data,"o1",&t0)) t0=0.; 


	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&x0)) x0=0.; 


	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* depth samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putfloat(image,"o1",0.);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putfloat(image,"o2",x0);
	sf_putstring(image,"label2","Distance");
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");
	sf_settype(data,SF_COMPLEX);

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1= in input");

	if (!sf_histint(image,"n2",&nx))  sf_error("No n2= in input");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(image,"o2",&x0)) x0=0.; 	

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	sf_putint(data,"n1",nt);
	sf_putfloat(data,"d1",dt);
	sf_putfloat(data,"o1",0.);
	sf_putstring(data,"label1","Time");
	sf_putstring(data,"unit1","s");

	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putfloat(data,"o2",x0);
	sf_putstring(data,"label2","Distance");
	}
    nk = cfft2_init(pad1,nx,nz,&nx2,&nz2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    ntx = nt*nx;
    img = sf_complexalloc2(nz,nx);
    dat = sf_complexalloc2(nt,nx);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("No n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lft = sf_complexalloc2(nzx,m2);
    rht = sf_complexalloc2(m2,nk);

    sf_complexread(lft[0],nzx*m2,left);
    sf_complexread(rht[0],m2*nk,right);

    curr = sf_complexalloc(nzx2);
    currm = sf_complexalloc(nzx2);    

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
	
	wave = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);


    if (adj) icfft2_allocate(cwave);
    else icfft2_allocate(cwavem);

    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    /* Main program */


    if (adj) { /* migration */
   
 	sf_complexread(dat[0],nt*nx,data);

    /* time stepping */
    for (it=nt-1; it>-1; it--) {
	sf_warning("it=%d;",it); 


	/* matrix multiplication */
	for (im = 0; im < m2; im++) {
	for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
			i = ix+iz*nx;  /* original grid */
			j = ix+iz*nx2; /* padded grid */
	#ifdef SF_HAS_COMPLEX_H
				currm[j] = conjf(lft[im][i])*curr[j];
	#else
				currm[j] = sf_cmul(conjf(lft[im][i]), curr[j]);
	#endif
		    }
		}
		cfft2(currm,wave[im]);
	    }
	    
	    for (ik = 0; ik < nk; ik++) {
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
	#ifdef SF_HAS_COMPLEX_H
			c += wave[im][ik]*conjf(rht[ik][im]);
	#else
			c = sf_cadd(c,sf_cmul(wave[im][ik],conjf(rht[ik][im])));
	#endif
		}
		cwave[ik] = c;
	    }
	    icfft2(curr,cwave);

	/* Inject data */
	for (ix=0; ix < nx; ix++) {
	#ifdef SF_HAS_COMPLEX_H
			curr[ix] += dat[ix][it];
	#else
			curr[ix] = sf_cadd(curr[ix],dat[ix][it]);
	#endif
	    }

    }/* End of time step */

	/* output image */
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[ix+iz*nx2];
	    }
	}

	sf_complexwrite(img[0],nzx,image);
    
/////////////////////////////////////////////////////
    } else { /* modeling */
	
	sf_complexread(img[0],nzx,image);

	/* transpose & initialize exploding refl*/
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[ix+iz*nx2]=img[ix][iz];
	    }
	}

	/* Initialize recorded data */
	for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
		dat[ix][it] =sf_cmplx(0.,0.);
	    }
	}


    /* time stepping */
    for (it=0; it < nt; it ++) {
		sf_warning("it=%d;",it);

	/* record data on the surface */
	for (ix=0; ix < nx; ix++) {
			dat[ix][it] = curr[ix];
	}

	/* matrix multiplication */
	    cfft2(curr,cwave);
	    
	    for (im = 0; im < m2; im++) {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*rht[ik][im];
#else
		    cwavem[ik] = sf_cmul(cwave[ik],rht[ik][im]);
#endif
		}
		icfft2(wave2[im],cwavem);
	    }
	    
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = ix+iz*nx;  /* original grid */
		    j = ix+iz*nx2; /* padded grid */
		    c = sf_cmplx(0.,0.);
		    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
			c += lft[im][i]*wave2[im][j];
#else
			c = sf_cadd(c,sf_cmul(lft[im][i], wave2[im][j]));
#endif
		    }
		    curr[j] = c;
		}
	    }

	}/* End of time step */
	
	sf_complexwrite(dat[0],ntx,data);

   
	}
    cfft2_finalize();
    exit(0);
}
