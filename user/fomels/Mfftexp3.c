/* 3-D FFT-based zero-offset exploding reflector modeling/migration  */
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

#include "fft3.h"

int main(int argc, char* argv[])
{
    bool mig;
    int it, nt, ix, nx, iy, ny, iz, nz, nx2, ny2, nz2, nxy, nzx, nzx2;
    int im, i, j, m2, it1, it2, its, ik, n2, nk, ompchunk;
    float dt, dx, dy, dz, c, old, x0, y0;
    float *curr, *prev, ***img, **dat, **lft, **rht, **wave;
    sf_complex *cwave, *cwavem;
    sf_file data, image, left, right;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */

    if(!sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */

    if (mig) { /* migration */
	data = sf_input("in");
	image = sf_output("out");

	if (!sf_histint(data,"n1",&nx)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dx)) sf_error("No d1= in input");
	if (!sf_histfloat(data,"o1",&x0)) x0=0.; 

	if (!sf_histint(data,"n2",&ny)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dy)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&y0)) y0=0.; 

	if (!sf_histint(data,"n3",&nt)) sf_error("No n3= in input");
	if (!sf_histfloat(data,"d3",&dt)) sf_error("No d3= in input");

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* time samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* time sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putfloat(image,"o1",0.);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putfloat(image,"o2",x0);
	sf_putstring(image,"label2","Inline");

	sf_putint(image,"n3",ny);
	sf_putfloat(image,"d3",dy);
	sf_putfloat(image,"o3",y0);
	sf_putstring(image,"label3","Crossline");
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1= in input");

	if (!sf_histint(image,"n2",&nx))  sf_error("No n2= in input");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(image,"o2",&x0)) x0=0.; 	

	if (!sf_histint(image,"n3",&ny)) sf_error("No n3= in input");
	if (!sf_histfloat(image,"d3",&dy)) sf_error("No d3= in input");
	if (!sf_histfloat(image,"o3",&y0)) y0=0.; 

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	sf_putint(data,"n1",nx);
	sf_putfloat(data,"d1",dx);
	sf_putfloat(data,"o1",x0);
	sf_putstring(data,"label1","Inline");

	sf_putint(data,"n2",ny);
	sf_putfloat(data,"d2",dy);
	sf_putfloat(data,"o2",y0);
	sf_putstring(data,"label2","Crossline");

	sf_putint(data,"n3",nt);
	sf_putfloat(data,"d3",dt);
	sf_putfloat(data,"o3",0.);
	sf_putstring(data,"label3","Time");
	sf_putstring(data,"unit3","s");
    }

    nk = fft3_init(false,1,nx,ny,nz,&nx2,&ny2,&nz2);

    nxy = nx*ny;
    nzx = nz*nx*ny;
    nzx2 = nz2*nx2*ny2;

    img = sf_floatalloc3(nz,nx,ny);
    dat = sf_floatalloc2(nx,ny);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lft = sf_floatalloc2(nzx,m2);
    rht = sf_floatalloc2(m2,nk);

    sf_floatread(lft[0],nzx*m2,left);
    sf_floatread(rht[0],m2*nk,right);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    ifft3_allocate(cwavem);

    for (iz=0; iz < nzx; iz++) {
	prev[iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
    }


    if (mig) { /* migration */
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* modeling */
	sf_floatread(img[0][0],nzx,image);

	/* transpose */
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    curr[ix+nx2*(iy+iz*ny2)]=img[iy][ix][iz];
		}
	    }
	}
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }


    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nxy,data);
	} else {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    dat[iy][ix] = 0.;
		}
	    }
	}

	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		if (mig) {
		    curr[ix+nx2*iy] += dat[iy][ix];
		} else {
		    dat[iy][ix] = curr[ix+nx2*iy];
		}
	    }
	}

	/* matrix multiplication */
	fft3(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rht[ik][im];
#else
		cwavem[ik] = sf_crmul(cwave[ik],rht[ik][im]);
#endif
	    }
	    ifft3(wave[im],cwavem);
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(iy,ix,iz,i,j,im,old,c) shared(curr,prev,lft,wave)
#endif	
	for (iy=0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = ix+nx*(iy+iz*ny);   /* original grid */
		    j = ix+nx2*(iy+iz*ny2); /* padded grid */
		    
		    old = c = curr[j];
		    c += c - prev[i];
		    prev[i] = old;

		    for (im = 0; im < m2; im++) {
			c += lft[im][i]*wave[im][j];
		    }
		    
		    curr[j] = c;
		}
	    }
	}
	
	if (!mig) { /* modeling -> write out data */
	    sf_floatwrite(dat[0],nxy,data);
	}
    }
    sf_warning(".");

    if (mig) {
	/* transpose */
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    img[iy][ix][iz] = curr[ix+nx2*(iy+iz*ny2)];
		}
	    }
	}

	sf_floatwrite(img[0][0],nzx,image);
    }
    
    exit(0);
}
