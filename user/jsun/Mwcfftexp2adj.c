/* 2-D FFT-based zero-offset exploding reflector modeling/migration  */
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

#include "cfft2nsps.h"
#include "timer.h"

int main(int argc, char* argv[])
{
    bool mig,timer;
    int it, nt, ix, nx, iz, nz, nx2, nz2, nzx, nzx2, ntx, pad1, snap, n0;
    int im, i, j, m2, it1, it2, its, ik, n2, nk, nth=1;
    float dt, dx, dz,x0;
    sf_complex *curr, **img, **dat, **lft, **rht, **wave, *cwave, *cwavem, c;
    sf_complex *currm;
    sf_file data, image, left, right, snaps;
    double time=0.,t0=0.,t1=0.;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&mig)) mig=false;
    /* if n, modeling; if y, migration */
    if(! sf_getbool("timer",&timer)) timer=false;
    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */
    if(!sf_getint("n0",&n0)) n0=0;
    /* geophone surface */

    if (mig) { /* migration */
	data = sf_input("in");
	image = sf_output("out");
	sf_settype(image,SF_COMPLEX);

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&x0)) x0=0.; 

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* depth samples */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling */

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
	/* time samples */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling */

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

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    ntx = nt*nx;

    img = sf_complexalloc2(nz,nx);
    dat = sf_complexalloc2(nt,nx);
    
    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	sf_settype(snaps,SF_COMPLEX);
	sf_putint(snaps,"n1",nz);
	sf_putfloat(snaps,"d1",dz);
	sf_putfloat(snaps,"o1",0.);
	sf_putstring(snaps,"label1","Depth");
	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putfloat(snaps,"o2",x0);
	sf_putstring(snaps,"label2","Distance");
	sf_putint(snaps,"n3",nt/snap);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.);
	sf_putstring(snaps,"label3","Time");
    } else {
	snaps=NULL;
    }
    
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
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);
    if (mig) {
	currm  = sf_complexalloc(nzx2);
	icfft2_allocate(cwave);
    } else {
	cwavem = sf_complexalloc(nk);
	icfft2_allocate(cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel
{   
    nth = omp_get_num_threads();
}
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif


#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    if (mig) { /* migration */
	sf_complexread(dat[0],ntx,data);
	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;	
    } else { /* modeling */
	sf_complexread(img[0],nzx,image);

	/*exploding reflector*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[iz+ix*nz2]=img[ix][iz];
	    }
	}
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    if (timer) t0 = gtod_timer();

    if (mig) { /* migration <- read data */
	/* time stepping */
	for (it=it1; it != it2; it += its) {
	    sf_warning("it=%d;",it);
	
	    /* matrix multiplication */
	    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j) shared(currm,lft,curr)
#endif
		for (ix = 0; ix < nx; ix++) {
		    for (iz=0; iz < nz; iz++) {
			i = iz+ix*nz;  /* original grid */
			j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
			currm[j] = conjf(lft[im][i])*curr[j];
#else
			currm[j] = sf_cmul(conjf(lft[im][i]), curr[j]);
#endif
		    }
		}
		cfft2(currm,wave[im]);
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ik,im,c) shared(wave,rht,cwave)
#endif	    
	    for (ik = 0; ik < nk; ik++) {
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += wave[im][ik]*conjf(rht[ik][im]);
#else
		    c += sf_cmul(wave[im][ik],conjf(rht[ik][im])); //complex multiplies complex
#endif
		}
		cwave[ik] = c;
	    }

	    icfft2(curr,cwave);

#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* data injection */
	    for (ix=0; ix < nx; ix++) {
		curr[n0+ix*nz2] += dat[ix][it];
	    }

	    if (NULL != snaps && 0 == it%snap) {
		for (ix = 0; ix < nx; ix++) {
		    sf_complexwrite(curr+ix*nz2,nz,snaps);
		}
	    }
	}
    } else { /* modeling -> write data */

	/* time stepping */
	for (it=it1; it != it2; it += its) {
	    sf_warning("it=%d;",it);
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* record data */
	    for (ix=0; ix < nx; ix++) {
		dat[ix][it] = curr[n0+ix*nz2];
	    }
	    /* matrix multiplication */
	    cfft2(curr,cwave);
	    
	    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*rht[ik][im];
#else
		    cwavem[ik] = sf_cmul(cwave[ik],rht[ik][im]);
#endif
		}
		icfft2(wave[im],cwavem);
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j,im,c) shared(curr,lft,wave)
#endif
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
		    
		    c = sf_cmplx(0.,0.); /* initialize */
		    
		    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
			c += lft[im][i]*wave[im][j];
#else
			c += sf_cmul(lft[im][i], wave[im][j]);
#endif	    
		    }
		    curr[j] = c;
		}
	    }
	    if (NULL != snaps && 0 == it%snap) {
		for (ix = 0; ix < nx; ix++) {
		    sf_complexwrite(curr+ix*nz2,nz,snaps);
		}
	    }
	}
    }
    sf_warning(".");

    if (timer)
      {
        t1 = gtod_timer();
        time = t1-t0;
        sf_warning("Time = %lf\n",time);
      }

    if (mig) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[iz+ix*nz2];
	    }
	}

	sf_complexwrite(img[0],nzx,image);
    } else {
	sf_complexwrite(dat[0],ntx,data);
    }

    cfft2_finalize();
    exit(0);
}
