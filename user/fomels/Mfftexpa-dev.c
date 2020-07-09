/* Development stage Two-step Lowrank - 2-D FFT-based zero-offset exploding reflector modeling/migration with adjoint  */
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
#include "fft2_new.h"//download raw from github(July 9) chaing line149 out->inp

int main(int argc, char* argv[])
{
    bool adj, cmplx;
    int it, nt, ix, nx, iz, nz, nx2, nz2, ntx, nzx, nzx2, pad1;
    int im, i, j, m2, ik, n2, nk, snap;
    float dt, dx, dz, old, x0, t0, cc;
    float *curr, *prev, **img, **dat, **lft, **rht;
    sf_complex *cwave, *cwavem,**wave;
    float *currm, **wave2;

    sf_file data, image, left, right, snaps;

    sf_init(argc,argv);
    if (!sf_getbool("adj",&adj)) adj=false;
    /* if n, modeling; if y, migration */

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */


    /* I/O Setup */

    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");



	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(data,"o1",&t0)) t0=0.; 


	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&x0)) x0=0.; 


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
	sf_putstring(image,"label2","Distance");
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1= in input");

	if (!sf_histint(image,"n2",&nx))  sf_error("No n2= in input");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(image,"o2",&x0)) x0=0.; 	

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */
	if (!sf_getfloat("t0",&t0)) t0=0.0f;
	/* time origin (if modeling) */

	sf_putint(data,"n1",nt);
	sf_putfloat(data,"d1",dt);
	sf_putfloat(data,"o1",t0);
	sf_putstring(data,"label1","Time");
	sf_putstring(data,"unit1","s");

	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putfloat(data,"o2",x0);
	sf_putstring(data,"label2","Distance");


    }

    /* Wavefield snapshots setup */

    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */

    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	
	sf_putint(snaps,"n1",nz);
	sf_putfloat(snaps,"d1",dz);
	sf_putfloat(snaps,"o1",0.);
	sf_putstring(snaps,"label1","Depth");

	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putfloat(snaps,"o2",x0);
	sf_putstring(snaps,"label2","Distance");

	if(adj){
		sf_putint(snaps,"n3",nt/snap);
		sf_putfloat(snaps,"d3",-dt*snap);
		sf_putfloat(snaps,"o3",(nt-1)*dt);
	} else{
		sf_putint(snaps,"n3",nt/snap);
		sf_putfloat(snaps,"d3",dt*snap);
		sf_putfloat(snaps,"o3",t0);
	}

    } else {
	snaps = NULL;
    }

    /* Parameters and arrays setup */

    nk = fft2_init(cmplx,pad1,nx,nz,&nx2,&nz2);

    ntx = nt*nx;
    nzx = nz*nx;
    nzx2 = nz2*nx2;

    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nt,nx);

    /* propagator matrices */
    left = sf_input("left"); 
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("No n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lft = sf_floatalloc2(nzx,m2);
    rht = sf_floatalloc2(m2,nk);

    sf_floatread(lft[0],nzx*m2,left);
    sf_floatread(rht[0],m2*nk,right);

    curr = sf_floatalloc(nzx2);
    currm = sf_floatalloc(nzx2);    

    prev = sf_floatalloc(nzx);


    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

    wave = sf_complexalloc2(nk,m2);
    wave2  = sf_floatalloc2(nzx2,m2);

/*    if (adj) fft2_allocate(cwave);
    else fft2_allocate(cwavem);*/

	for (ik = 0; ik < nk; ik++) {
		for (im = 0; im < m2; im++) {
			wave[im][ik] = sf_cmplx(0.,0.);
		}
	}

    for (iz=0; iz < nzx; iz++) {
	prev[iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
    }




    /* Main program */

    if(adj){ /* migration */

    /* scheme: p(t-dt) = 2p(t) - p(t+dt) + A'p(t) + src(t) */

    float *c; /* c is used to store 2*p(t) - p(t+dt) */

    c = sf_floatalloc(nzx2);
    sf_floatread(dat[0],ntx,data);


    /* time stepping */
    for (it=nt-1;it>-1; it--) {
	
	sf_warning("it=%d;",it);



	/* update P(t) and P(t+dt) */
	for (ix = 0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
		i = ix+iz*nx;  /* original grid */
		j = ix+iz*nx2; /* padded grid */
			
		c[j] = curr[j];
		old = curr[j];
		c[j] += c[j] - prev[i];
		prev[i] = old;

		//c is 2*p(t) - p(t+dt)
	}
	}

///////////////////////////////////////////////////////////

	sf_complex *ctmp;
	ctmp = sf_complexalloc(nk);
	
	/* for forward ffts of each rank column */
	fft2_allocate(ctmp);

	/* matrix multiplication */
	for (im = 0; im < m2; im++) {
	for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
			i = ix+iz*nx;  /* original grid */
			j = ix+iz*nx2; /* padded grid */

			currm[j] = lft[im][i]*curr[j];
		    }
		}
		fft2(currm,ctmp);
		for (ik=0;ik<nk;ik++){
			/* copy data to wave*/
			wave[im][ik] =  ctmp[ik];
		}
	}

	/* for inverse fft */

	fft2_allocate(cwave);

	for (ik = 0; ik < nk; ik++) {
		cc = 0.0;
		for (im = 0; im < m2; im++) {

			cc += wave[im][ik]*rht[ik][im];
		}	 
		cwave[ik] = cc;
	}

	ifft2(curr,cwave);
///////////////////////////////////////////////////////////			
		

	for (ix = 0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
		i = ix+iz*nx;  /* original grid */
		j = ix+iz*nx2; /* padded grid */

	/* p(t-dt) = 2p(t) - p(t+dt) + A'p(t)  */
		curr[j] += c[j];
	}
	}



	/* inject data */
	for (ix=0; ix < nx; ix++) {
		curr[ix] += dat[ix][it];
	}


	/* output snap wavefield */
	if (NULL != snaps && 0 == it%snap) {
	    for (ix=0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    img[ix][iz] = curr[ix+iz*nx2];
		}
	    }

	    sf_floatwrite(img[0],nzx,snaps);
	}


	}/* End of time stepping */

	/* collect imag at time 0 */
	for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[ix+iz*nx2];
	    }
	}



	/* output image */

	sf_floatwrite(img[0],nzx,image);

/////////////////////////////////////////////////
    } else{ /* modeling */
	
	float c;

	sf_floatread(img[0],nzx,image);

	/* transpose & initialize exploding refl*/
	for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
		curr[ix+iz*nx2]=img[ix][iz];
	    }
	}
	/* Initialize recorded data */

	for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
		dat[ix][it] = 0.0f;
	    }
	}

	/* time stepping */

	for (it=0; it<nt; it++) {
		sf_warning("it=%d;",it);

	/* record data on the surface */
	for (ix=0; ix < nx; ix++) {
			dat[ix][it] = curr[ix];
	}

	/* matrix multiplication */

/////////////////My mess of fft ///////////////////////
	/* Alloc for forward (cwave) only */
	fft2_allocate(cwave);
	fft2(curr,cwave);

	/* I use separated variables to do fft / ifft */

	/* Alloc for inverse (cwavem) only */
	fft2_allocate(cwavem);
//////////////////////////////////////////////////////


	for (im = 0; im < m2; im++) {
	for (ik = 0; ik < nk; ik++) {

		cwavem[ik] = cwave[ik]*rht[ik][im];

	    }

	    ifft2(wave2[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
		i = ix+iz*nx;  /* original grid */
		j = ix+iz*nx2; /* padded grid */
			
		old = c = curr[j];
		c += c - prev[i];
		prev[i] = old;

		for (im = 0; im < m2; im++) {
			c += lft[im][i]*wave2[im][j];		    
		}
		curr[j] = c;
	}
	}


	}/* End of time stepping */

	sf_floatwrite(dat[0],nt*nx,data);

    }

    exit(0);
}
