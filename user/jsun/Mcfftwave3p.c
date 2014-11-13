/* Simple 3-D lowrank onestep wave propagation */
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

#include "cfft3w.h"

int main(int argc, char* argv[])
{
    bool verb;        
    int it,iz,im,ik,ix,iy,i,j, snap;     /* index variables */
    int nt,nz,nx,ny, m2, nk, nzx, nz2, nx2, ny2, nzx2, n2, pad1;
    float dt;
    sf_complex c;

    float  *rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem, *ww;
    sf_complex **wave, *curr;
    float *rcurr;

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax,ay;    /* cube axes */

    sf_complex **lt, **rt;
    sf_file left, right, snaps;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=true; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at); 
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 
    ay = sf_iaxa(Fr,3); ny = sf_n(ay); 

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2);
    sf_oaxa(Fo,ay,3);
    
    sf_settype(Fo,SF_FLOAT);

    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    
    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	
	sf_oaxa(snaps,az,1); 
	sf_oaxa(snaps,ax,2);
	sf_oaxa(snaps,ay,3);
	sf_oaxa(snaps,at,4);
	sf_settype(snaps,SF_FLOAT);
	sf_putint(snaps,"n4",nt/snap);
	sf_putfloat(snaps,"d4",dt*snap);
	sf_putfloat(snaps,"o4",0.);
    } else {
	snaps = NULL;
    }

    nk = cfft3_init(pad1,nz,nx,ny,&nz2,&nx2,&ny2);

    nzx = nz*nx*ny;
    nzx2 = nz2*nx2*ny2;

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2=%d in left",m2);
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_complexalloc2(nzx,m2);
    rt = sf_complexalloc2(m2,nk);

    sf_complexread(lt[0],nzx*m2,left);
    sf_complexread(rt[0],m2*nk,right);

    /* read wavelet & reflectivity */
    ww=sf_complexalloc(nt);  sf_complexread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curr = sf_complexalloc(nzx2);
    rcurr= sf_floatalloc(nzx2);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);

    icfft3_allocate(cwavem);

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=sf_cmplx(0.,0.);
	rcurr[iz]=0.;
    }


    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* matrix multiplication */
	cfft3(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    icfft3(wave[im],cwavem);
	}

	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+nz *(ix+nx *iy); /* original grid */
		    j = iz+nz2*(ix+nx2*iy); /* padded grid */
#ifdef SF_HAS_COMPLEX_H		
		    c = ww[it] * rr[i];
#else
		    c = sf_crmul(ww[it],rr[i]);
#endif

		    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
			c += lt[im][i]*wave[im][j];
#else
			c += sf_cmul(lt[im][i],wave[im][j]);
#endif
		    }
		    
		    curr[j] = c;
		    rcurr[j]= crealf(c);
		}
	    }
	}

	if (NULL != snaps && 0 == it%snap) {
	    for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
		    sf_floatwrite(rcurr+nz2*(ix+nx2*iy),nz,snaps);
		}
	    }
	}
    }
    if(verb) sf_warning(".");    
	    	
    /* write wavefield to output */
    
    for (iy = 0; iy < ny; iy++) {
	for (ix = 0; ix < nx; ix++) {
	    sf_floatwrite(rcurr+nz2*(ix+nx2*iy),nz,Fo);
	}
    }
    
    cfft3_finalize();
    exit (0);
}
