/* 1-D FFT wave extrapolation */
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

int main(int argc, char* argv[]) 
{
    int nx, nx2, nk, nt, nm, ix, ik, it, im;
    float dt, f, *curr, *prev, **lft, **rht, **mid, **wave;
    sf_complex **mat, *cwave, *cwavem;
    kiss_fftr_cfg cfg, icfg;
    sf_file inp, out, prop, left, right;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");

    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.);
    sf_putstring(out,"label2","Time");
    sf_putstring(out,"unit2","s");

    nk = kiss_fft_next_fast_size((nx+1)/2)+1;
    nx2 = 2*(nk-1);

    curr = sf_floatalloc(nx2);
    prev = sf_floatalloc(nx);

    cwave = sf_complexalloc(nk);

    cfg = kiss_fftr_alloc(nx2,0,NULL,NULL);

    prop = sf_input("prop");

    if (NULL != sf_getstring("left")) {
	if (SF_FLOAT != sf_gettype(prop)) sf_error("Need float prop");

	left = sf_input("left");
	right = sf_input("right");

	if (!sf_histint(prop,"n1",&nm)) sf_error("No n1= in prop");

	lft = sf_floatalloc2(nk,nm);
	mid = sf_floatalloc2(nm,nm);
	rht = sf_floatalloc2(nm,nx);

	sf_floatread(lft[0],nk*nm,left);
	sf_floatread(rht[0],nx*nm,right);
	sf_floatread(mid[0],nm*nm,prop);

	sf_fileclose(left);
	sf_fileclose(right);
	sf_fileclose(prop);

	mat = NULL;

	cwavem = sf_complexalloc(nk);
	wave = sf_floatalloc2(nx,nm);

	icfg = kiss_fftr_alloc(nx2,1,NULL,NULL);
    } else {
	if (SF_COMPLEX != sf_gettype(prop)) sf_error("Need complex prop");

	lft = NULL;
	rht = NULL;
	mid = NULL;

	mat = sf_complexalloc2(nk,nx);
	sf_complexread(mat[0],nx*nk,prop);
	sf_fileclose(prop);

	cwavem = NULL;
	wave = NULL;
	
	icfg = NULL;
    }
    
    sf_floatread (curr,nx,inp);
    for (ix = 0; ix < nx; ix++) {
	prev[ix] = 0;
    }
    for (ix = nx; ix < nx2; ix++) {
	curr[ix] = 0.;
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
	kiss_fftr (cfg,curr,(kiss_fft_cpx *) cwave);

	if (NULL == mat) {
	    for (im = 0; im < nm; im++) {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*lft[im][ik]/nx2;
#else
		    cwavem[ik] = sf_crmul(cwave[ik],lft[im][ik]/nx2);
#endif
		}
		kiss_fftri(icfg,(kiss_fft_cpx *) cwavem,wave[im]);
	    }
	}

	for (ix = 0; ix < nx; ix++) {
	    f = curr[ix];
	    curr[ix] = -prev[ix];
	    prev[ix] = f;

	    if (NULL == mat) {
		for (im = 0; im < nm; im++) {
		    for (ik = 0; ik < nm; ik++) {
			curr[ix] += rht[ix][ik]*mid[ik][im]*wave[im][ix];
		    }
		}
	    } else {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    curr[ix] += crealf(mat[ix][ik] * cwave[ik]);
#else
		    curr[ix] += crealf(sf_cmul(mat[ix][ik],cwave[ik]));
#endif
		}
	    }
	} 

	sf_floatwrite(curr,nx,out);
    }
    
    exit(0);
}
