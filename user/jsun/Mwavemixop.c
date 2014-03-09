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

#include "cfft2nsps.h"

int prop(float *rr, sf_complex *ww, sf_complex **lt, sf_complex **rt, int nz, int nx, int nt, int m2, int nkzx, char *mode, int pad1, int snap, sf_complex **cc, sf_complex ***wvfld)
/*^*/
{
    bool verb=true;
    /* index variables */
    int it,iz,ix,im,ik,i,j,wfit;
    int nz2,nx2,nk,nzx2;
    sf_complex c;
    /* wavefield */
    sf_complex **wave,**wave2, *curr, *currm, *cwave, *cwavem;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    nzx2 = nz2*nx2;

    if (nk!=nkzx) sf_error("nk discrepancy!");

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    if (mode[0]=='n') {
	icfft2_allocate(cwave);
    } else {
	icfft2_allocate(cwavem);
    }

    /* initialization */
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }
    wfit = 0;

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	if (mode[0]!='p') {
	    if (mode[0]=='n') {
		/* source injection */
		for (ix = 0; ix < nx; ix++) {
		    for (iz=0; iz < nz; iz++) {
			i = iz+ix*nz;  /* original grid */
			j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
			curr[j] += ww[it] * rr[i]; // source term
#else
			curr[j] += sf_crmul(ww[it], rr[i]); // source term
#endif
		    }
		}
	    }

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
	    if (mode[0]=='n') {
		icfft2(curr,cwave);
	    }
	}

	if (mode[0]!='n') {
	    if (mode[0]=='p') {
		cfft2(curr,cwave);
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
		}
	    }
	}
	/* outout wavefield */
	if(snap>0) {
	    if(it%snap==0 && wfit<=(int)(nt-1)/snap) {
		for (ix=0; ix<nx; ix++)
		    for (iz=0; iz<nz; iz++)
			wvfld[wfit][ix][iz] = curr[iz+ix*nz2];
		wfit++;
	    }
	}
    } /* time stepping */
    if(verb) sf_warning("."); 
    /* output final result*/
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	    cc[ix][iz] = curr[iz+ix*nz2];
    
    cfft2_finalize();
    return 0;
}

int xjj1(float *rr, sf_complex *ww, sf_complex **lt, sf_complex **rt, sf_complex *alpha, sf_complex *beta, int nz, int nx, int nt, int m2, int nkzx, int pad1, int snap, sf_complex **cc, sf_complex ***wvfld)
/*^*/
{
    bool verb=true;
    /* index variables */
    int it,iz,ix,im,ik,i,j,wfit;
    int nz2,nx2,nk,nzx2;
    sf_complex c;
    /* wavefield */
    sf_complex **wave,**wave2, *curr, *currm, *cwave, *cwavem;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    nzx2 = nz2*nx2;

    if (nk!=nkzx) sf_error("nk discrepancy!");

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

    /* initialization */
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }
    wfit = 0;

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
	    }
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		currm[j] = curr[j]/alpha[i];
#else
		currm[j] = sf_cdiv(curr[j],alpha[i]);
#end
	    }
	}
	cfft2(currm,cwave);

	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	    cwavem[ik] = cwave[ik]/beta[ik];
#else
	    cwavem[ik] = sf_cdiv(cwave[ik],beta[ik]);
#endif
	}
	icfft2(curr,cwavem);
	
	/* outout wavefield */
	if(snap>0) {
	    if(it%snap==0 && wfit<=(int)(nt-1)/snap) {
		for (ix=0; ix<nx; ix++)
		    for (iz=0; iz<nz; iz++)
			wvfld[wfit][ix][iz] = curr[iz+ix*nz2];
		wfit++;
	    }
	}
    } /* time stepping */
    if(verb) sf_warning("."); 
    /* output final result*/
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	    cc[ix][iz] = curr[iz+ix*nz2];
    
    cfft2_finalize();
    return 0;
}

int prop1(sf_complex **ini, sf_complex **lt, sf_complex **rt, int nz, int nx, int nt, int m2, int nkzx, int pad1, int snap, sf_complex **cc, sf_complex ***wvfld, int offset)
/* nsps(-) | pspi(+) */
{
    bool verb=true;
    /* index variables */
    int it,iz,ix,im,ik,i,j,wfit;
    int nz2,nx2,nk,nzx2;
    sf_complex c;
    /* wavefield */
    sf_complex **wave,**wave2, *curr, *currm, *cwave, *cwavem;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    nzx2 = nz2*nx2;

    if (nk!=nkzx) sf_error("nk discrepancy!");

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

    /* initialization */
    for (ix = 0; ix < nx2; ix++) {
	for (iz=0; iz < nz2; iz++) {
	    j = iz+ix*nz2;
	    if (ix<nx && iz<nz)
		curr[j] = ini[ix][iz];
	    else 
		curr[j] = sf_cmplx(0.,0.);
	}
    }

    wfit = offset;

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* NSPS */
	
	/* matrix multiplication */
	for (im = 0; im < m2; im++) {
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
		    currm[j] = conjf(lt[im][i])*curr[j];
#else
		    currm[j] = sf_cmul(conjf(lt[im][i]), curr[j]);
#endif
		}
	    }
	    cfft2(currm,wave[im]);
	}
	
	for (ik = 0; ik < nk; ik++) {
	    c = sf_cmplx(0.,0.);
	    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		c += wave[im][ik]*conjf(rt[ik][im]);
#else
		c += sf_cmul(wave[im][ik],conjf(rt[ik][im]));
#endif
	    }
	    cwave[ik] = c;
	}
	
	/* saving a pair of FFTs */
	// icfft2(curr,cwave);
	
	/* PSPI */

	//  cfft2(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    icfft2(wave2[im],cwavem);
	}
	
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave2[im][j];
#else
		    c += sf_cmul(lt[im][i], wave2[im][j]);
#endif
		}
		curr[j] = c;
	    }
	}
	/* outout wavefield */
	if(snap>0) {
	    if(it%snap==0 && wfit-offset<=(int)(nt-1)/snap) {
		for (ix=0; ix<nx; ix++)
		    for (iz=0; iz<nz; iz++)
			wvfld[wfit][ix][iz] = curr[iz+ix*nz2];
		wfit++;
	    }
	}
    } /* time stepping */
    if(verb) sf_warning("."); 
    /* output final result*/
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	    cc[ix][iz] = curr[iz+ix*nz2];
    
    cfft2_finalize();
    return 0;
}

int prop2(sf_complex **ini, sf_complex **lt, sf_complex **rt, int nz, int nx, int nt, int m2, int nkzx, int pad1, int snap, sf_complex **cc, sf_complex ***wvfld, int offset)
/* pspi(-) + nsps(+) */
{
    bool verb=true;
    /* index variables */
    int it,iz,ix,im,ik,i,j,wfit;
    int nz2,nx2,nk,nzx2;
    sf_complex c;
    /* wavefield */
    sf_complex **wave,**wave2, *curr, *currm, *cwave, *cwavem;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    nzx2 = nz2*nx2;

    if (nk!=nkzx) sf_error("nk discrepancy!");

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

    /* initialization */
    for (ix = 0; ix < nx2; ix++) {
	for (iz=0; iz < nz2; iz++) {
	    j = iz+ix*nz2;
	    if (ix<nx && iz<nz)
		curr[j] = ini[ix][iz];
	    else 
		curr[j] = sf_cmplx(0.,0.);
	}
    }

    wfit = offset;

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);
	
	/* PSPI */

	cfft2(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	        cwavem[ik] = cwave[ik]*conjf(rt[ik][im]);
#else
		cwavem[ik] = sf_cmul(cwave[ik],conjf(rt[ik][im]));
#endif
	    }
	    icfft2(wave2[im],cwavem);
	}
	
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += conjf(lt[im][i])*wave2[im][j];
#else
		    c += sf_cmul(conjf(lt[im][i]), wave2[im][j]);
#endif
		}
		curr[j] = c;
	    }
	}

	/* NSPS */

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
	
	for (ik = 0; ik < nk; ik++) {
	    c = sf_cmplx(0.,0.);
	    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		c += wave[im][ik]*rt[ik][im];
#else
		c += sf_cmul(wave[im][ik],rt[ik][im]);
#endif
	    }
	    cwavem[ik] = c;
	}
	
	icfft2(curr,cwavem);

	/* outout wavefield */
	if(snap>0) {
	    if(it%snap==0 && wfit-offset<=(int)(nt-1)/snap) {
		for (ix=0; ix<nx; ix++)
		    for (iz=0; iz<nz; iz++)
			wvfld[wfit][ix][iz] = curr[iz+ix*nz2];
		wfit++;
	    }
	}
    } /* time stepping */
    if(verb) sf_warning("."); 
    /* output final result*/
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	    cc[ix][iz] = curr[iz+ix*nz2];
    
    cfft2_finalize();
    return 0;
}


int main(int argc, char* argv[])
{
    bool verb;
    char *mode;
    int nt,nz,nx,m2,nk,nzx,nz2,nx2,n2,pad1;
    int snap,wfnt;
    float dt,wfdt;

    float  *rr;          /* I/O arrarys */
    sf_complex *ww;
    sf_complex **cc, ***wvfld;
    
    sf_file Fw,Fr,Fo,Fs,Fa,Fb;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    sf_complex **lt, **rt, *af, *bt;
    sf_file left, right;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if (!sf_getint("snap",&snap)) snap=0;     /* interval for snapshots */

    mode=sf_getstring("mode"); /* default mode is pspi */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");
    if (SF_FLOAT != sf_gettype(Fr)) sf_error("Need float ref");

    sf_settype(Fo,SF_COMPLEX);

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt=sf_d(at);
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

    sf_setn(at,1);
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);
    if (snap > 0) {
      Fs = sf_output("snaps");	/* (optional) snapshot file */
      wfnt = (int)(nt-1)/snap + 1;
      wfdt = dt*snap;
      sf_setn(at,wfnt);
      sf_setd(at,wfdt);
      sf_oaxa(Fs,az,1); 
      sf_oaxa(Fs,ax,2); 
      sf_oaxa(Fs,at,3);
    } else Fs=NULL;

    if (mode[0]=='p') {sf_warning(">>>>> Using PSPI! <<<<< \n");}
    else if (mode[0]=='n') {sf_warning(">>>>> Using NSPS! <<<<< \n");}
    else if (mode[0]=='m') {sf_warning(">>>>> Using MIXED! <<<<< \n");}
    else { mode[0]='m'; sf_warning(">>>>> Default mode: Using MIX! <<<<< \n"); }

    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    nz2 = kiss_fft_next_fast_size(nz*pad1);
    nx2 = kiss_fft_next_fast_size(nx);
    nk = nz2*nx2; /*wavenumber*/
    nzx = nz*nx;

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

    alpha = sf_input("alpha");
    beta = sf_input("beta");

    if (!sf_histint(alpha,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in alpha",nzx);
    if (!sf_histint(beta,"n1",&n2) || n2 != nk) sf_error("Need n1=%d in beta",nk);

    af = sf_complexalloc(nzx);
    bt = sf_complexalloc(nk);

    sf_complexread(af,nzx,alpha);
    sf_complexread(bt,nk,beta);

    sf_fileclose(alpha);
    sf_fileclose(beta);

    /* read wavelet & reflectivity */
    ww=sf_complexalloc(nt);  
    sf_complexread(ww,nt ,Fw);

    rr=sf_floatalloc(nzx); 
    sf_floatread(rr,nzx,Fr);

    cc=sf_complexalloc2(nz,nx);
    if (snap>0) wvfld = sf_complexalloc3(nz,nx,wfnt);
    else wvfld = NULL;

    /* wave propagation*/
    int xjj1(rr, ww, lt, rt, af, bt, nz, nx, nt, m2, nkzx, pad1, snap, cc, wvfld);
    //    prop(rr, ww, lt, rt, nz, nx, nt, m2, nk, mode, 1, snap, cc, wvfld);
    /*
    int nt1 = nt/2;
    prop(rr, ww, lt, rt, nz, nx, nt1, m2, nk, mode, 1, snap, cc, wvfld);
    int offset = (int)(nt1-1)/snap + 1;
    sf_warning("offset=%d",offset);
    int nt2 = nt-nt1;
    sf_complex **dd;
    dd=sf_complexalloc2(nz,nx);
    prop2(cc, lt, rt, nz, nx, nt2, m2, nk, 1, snap, dd, wvfld, offset);
    */

    /* output result */
    sf_complexwrite(cc[0], nzx, Fo);
    if (snap>0)
      sf_complexwrite(wvfld[0][0], nzx*wfnt, Fs);
    exit (0);
}

