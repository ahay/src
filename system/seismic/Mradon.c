/* High-resolution Radon transform. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>

#include <rsf.h>

#include "radon.h"       
#include "ctoeplitz.h"

int main (int argc, char **argv)
{
    bool verb, cmplx;         /* verbosity flag */
    int nx2;                  /* number of entries in the offset file */
    int nt,it;                /* number of time samples, time counter */
    int nt2;                  /* extended number of time samples */
    int nx,ix;                /* number of offsets, offset counter */ 
    int np,ip;                /* number of slopes, slope counter */
    int nw,iw;                /* number of frequencies, frequency counter */
    int ns,is;                /* number of spiking iterations, counter */
    int nc,ic;                /* number of CMP gathers, CMP gather counter */
    int niter;                /* number of iterations */
    bool adj, inv, spk, par;  /* adjoint, inversion, spiking, parabolic */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float dw,w0,w;            /* frequency increment, frequency */
    float eps,tol;            /* damping and tolerance for inversion */ 
    float x0, dx, ox;         /* reference offset, increment, origin */
    sf_complex *dd;           /* data (CMP gather)    */
    sf_complex *mm;           /* model (Radon gather) */
    sf_complex *pp=NULL;      /* preconditioned model */
    sf_complex *qq=NULL;      /* work array */
    sf_complex **cm, **cd;    /* model and data storage */
    float *xx;                /* offset header */
    float *tt;                /* trace */
    sf_complex *cc;           /* complex trace */
    float perc;               /* percentage for sharpening */
    kiss_fftr_cfg forw=NULL, invs=NULL;
    kiss_fft_cfg cfor=NULL, cinv=NULL;
    sf_file in, out, offset;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    cmplx = (SF_COMPLEX == sf_gettype(in))? true: false;

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, perform adjoint operation */
    if (!sf_getbool("inv",&inv)) inv=adj;
    /* if y, perform inverse operation */
    if (!sf_getbool("spk",&spk)) spk=inv;
    /* if y, use spiking (hi-res) inversion */
    if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity flag */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    if (adj) { 
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

	/* specify slope axis */

	if (!sf_getint  ("np",&np)) sf_error("Need np=");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
        /* p origin (if adj=y) */

	sf_putint(  out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);

    } else { /* modeling */

	if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");

	if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
	/* number of offsets (if adj=n) */
	
	sf_putint(out,"n2",nx);
    }

    nc = sf_leftsize(in,2);

    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nx2 = sf_filesize(offset);
	if (nx2 != nx && nx2 != nc*nx) sf_error("Wrong dimensions in offset");
    } else {
	offset = NULL;
	nx2 = nx;
    }
    
    if (cmplx) {
	/* determine frequency sampling (for complex to complex FFT) */
	nt2 = kiss_fft_next_fast_size(nx*2);
	nw = nt2;
	dw = 2.0*SF_PI/(nw*dt);
	w0 = -SF_PI/dt;

	cfor = kiss_fft_alloc(nw,0,NULL,NULL);
	cinv = kiss_fft_alloc(nw,1,NULL,NULL);
    } else {
	/* determine frequency sampling (for real to complex FFT) */
	nt2 = 2*kiss_fft_next_fast_size(nt);
	nw  = nt2/2+1;
	dw  = 2.0*SF_PI/(nt2*dt);
	w0  = 0.0f;
	
	forw = kiss_fftr_alloc(nt2,0,NULL,NULL);
	invs = kiss_fftr_alloc(nt2,1,NULL,NULL);
    }

    if (adj && inv) {
	if (!sf_getfloat("eps",&eps)) eps=1.;
	if (spk) {
	    if (!sf_getint("ns",&ns)) ns=1;
	    /* number of sharpening cycles */
	    if (!sf_getfloat("tol",&tol)) tol=1.e-6;
	    /* inversion tolerance */
	    if (!sf_getfloat("perc",&perc)) perc=50.0;
	    /* percentage for sharpening */
	    sf_sharpen_init(np,perc);
	}
    }

    dd = sf_complexalloc(nx);
    mm = sf_complexalloc(np);
    xx = sf_floatalloc(nx2);
    
    if (NULL != offset) {
	sf_floatread(xx,nx2,offset);
	sf_fileclose(offset);
    } else {
	if (adj) {
	    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
	    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	} else {
	    if (!sf_getfloat("ox",&ox)) sf_error("Need ox=");
	    if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
	}
	
	for (ix=0; ix < nx; ix++) {
	    xx[ix] = ox + ix*dx;
	}
	
	if (!adj) {
	    sf_putfloat(out,"o2",ox);
	    sf_putfloat(out,"d2",dx);
	}
    }

    if (!sf_getbool("parab",&par)) par=false;
    /* if y, parabolic Radon transform */
    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* reference offset */

    /* normalize offsets */
    for (ix=0; ix < nx2; ix++) {
	if (par) {
	    xx[ix] *= xx[ix]/(x0*x0); 
	} else if (1. != x0) {
	    xx[ix] /= x0;
	}
    }

    radon_init (nx, np, dp, p0); /* initiliaze radon operator */

    if (adj && inv) {
	qq = sf_complexalloc(np);
	ctoeplitz_init (np); /* initialize toeplitz inversion */
	if (spk) {
	    pp = sf_complexalloc(np);

	    /* cg2_init (np, nx, false, tol); initialize CG */
	    sf_cconjgrad_init(np, np, nx, nx, eps, tol, verb, false);
	    if (!sf_getint("niter",&niter)) niter=100;
	}
    }

    tt = cmplx? NULL: sf_floatalloc(nt2);
    cc = cmplx? sf_complexalloc(nt2): NULL;
    cm = sf_complexalloc2 (nw,np);
    cd = sf_complexalloc2 (nw,nx);

    for (ic = 0; ic < nc; ic++) { /* loop over CMPs */
	if(verb) sf_warning("i=%d of %d;",ic+1,nc);

	if (adj) {
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
		if (cmplx) {
		    sf_complexread(cc,nt,in);
		    for (it=nt; it < nt2; it++) {
			cc[it]=sf_cmplx(0.,0.);
		    }
		    
		    /* FFT to frequency */
		    kiss_fft(cfor,(kiss_fft_cpx *) cc, (kiss_fft_cpx *) cd[ix]);
		} else {
		    sf_floatread(tt,nt,in);
		    for (it=nt; it < nt2; it++) {
			tt[it]=0.;
		    }
		    
		    /* FFT to frequency */
		    kiss_fftr(forw,tt, (kiss_fft_cpx *) cd[ix]);
		}
	    }
	} else { /* modeling */
	    for (ip=0; ip < np; ip++) { /* loop over slopes */
		if (cmplx) {
		    sf_complexread(cc,nt,in);
		    for (it=nt; it < nt2; it++) {
			cc[it]=sf_cmplx(0.,0.);
		    }
		    
		    /* FFT to frequency */
		    kiss_fft(cfor,(kiss_fft_cpx *) cc, (kiss_fft_cpx *) cm[ip]);
		} else {
		    sf_floatread(tt,nt,in);
		    for (it=nt; it < nt2; it++) {
			tt[it]=0.;
		    }
		    
		    /* FFT to frequency */
		    kiss_fftr(forw,tt, (kiss_fft_cpx *) cm[ip]);
		}
	    }
	}
	
	for (iw=0; iw < nw; iw++) { /* loop over frequencies */
	    w = w0+iw*dw;

	    if (adj) {
		for (ix=0; ix < nx; ix++) { /* loop over offsets */
#ifdef SF_HAS_COMPLEX_H
		    dd[ix] = cd[ix][iw]/nt2; /* transpose */
#else
		    dd[ix] = sf_crmul(cd[ix][iw],1.0/nt2); 
#endif
		} 
	    } else {
		for (ip=0; ip < np; ip++) { /* loop over slopes */
#ifdef SF_HAS_COMPLEX_H
		    mm[ip] = cm[ip][iw]/nt2; /* transpose */
#else
		    mm[ip] = sf_crmul(cm[ip][iw],1.0/nt2); /* transpose */
#endif
		}
	    }
	    
	    if (nx2 == nx) {
		radon_set (w, xx);
	    } else {
		radon_set (w, xx+ic*nx);
	    }

	    radon_lop (adj,false,np,nx,mm,dd); /* apply Radon */
        
	    if (adj && inv) {
		radon_toep (qq,eps); /* fill Toeplitz matrix */
		ctoeplitz_solve(qq,mm); /* Toeplitz inversion */

		if (spk) {
		    for (is=0; is < ns; is++) {
			sf_csharpen(mm);
			sf_cconjgrad(NULL, radon_lop, sf_cweight_lop, 
				     pp, mm, dd, niter);
		    }
		}
	    }
  
	    if (adj) {
		for (ip=0; ip < np; ip++) { /* loop over slopes */
		    cm[ip][iw] = mm[ip]; /* transpose */
		}
	    } else {
		for (ix=0; ix < nx; ix++) { /* loop over offsets */
		    cd[ix][iw] = dd[ix]; /* transpose */
		}
	    }
	} /* loop over frequencies */

	if (adj) {
	    for (ip=0; ip < np; ip++) { /* loop over slopes */
		if (cmplx) {
		    /* FFT to time */
		    kiss_fft(cinv,(const kiss_fft_cpx *) cm[ip], (kiss_fft_cpx *) cc);
		    
		    sf_complexwrite(cc,nt,out);
		} else {
		    /* FFT to time */
		    kiss_fftri(invs,(const kiss_fft_cpx *) cm[ip], tt);
		
		    sf_floatwrite(tt,nt,out);
		}
	    }
	} else { /* modeling */
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
		if (cmplx) {
		    /* FFT to time */
		    kiss_fft(cinv,(const kiss_fft_cpx *) cd[ix], (kiss_fft_cpx *) cc);
		    
		    sf_complexwrite(cc,nt,out);
		} else {
		    /* FFT to time */
		    kiss_fftri(invs,(const kiss_fft_cpx *) cd[ix], tt);
		    
		    sf_floatwrite(tt,nt,out);
		}
	    }
	}
    } /* loop over CMPs */
    if(verb) sf_warning(".");

    exit (0);
}

/* 	$Id: Mradon.c 10624 2013-07-24 13:34:01Z sfomel $	 */

