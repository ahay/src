/* High-resolution Radon transform. 

Takes: < in.rsf > out.rsf
*/

#include <math.h>

#include <rsf.h>

/* #include "cg2.h" */
#include "cweight.h"
#include "radon.h"       
#include "ctoeplitz.h"

int main (int argc, char **argv)
{
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
    float dw,w;               /* frequency increment, frequency */
    float eps,tol;            /* damping and tolerance for inversion */ 
    float x0, dx, ox;         /* reference offset, increment, origin */
    float maxw;               /* maximum weight */
    float complex *dd;        /* data (CMP gather)    */
    float complex *mm;        /* model (Radon gather) */
    float complex *pp;        /* preconditioned model */
    float complex *qq;        /* work array */
    float complex **cm, **cd; /* model and data storage */
    float *xx;                /* offset header */
    float *ww;                /* weight */
    float *tt;                /* trace */
    sf_file in, out, offset;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
    } else {
	offset = NULL;
    }

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, perform adjoint operation */
    if (!sf_getbool("inv",&inv)) inv=adj;
    /* if y, perform inverse operation */
    if (!sf_getbool("spk",&spk)) spk=inv;
    /* if y, use spiking (hi-res) inversion */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    if (adj) { 
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

	/* specify slope axis */
	if (!sf_getint("np",&np)) sf_error("Need np=");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
	/* p origin (if adj=y) */

	sf_putint(out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);
    } else { /* modeling */
	if (!sf_histint(in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");

	/* find number of offsets */
	if (NULL != offset) {
	    if (!sf_histint(offset,"n1",&nx)) sf_error ("No n1= in offset");
	} else {
	    if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
	}

	sf_putint(out,"n2",nx);
    }

    nc = sf_leftsize(in,2);

    /* determine frequency sampling (for real to complex FFT) */
    nt2 = sf_npfaro(nt,2*nt);
    nw = nt2/2+1;
    dw = 2.0*SF_PI/(nt2*dt);

    if (adj && inv) {
	if (!sf_getfloat("eps",&eps)) eps=1.;
	if (spk) {
	    if (!sf_getint("ns",&ns)) ns=1;
	    if (!sf_getfloat("tol",&tol)) tol=1.e-5;
	}
    }

    dd = sf_complexalloc(nx);
    mm = sf_complexalloc(np);
    xx = sf_floatalloc(nx);

    if (NULL != offset) {
	sf_floatread(xx,nx,offset);
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
    for (ix=0; ix < nx; ix++) {
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
	    ww = sf_floatalloc(np);
	    pp = sf_complexalloc(np);

	    /* cg2_init (np, nx, false, tol); initialize CG */
	    sf_cconjgrad_init(np, np, nx, nx, eps, tol, false, false);
	    cweight_init(ww);

	    if (!sf_getint("niter",&niter)) niter=100;
	}
    }

    tt = sf_floatalloc (nt2);
    cm = sf_complexalloc2 (nw,np);
    cd = sf_complexalloc2 (nw,nx);

    for (ic = 0; ic < nc; ic++) { /* loop over CMPs */
	if (adj) {
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
		sf_floatread(tt,nt,in);
		for (it=nt; it < nt2; it++) {
		    tt[it]=0.;
		}
		sf_pfarc(-1,nt2,tt,cd[ix]); /* FFT to frequency */
	    }
	} else { /* modeling */
	    for (ip=0; ip < np; ip++) { /* loop over slopes */
		sf_floatread(tt,nt,in);
		for (it=nt; it < nt2; it++) {
		    tt[it]=0.;
		}
		sf_pfarc(-1,nt2,tt,cm[ip]); /* FFT to frequency */
	    }
	}
	
	for (iw=0; iw < nw; iw++) { /* loop over frequencies */
	    w = iw*dw;

	    if (adj) {
		for (ix=0; ix < nx; ix++) { /* loop over offsets */
		    dd[ix] = cd[ix][iw]/nt2; /* transpose */
		} 
	    } else {
		for (ip=0; ip < np; ip++) { /* loop over slopes */
		    mm[ip] = cm[ip][iw]/nt2; /* transpose */
		}
	    }
	    
	    radon_set (w, xx);
	    radon_lop (adj,false,np,nx,mm,dd); /* apply Radon */
        
	    if (adj && inv) {
		radon_toep (qq,eps); /* fill Toeplitz matrix */
		ctoeplitz_solve(qq,mm); /* Toeplitz inversion */

		if (spk) {
		    for (is=0; is < ns; is++) {
			maxw = 0.;
			for (ip=0; ip < np; ip++) { 
			    ww[ip] = cabsf(mm[ip]); /* fill weight */
			    if (ww[ip] > maxw) maxw = ww[ip];
			}
			if (maxw == 0.) break;
			for (ip=0; ip < np; ip++) { 
			    ww[ip] /= maxw; /* normalize weight */
			}

			sf_cconjgrad(NULL, radon_lop, cweight_lop, 
				     pp, mm, dd, niter);

			/* cg (radon_lop, ww, mm, dd, nx, eps); CG inversion */
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
		sf_pfacr(1,nt2,cm[ip],tt); /* FFT to time */
		
		sf_floatwrite(tt,nt,out);
	    }
	} else { /* modeling */
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
		sf_pfacr(1,nt2,cd[ix],tt); /* FFT to time */
		
		sf_floatwrite(tt,nt,out);
	    }
	}
    } /* loop over CMPs */

    sf_close();
    exit (0);
}

/* 	$Id: Mradon.c,v 1.1 2004/05/13 22:27:10 fomels Exp $	 */
