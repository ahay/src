/* 2-D tau-p transform.

Takes: < input.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

#include "radon.h"       
#include "ctoeplitz.h"

int main (int argc, char **argv)
{
    int nt,it;                /* number of time samples, time counter */
    int nt2;                  /* extended number of time samples */
    int nx,ix;                /* number of offsets, offset counter */ 
    int np,ip;                /* number of slopes, slope counter */
    int nw,iw;                /* number of frequencies, frequency counter */
    int nc,ic;                /* number of CMP gathers, CMP gather counter */
    bool adj;                 /* adjoint flag */
    bool inv;                 /* inversion flag */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float dw,w;               /* frequency increment, frequency */
    float eps;                /* damping for inversion */ 
    float x0;                 /* reference offset */
    float complex *dd;        /* data (CMP gather)    */
    float complex *mm;        /* model (Radon gather) */
    float complex *qq=NULL;   /* work array */
    float complex **cm, **cd; /* model and data storage */
    float *xx;                /* offset header */
    float *tt;                /* trace */
    sf_file in, out, offset;

    sf_init(argc,argv);
    in = sf_input ("in");
    offset = sf_input("offset");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, perform adjoint operation */
    if (!sf_getbool("inv",&inv)) inv=adj;
    /* if y, perform inverse operation */

    /* read input file parameters */
    if (!sf_histint  (in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    nc = sf_leftsize(in,2);

    if (adj) { 
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

	/* specify slope axis */
	if (!sf_getint  ("np",&np)) sf_error("Need np");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0");
	/* p origin (if adj=y) */

	sf_putint(out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);
    } else { /* modeling */
	if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");

	/* find number of offsets */
	nx = sf_filesize(offset);
	sf_putint(out,"n2",nx);
    }

/* determine frequency sampling (for real to complex FFT) */
    nt2 = sf_npfaro(nt,2*nt);
    nw = nt2/2+1;
    dw = 2.0*SF_PI/(nt2*dt);

    if (adj && inv && !sf_getfloat("eps",&eps)) eps=1.;
    /* smoothness for inversion (if adj=inv=y) */
    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* reference offset */

    dd = sf_complexalloc(nx);    
    mm = sf_complexalloc(np);
    xx = sf_floatalloc(nx); 

    sf_floatread (xx,nx,offset);

    if (1. != x0) { 
	for (ix=0; ix < nx; ix++) {
	    xx[ix] /= x0; /* normalize offset */
	}
    }

    radon_init (nx, np, dp, p0); /* initiliaze radon operator */

    if (adj && inv) {
	qq = sf_complexalloc(np);
	ctoeplitz_init (np); /* initialize toeplitz inversion */
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
        
	    if (adj) {
		if (inv) {
		    radon_toep (qq,eps); /* fill Toeplitz matrix */
		    ctoeplitz_solve(qq,mm); /* Toeplitz inversion */
		}
  
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

/* 	$Id: Mtaup.c,v 1.4 2004/04/19 21:51:46 fomels Exp $	 */
