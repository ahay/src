/* Multicomponent data registration by 1-D warping.

Takes: < input.rsf > warped.rsf
*/

#include <string.h>
#include <math.h>

#include <rsf.h> 

#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"
#include "div2.h"
/* #include "divlap1.h" */
/* #include "divide1.h" */

int main(int argc, char* argv[])
{ 
    int i, i1, n1, i2, m2, n2, n, order, iter, nliter, niter;
    float *coord, *inp, *out, *oth, *der, *warp;
    float *ampl=NULL, *damp=NULL;
    float o1, d1, o2, d2, error, mean, eps, lam, eps2, lam2, *num, *den;
    bool verb, noamp, gauss;
/*    divlap divw, diva=NULL; */
/*    div1 divw, diva=NULL; */
    sf_file in, warped, other, warpin, warpout, amplout=NULL;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;
    m2 = sf_leftsize(in,1);

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    sf_putint  (warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    n = n2*m2;

    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getbool("noamp",&noamp)) noamp = false;
    /* if y, don't correct amplitudes */

    if(!sf_getint("accuracy",&order)) {
	/* [1-4] interpolation accuracy */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    if (!sf_getint("nliter",&nliter)) nliter = 10;
    /* number of non-linear iterations */
    if (!sf_getint("niter",&niter)) niter = 100;
    /* maximum number of linear iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.; 
    /* vertical smoothness */
    if (!sf_getfloat("lam",&lam)) lam=1.; 
    /* horizontal smoothness */

    eps = sqrtf(12.*eps*eps+1.);
    lam = sqrtf(12.*lam*lam+1.);

    if (1==m2) lam=1.;

    if (!noamp) {
	if (!sf_getfloat("eps2",&eps2)) eps2=10.; 
	/* vertical smoothness for amplitudes (if noamp=n) */
	if (!sf_getfloat("lam2",&lam2)) lam2=10.; 
	/* horizontal smoothness for amplitudes (if noamp=n) */

	eps2 = sqrtf(12.*eps2*eps2+1.);
	lam2 = sqrtf(12.*lam2*lam2+1.);

	if (1==m2) lam2=1.;
    }

    if (!sf_getbool("gauss",&gauss)) gauss=true;
    /* if y, use Gaussian smoothing; if n, triangle smoothing for shaping */

    warpout = sf_output("warpout");
    sf_putint(warpout,"n1",n2);
    sf_putfloat(warpout,"d1",d2);
    sf_putint(warpout,"n2",m2);

    if (!noamp) {
	amplout = sf_output("amplout");
	sf_putint(amplout,"n1",n2);
	sf_putfloat(amplout,"d1",d2);
	sf_putint(amplout,"n2",m2);
	sf_putfloat(amplout,"o2",1.);
    }
    
    coord = sf_floatalloc (n); 
    inp =   sf_floatalloc (n1*m2);
    out =   sf_floatalloc (n);
    oth =   sf_floatalloc (n);
    der =   sf_floatalloc (n);
    warp =  sf_floatalloc (n);
    if (!noamp) {
	ampl =  sf_floatalloc (n);
	damp =  sf_floatalloc (n);
    }
    num =   sf_floatalloc (n);
    den =   sf_floatalloc (n);

    prefilter_init (order, n1, order*10);     
    for (i2=0; i2 < m2; i2++) {
	sf_read(inp+i2*n1,sizeof(float),n1,in);
	prefilter_apply (n1, inp+i2*n1);
    }
    prefilter_close();

    sf_read(oth,sizeof(float),n,other);
    sf_fileclose(other);

    if (NULL != sf_getstring ("warpin")) {
	/* optional initial warp file */
	warpin = sf_input("warpin");
	sf_read(coord,sizeof(float),n,warpin);
	sf_fileclose(warpin);
    } else {
	for (i=0; i < n; i++) {
	    coord[i] = 0.;
	}
    }
    for (i2=0; i2 < m2; i2++) {
	for (i1=0; i1 < n2; i1++) {
	    coord[i2*n2+i1] += (o2+i1*d2);
	}
    }
    
    if (verb) sf_warning("Initialization completed");
  
/*    divw = divlap1_init(n2, eps, lam);
      if (!noamp) diva = divlap1_init(n2, eps2, lam2); */
/*    divw = divide1_init(n2, eps, lam);
      if (!noamp) diva = divide1_init(n2, eps2, lam2); */

    div2_init(n2, m2, eps, lam, niter, gauss);

    for (iter=0; iter < nliter; iter++) {
	for (i2=0; i2 < m2; i2++) {
	    int1_init (coord+i2*n2, o1, d1, n1, spline_int, order, n2);
	    int1_lop (false,false,n1,n2,inp+i2*n1,out+i2*n2);
	    
	    int1_init (coord+i2*n2, o1, d1, n1, spline_der, order, n2);
	    int1_lop (false,false,n1,n2,inp+i2*n1,der+i2*n2);
	}

	if (!noamp) {
	    if (eps2 != eps || lam2 != lam) {
		div2_close();
		div2_init(n2, m2, eps2, lam2, niter, gauss);
	    }

	    mean = 0.;
	    for (i=0; i < n; i++) {
		mean  += out[i]*out[i];
	    }
	    mean = n/mean;
	    
	    for (i=0; i < n; i++) {
		num[i] = oth[i]*out[i]*mean;
		den[i] = out[i]*out[i]*mean;
	    }
	    
/*	    divlap2 (diva, m2, num, den, NULL, ampl); */
/*	    divide2 (diva, m2, num, den, ampl); */

	    div2 (num, den, ampl);
	    
	    for (i=0; i < n; i++) {
		num[i] = (oth[i]-2.*ampl[i]*out[i])*der[i]*mean;
	    }

/*	    divlap2 (diva, m2, num, den, NULL, damp); */
/*	    divide2 (diva, m2, num, den, damp); */

	    div2 (num, den, ampl);
	
	    for (i=0; i < n; i++) {
		der[i] = ampl[i]*der[i] + out[i]*damp[i];
	    }

	    if (eps2 != eps || lam2 != lam) {
		div2_close();
		div2_init(n2, m2, eps, lam, niter, gauss);
	    }
	} /* if amp */

	error = 0.;
	mean = 0.;

	for (i=0; i < n; i++) {
	    if (noamp) {
		out[i] -= oth[i];
	    } else {
		out[i] = ampl[i]*out[i] - oth[i];
	    }
	    error += out[i]*out[i];
	    mean  += der[i]*der[i];
	}
	error = sqrt (error/n);
	mean = n/mean;

	if (verb) fprintf(stderr,"%d\t%f\t%f\n",iter,error,sqrt(mean));

	for (i=0; i < n; i++) {
	    out[i] *= der[i]*mean;
	    der[i] *= der[i]*mean;
	}

	/* warp <- out/der */
/*	divlap2 (divw, m2, out, der, NULL, warp); */
/*	divide2 (divw, m2, out, der, warp); */
	
	div2(out, der, warp);

	for (i=0; i < n; i++) {
	    coord[i] -= warp[i]*d2;
	}
    }

/*    if (!noamp) divlap1_close(diva);
      divlap1_close(divw);  */
/*    if (!noamp) divide1_close(diva);
      divide1_close(divw); */

    for (i2=0; i2 < m2; i2++) {
	int1_init (coord+i2*n2, o1, d1, n1, spline_int, order, n2);
	int1_lop (false,false,n1,n2,inp+i2*n1,out+i2*n2);

	for (i1=0; i1 < n2; i1++) {
	    warp[i2*n2+i1] = coord[i2*n2+i1] - (o2+i1*d2);
	}
    }

    if (nliter > 0 && !noamp) {
	for (i=0; i < n; i++) {
	    out[i] *= ampl[i];
	}

	sf_write(ampl,sizeof(float),n,amplout);
    }

    sf_write(out, sizeof(float),n,warped);
    sf_write(warp,sizeof(float),n,warpout);
    
    sf_close();
    exit (0);
}

/* 	$Id: Mwarp1.c,v 1.8 2004/04/13 06:02:53 fomels Exp $	 */
