/**************** System includes **************/
#include <string.h>
#include <math.h>

/**************** RSF includes *****************/
#include <rsf.h> 

/**************** Local includes ***************/
#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"
#include "divlap1.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, m2, n2, n, order, iter, nliter;
    float **coord, **inp, **out, **oth, **der, **warp, **ampl, **damp;
    float o1, d1, o2, d2, error, mean, eps, lam, eps2, lam2, **num, **den;
    bool verb;
    divlap divw, diva;
    sf_file in, warped, other, warpin, warpout, amplout;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m2)) m2 = 1;

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    sf_putint  (warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    n = n2*m2;

    if(!sf_getbool("verb",&verb)) verb = false;

    if(!sf_getint("accuracy",&order)) {
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    if (!sf_getint("nliter",&nliter)) nliter = 10;
    if (!sf_getfloat("eps",&eps)) eps=1.; eps = eps*eps;
    if (!sf_getfloat("lam",&lam)) lam=1.; lam = lam*lam;
    
    if (!sf_getfloat("eps2",&eps2)) eps2=10.; eps2 = eps2*eps2;
    if (!sf_getfloat("lam2",&lam2)) lam2=10.; lam2 = lam2*lam2;

    warpout = sf_output("warpout");
    sf_putint(warpout,"n1",n2);
    sf_putfloat(warpout,"d1",d2);
    sf_putint(warpout,"n2",m2);
    sf_putfloat(warpout,"o2",1.);

    amplout = sf_output("amplout");
    sf_putint(amplout,"n1",n2);
    sf_putfloat(amplout,"d1",d2);
    sf_putint(amplout,"n2",m2);
    sf_putfloat(amplout,"o2",1.);
    
    coord = sf_floatalloc2 (n2,m2); 
    inp =   sf_floatalloc2 (n1,m2);
    out =   sf_floatalloc2 (n2,m2);
    oth =   sf_floatalloc2 (n2,m2);
    der =   sf_floatalloc2 (n2,m2);
    warp =  sf_floatalloc2 (n2,m2);
    ampl =  sf_floatalloc2 (n2,m2);
    damp =  sf_floatalloc2 (n2,m2);
    num =   sf_floatalloc2 (n2,m2);
    den =   sf_floatalloc2 (n2,m2);

    prefilter_init (order, n1, order*10);     
    for (i2=0; i2 < m2; i2++) {
	sf_read(inp[i2],sizeof(float),n1,in);
	prefilter_apply (n1, inp[i2]);
    }
    prefilter_close();

    sf_read(oth[0],sizeof(float),m2*n2,other);
    sf_fileclose(other);

    if (NULL != sf_getstring ("warpin")) {
	warpin = sf_input("warpin");
	sf_read(coord[0],sizeof(float),m2*n2,warpin);
	sf_fileclose(warpin);
    } else {
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		coord[i2][i1] = 0.;
	    }
	}
    }
    for (i2=0; i2 < m2; i2++) {
	for (i1=0; i1 < n2; i1++) {
	    coord[i2][i1] += (o2+i1*d2);
	}
    }
    
    if (verb) sf_warning("Initialization completed");
  
    divw = divlap1_init(n2, eps, lam);
    diva = divlap1_init(n2, eps2, lam2);

    for (iter=0; iter < nliter; iter++) {
	for (i2=0; i2 < m2; i2++) {
	    int1_init (coord[i2], o1, d1, n1, spline_int, n2, order);
	    int1_lop (false,false,n1,n2,inp[i2],out[i2]);
	    
	    int1_init (coord[i2], o1, d1, n1, spline_der, n2, order);
	    int1_lop (false,false,n1,n2,inp[i2],der[i2]);
	}

	mean = 0.;
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		mean  += out[i2][i1]*out[i2][i1];
	    }
	}
	mean = n/mean;
	
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		num[i2][i1] = oth[i2][i1]*out[i2][i1]*mean;
		den[i2][i1] = out[i2][i1]*out[i2][i1]*mean;
	    }
	}

	divlap2 (diva, m2, num, den, NULL, ampl);
	
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		num[i2][i1] = (oth[i2][i1]*der[i2][1]
			       -2.*ampl[i2][i1]*der[i2][i1]*out[i2][i1])*mean;
	    }
	}

	divlap2 (diva, m2, num, den, NULL, damp);
	
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		der[i2][i1] = ampl[i2][i1]*der[i2][i1] 
		    + out[i2][i1]*damp[i2][i1];
	    }
	}


	error = 0.;
	mean = 0.;

	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		out[i2][i1] = ampl[i2][i1]*out[i2][i1] - oth[i2][i1];
		error += out[i2][i1]*out[i2][i1];
		mean  += der[i2][i1]*der[i2][i1];
	    }
	}
	error = sqrt (error/n);
	mean = n/mean;

	if (verb) fprintf(stderr,"%d\t%f\t%f\n",iter,error,sqrt(mean));

	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		out[i2][i1] *= der[i2][i1]*mean;
		der[i2][i1] *= der[i2][i1]*mean;
	    }
	}

	/* warp <- out/der */
	divlap2 (divw, m2, out, der, NULL, warp);
	
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < n2; i1++) {
		coord[i2][i1] -= warp[i2][i1]*d2;
	    }
	}
    }

    divlap1_close(diva);
    divlap1_close(divw);

    for (i2=0; i2 < m2; i2++) {
	int1_init (coord[i2], o1, d1, n1, spline_int, n2, order);
	int1_lop (false,false,n1,n2,inp[i2],out[i2]);

	if (nliter > 0) {
	    for (i1=0; i1 < n2; i1++) {
		out[i2][i1] *= ampl[i2][i1];
	    }
	}

	sf_write(out[i2],sizeof(float),n2,warped);

	for (i1=0; i1 < n2; i1++) {
	    warp[i2][i1] = coord[i2][i1] - (o2+i1*d2);
	}

	sf_write(warp[i2],sizeof(float),n2,warpout);

	if (nliter > 0) 
	    sf_write(ampl[i2],sizeof(float),n2,amplout);
    }

    exit (0);
}
