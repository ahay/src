/* Random plane wave modeling. */
/*
  Copyright (C) 2000 The Board of Trustees of Stanford University
  
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
#include <time.h>

#include <rsf.h>

int main (int argc, char* argv[]) 
{
    int n1,n2, i1,i2, np,ip,gauss, type,seed;
    float rc1,h1,v1,rc2,h2,v2,pmax,ampmax, p,x,w, o1,d1,o2,d2;
    float rnd1, delay1,delay2,z1,z2, theta, phi, xloc,amp;
    sf_complex **data, refl;
    sf_file inp, out;

    sf_init(argc, argv);

    inp = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"d2",&d2)) sf_error("No d2= in input");
 
    if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"o2",&o2)) sf_error("No o2= in input");

    if (!sf_getint("np",&np)) np=1;
    if (!sf_getint("gauss",&gauss)) gauss=0;
    if (!sf_getint("type",&type)) type=1;
    /* 1 single plane layer
       2 two plane layers
       3 point diffractor */


    if (!sf_getfloat("ampmax",&ampmax)) ampmax=1.;
    if (!sf_getfloat("rc1",&rc1)) rc1=0.2;
    if (!sf_getfloat("rc2",&rc2)) rc2=0.2;
    if (!sf_getfloat("h1",&h1)) h1=200.;
    if (!sf_getfloat("h2",&h2)) h2=150.;
    if (!sf_getfloat("v1",&v1)) v1=2000.;
    if (!sf_getfloat("v2",&v2)) v2=3000.;
    if (!sf_getfloat("pmax",&pmax)) pmax=0.000332;
    if (!sf_getfloat("phi",&phi)) phi=0.1;
    if (!sf_getfloat("xloc",&xloc)) xloc=200.;

    data = sf_complexalloc2(n1,n2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    data[i2][i1] = sf_cmplx(0.,0.);
	}
    }
    
    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);


    /******************************************
      *                                       *
      * Choose random dip for each plane wave *
      *                                       *
      *****************************************/

    for (ip=0; ip < np; ip++) { /* For each plane wave */
	sf_warning("plane wave number %d;",ip+1);
	rnd1 = genrand_real1();
	rnd1 = (2. * rnd1 - 1.);
	p = pmax * rnd1;

        /*****************************************************
	 *                                                   *
	 *  Calculate delays for a one or two layered earth  *
	 *                                                   *
	 ****************************************************/

	if (1==type || 2 == type) {
	    z1 = h1 / (v1 * sqrtf(1. - p*p * v1*v1));
	    delay1 = 2. * z1 - 2. * p*p * (z1 * v1*v1);
	} else {
	    delay1 = 0.0;
	}
	if (2==type) {
	    z2 = h2 / (v2 * sqrtf(1. - p*p * v2*v2));
	    delay2 = 2. * (z1 + z2) - 2. * p*p * (z1 * v1*v1 + z2 * v2*v2);
	} else {
	    delay2 = 0.0;
	}

        /***********************************
	 *                                 *
	 * Set random frequencies in 'amp' *
	 *                                 *
	 **********************************/

	for (i1=0; i1 < n1; i1++) { /* For each frequency */
	    w = 2. * SF_PI * (o1 + i1 * d1);

	    if (gauss == 0) {
		rnd1 = genrand_real1();
		rnd1 = (2. * rnd1 - 1.);
	    } else {
		rnd1 = genrand_real1();
	    }
	    amp = rnd1  *  ampmax;

	    if (1==type) {
#ifdef SF_HAS_COMPLEX_H
		refl =  amp * (1. - rc1 * cexpf(sf_cmplx(0.,w*delay1)));
#else
		refl =  sf_crmul(
		    sf_cadd(sf_cmplx(1.,0.),
			    sf_crmul(cexpf(sf_cmplx(0.,w*delay1)),-rc1)),ampl);
#endif
	    } else if (2==type) {
#ifdef SF_HAS_COMPLEX_H
		refl =  amp * (1. - rc1 * cexpf(sf_cmplx(0.,w*delay1)) - rc2 * cexpf(sf_cmplx(0.,w*delay2)));
#else
		refl = sf_crmul(
		    sf_cadd(sf_cmplx(1.,0.),
			    sf_cadd(
				sf_crmul(cexpf(sf_cmplx(0.,w*delay1)),-rc1),
				sf_crmul(cexpf(sf_cmplx(0.,w*delay2)),-rc2))),amp);
#endif
	    } else {
		refl = sf_cmplx(0.,0.);
	    }

	    theta = asinf(v1*p);

	    for (i2=0; i2 < n2; i2++) {
		x = o2 + i2 * d2;

                /**************************************************
		 *                                                *
		 *  Calculate delays for a dipping layered earth  *
		 *                                                *
		 *************************************************/

		/*********************************************
		 *                                           *
		 *  Calculate delays for a point diffractor  *
		 *                                           *
		 ********************************************/

		if (3==type) {
		    delay1 = -v1*p * (h1*tanf(theta) + x - xloc) + h1/cosf(theta) + hypotf(x-xloc,h1);


		    /****************
		     *               *
		     *  Apply delay  *
		     *               *
		     ****************/

		    delay1 = delay1 / v1;
#ifdef SF_HAS_COMPLEX_H
		    refl =  amp * (1. - rc1 * cexpf(sf_cmplx(0.,w*delay1)));
#else
		    refl =  sf_crmul(
			sf_cadd(sf_cmplx(1.,0.),
				sf_crmul(cexpf(sf_cmplx(0.,w*delay1)),-rc1)),ampl);
#endif
		}

		/*******************
		 *                  *
		 *  Write out data  *
		 *                  *
		 *******************/
#ifdef SF_HAS_COMPLEX_H
		data[i2][i1] += cexpf(sf_cmplx(0.,w*p*x))*refl;
#else
		data[i2][i1] = sf_cadd(data[i2][i1],
				       sf_cmul(cexpf(sf_cmplx(0.,w*p*x)),refl));
#endif
	    } 
	} 
    }

    sf_complexwrite(data[0],n1*n2,out);
    exit(0);
}

