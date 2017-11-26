/* 2-D synthetic model from J.F.Claerbout. 

   October 2014 program of the month:
   http://ahay.org/blog/2014/10/08/program-of-the-month-sfsigmoid/
*/
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

#include "random.h"

int main (int argc, char* argv[])
{
    bool taper, reflectivity;
    int large, n1, n2, i1, i2, it, is;
    float *imp1, *imp2;
    float **dipper, **earth, **refl, **sig1, **sig2, **fault;
    float o1, d1, o2, d2, rand, t, frac;
    sf_file mod;

    sf_init(argc, argv);
    mod = sf_output("out");

    if (!sf_getint ("n1",&n1)) n1=400;
    /* vertical axis */
    if (!sf_getint ("n2",&n2)) n2=100;
    /* horizontal axis */
    sf_putint(mod,"n1",n1);
    sf_putint(mod,"n2",n2);
    sf_setformat(mod,"native_float");

    if (!sf_getint ("large",&large)) large=5*n1;
    /* reflectivity series */

    if (!sf_getfloat("o1",&o1)) o1=0.; sf_putfloat(mod,"o1",o1);
    if (!sf_getfloat("o2",&o2)) o2=0.; sf_putfloat(mod,"o2",o2);
    /* origin */
    
    if (!sf_getfloat("d1",&d1)) d1=0.004; sf_putfloat(mod,"d1",d1);
    if (!sf_getfloat("d2",&d2)) d2=0.032; sf_putfloat(mod,"d2",d2);
    /* sampling */

    if (!sf_getbool("reflectivity",&reflectivity)) reflectivity=true;
    /* if output reflectivity (otherwise output impedance model) */

    if (!sf_getbool("taper",&taper)) taper=true;
    /* if taper the edges */

    sf_putstring(mod,"label1","Time");    
    sf_putstring(mod,"label2","Lateral");
    sf_putstring(mod,"unit1","s");
    sf_putstring(mod,"unit2","km");

    imp1 = sf_floatalloc(large);
    imp2 = sf_floatalloc(large);

    dipper = sf_floatalloc2(n1,n2);
    earth  = sf_floatalloc2(n1,n2);
    refl   = sf_floatalloc2(n1,n2);
    sig1   = sf_floatalloc2(n1,n2);
    sig2   = sf_floatalloc2(n1,n2);
    fault  = sf_floatalloc2(n1,n2);

    random_init (19921995);
    imp1[0] = 1.;
    for (i1=1; i1 < large; i1++) {
	rand = random0();
	if (rand > 0.2) {
	    imp1[i1] = imp1[i1-1];
	} else {	
	    rand = random0();
	    imp1[i1] = 1.0 + 0.1 * rand;
	}
    }

    imp2[0] = 1.;
    for (i1=1; i1 < large; i1++) {
	rand = random0();
	if (rand > 0.3) {
	    imp2[i1] = imp2[i1-1];
	} else {	
	    rand = random0();
	    imp2[i1] = 1.0 + 0.1 * rand;
	}
    }

    for (i2= 0; i2 < n2; i2++) {
	for (i1= 0; i1 < n1; i1++) {
	    t = i1 + (i2+1.) * 0.3 * (1.0*n1)/(1.0*n2);
	    it = t;
	    frac = t - it;
	    if (it >= 0 && it+1 < large) {
		dipper[i2][i1] = imp1[it] * (1.-frac) + imp1[it+1] * frac;
	    } else {
		dipper[i2][i1] = 0.;
	    }
	}
    }

    for (i2= 0; i2 < n2; i2++) {
	for (i1= 0; i1 < n1; i1++) {
	    t = i1 + n1 *0.15 * ( 1.1 + cos( ( 3 *3.14 * (i2+1.))/n2) );
	    it = t;
	    is = it + 10;
	    frac = t - it;
	    if (it >= 0 && it+1 < large && is >= 0 && is+1 < large) {
		sig1[i2][i1] = imp2[it] * (1.-frac) + imp2[it+1] * frac;
		sig2[i2][i1] = imp2[is] * (1.-frac) + imp2[is+1] * frac;
	    } else {
		sig1[i2][i1] = 0.;
		sig2[i2][i1] = 0.;
	    }
	}
    }

    for (i2= 0; i2 < n2; i2++) {
	for (i1= 0; i1 < n1; i1++) {
	    if ((i2+1.)/n2-0.9 < -0.75*pow((i1+2.)/n1,2)) {
		fault[i2][i1]=sig1[i2][i1];
	    }	else {
		fault[i2][i1]=sig2[i2][i1];
	    }
	}
    }

    for (i2= 0; i2 < n2; i2++) {
	for (i1= 0; i1 < n1; i1++) {
	    it = i1 + 1. + (i2 + 1.) * .3 * (1.*n1)/(1.*n2);
	    if (it > n1/2) {
		earth[i2][i1] = fault[i2][i1];
	    }	else {
		earth[i2][i1] = dipper[i2][i1];
	    }
	}
    }

    if(reflectivity){
	for (i2= 0; i2 < n2; i2++) {
	    refl[i2][0] = 0.;
	    for (i1= 1; i1 < n1; i1++) {
		refl[i2][i1] = 
		    (earth[i2][i1-1] - earth[i2][i1])/
		    (earth[i2][i1-1] + earth[i2][i1]);
	    }
	}

	if (taper) {

	    for (i2= 0; i2 < 10; i2++) {
		for (i1= 0; i1 < n1; i1++) {
		    refl[i2][i1] *= (i2/10.);
		    refl[n2-i2-1][i1] *= (i2/10.);
		}
	    }

	    for (i2= 0; i2 < 5; i2++) {
		for (i1= 0; i1 < n1; i1++) {
		    refl[i2][i1] *= (i2/5.);
		    refl[n2-i2-1][i1] *= (i2/5.);
		}
	    }
	}

	sf_floatwrite (refl[0],n1*n2,mod);

    }else{
	/* output the impedance model */
	sf_floatwrite (earth[0],n1*n2,mod);
    }

    exit (0);
}

