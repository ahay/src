/* Multicomponent data registration by 1-D warping. */
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

#include <string.h>
#include <math.h>

#include <rsf.h> 

int main(int argc, char* argv[])
{ 
    int i, i1, i2, m2, n1, n2, nd, dim, order, iter, nliter, niter;
    int n[SF_MAX_DIM], rect[SF_MAX_DIM], rect2[SF_MAX_DIM];
    float *coord, *inp, *out, *oth, *der, *warp;
    float *ampl=NULL, *damp=NULL, o1, d1, o2, d2, error, mean, *num, *den;
    bool verb, noamp;
    char key[6];
    sf_bands spl;
    sf_file in, warped, other, warpin, warpout, amplout=NULL;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);

    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	    snprintf(key,6,"arect%d",i+1);
	    if (!sf_getint(key,rect2+i)) rect2[i]=1;
	} else {
	    rect[i]=1;
	}
    }
    n1 = n[0];
    m2 = nd/n1;

    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    nd = m2*n2;

    sf_putint  (warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);
    
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

    warpout = sf_output("warpout");
    sf_putint(warpout,"n1",n2);
    sf_putfloat(warpout,"d1",d2);

    if (!noamp) {
	amplout = sf_output("amplout");
	sf_putint(amplout,"n1",n2);
	sf_putfloat(amplout,"d1",d2);
	sf_putfloat(amplout,"o2",1.);
    }
    
    coord = sf_floatalloc (nd); 
    inp =   sf_floatalloc (n1*m2);
    out =   sf_floatalloc (nd);
    oth =   sf_floatalloc (nd);
    der =   sf_floatalloc (nd);
    warp =  sf_floatalloc (nd);
    if (!noamp) {
	ampl =  sf_floatalloc (nd);
	damp =  sf_floatalloc (nd);
    }
    num =   sf_floatalloc (nd);
    den =   sf_floatalloc (nd);

    /* convert to spline coefficients */
    spl = sf_spline_init (order, n1);     
    for (i2=0; i2 < m2; i2++) {
	sf_floatread(inp+i2*n1,n1,in);
	sf_banded_solve(spl,inp+i2*n1);
    }
    sf_banded_close(spl);

    sf_floatread(oth,nd,other);
    sf_fileclose(other);

    if (NULL != sf_getstring ("warpin")) {
	/* optional initial warp file */
	warpin = sf_input("warpin");
	sf_floatread(coord,nd,warpin);
	sf_fileclose(warpin);
    } else {
	for (i=0; i < nd; i++) {
	    coord[i] = 0.;
	}
    }
    for (i2=0; i2 < m2; i2++) {
	for (i1=0; i1 < n2; i1++) {
	    coord[i2*n2+i1] += (o2+i1*d2);
	}
    }
    
    if (verb) sf_warning("Initialization completed");
  
    sf_divn_init(dim, nd, n, rect, niter,true);

    for (iter=0; iter < nliter; iter++) {
	for (i2=0; i2 < m2; i2++) {
	    sf_int1_init (coord+i2*n2, o1, d1, n1, sf_spline_int, order, n2, 0.0);
	    sf_int1_lop (false,false,n1,n2,inp+i2*n1,out+i2*n2);
	    
	    sf_int1_init (coord+i2*n2, o1, d1, n1, sf_spline_der, order, n2, 0.0);
	    sf_int1_lop (false,false,n1,n2,inp+i2*n1,der+i2*n2);
	}

	if (!noamp) {
	    sf_divn_close();
	    sf_divn_init(dim, nd, n, rect2, niter,true);

	    mean = 0.;
	    for (i=0; i < nd; i++) {
		mean  += out[i]*out[i];
	    }
	    mean = nd/mean;
	    
	    for (i=0; i < nd; i++) {
		num[i] = oth[i]*out[i]*mean;
		den[i] = out[i]*out[i]*mean;
	    }
	    
/*	    divlap2 (diva, m2, num, den, NULL, ampl); */
/*	    divide2 (diva, m2, num, den, ampl); */

	    sf_divn (num, den, ampl);
	    
	    for (i=0; i < nd; i++) {
		num[i] = (oth[i]-2.*ampl[i]*out[i])*der[i]*mean;
	    }

	    sf_divn (num, den, ampl);
	
	    for (i=0; i < nd; i++) {
		der[i] = ampl[i]*der[i] + out[i]*damp[i];
	    }

	    sf_divn_close();
	    sf_divn_init(dim, nd, n, rect, niter, true);
	} /* if amp */

	error = 0.;
	mean = 0.;

	for (i=0; i < nd; i++) {
	    if (noamp) {
		out[i] -= oth[i];
	    } else {
		out[i] = ampl[i]*out[i] - oth[i];
	    }
	    error += out[i]*out[i];
	    mean  += der[i]*der[i];
	}
	error = sqrt (error/nd);
	mean = nd/mean;

	if (verb) fprintf(stderr,"%d\t%f\t%f\n",iter,error,sqrt(mean));

	for (i=0; i < nd; i++) {
	    out[i] *= der[i]*mean;
	    der[i] *= der[i]*mean;
	}

	sf_divn(out, der, warp);

	for (i=0; i < nd; i++) {
	    coord[i] -= warp[i]*d2;
	}
    }

    for (i2=0; i2 < m2; i2++) {
	sf_int1_init (coord+i2*n2, o1, d1, n1, sf_spline_int, order, n2, 0.0);
	sf_int1_lop (false,false,n1,n2,inp+i2*n1,out+i2*n2);

	for (i1=0; i1 < n2; i1++) {
	    warp[i2*n2+i1] = coord[i2*n2+i1] - (o2+i1*d2);
	}
    }

    if (nliter > 0 && !noamp) {
	for (i=0; i < nd; i++) {
	    out[i] *= ampl[i];
	}

	sf_floatwrite(ampl,nd,amplout);
    }

    sf_floatwrite(out,nd,warped);
    sf_floatwrite(warp,nd,warpout);


    exit (0);
}

/* 	$Id: Mwarp1.c 13985 2015-03-26 13:56:59Z sfomel $	 */
