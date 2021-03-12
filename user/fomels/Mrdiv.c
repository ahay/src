/* Rough division. */
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

int main(int argc, char* argv[])
{
    bool verb;
    int id, nd, niter, niter2, iter, iter2;
    float *num, *den, *rat, *tmp, *num2, eps, perc, fact;
    double norm, norm0=1.0, norm2=1.0;
    sf_file fnum, fden, frat;

    sf_init(argc,argv);
    fnum = sf_input("in");
    fden = sf_input("den");
    frat = sf_output("out");

    if (SF_FLOAT != sf_gettype(fnum) ||
	SF_FLOAT != sf_gettype(fden)) sf_error("Need float input");

    nd = sf_filesize (fnum);
    num = sf_floatalloc(nd);
    den = sf_floatalloc(nd);
    rat = sf_floatalloc(nd);
    tmp = sf_floatalloc(nd);
    num2 = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("niter2",&niter2)) niter2=1;
    /* number of outer iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    eps *= eps;

    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (!sf_getfloat("fact",&fact)) fact=0.5;
    /* factor for sharpening */

    sf_floatread(num,nd,fnum);
    sf_floatread(den,nd,fden);

    sf_sharpen_init(nd,perc,fact);

    for (id=0; id < nd; id++) {
	rat[id] = 0.;
	num2[id] = 0.;
    }

    /* n2 = iter (n + (I - D F) n2 */
    for (iter2=0; iter2 < niter2; iter2++) {
	for (id=0; id < nd; id++) {
	    num2[id] += num[id] - den[id]*rat[id];
	    rat[id] = 0.0;
	}
	if (verb) {
	    if (iter2==0) {
		norm2=cblas_dsdot(nd,num2,1,num2,1);
		norm=1.0;
	    } else {
		norm=cblas_dsdot(nd,num2,1,num2,1)/norm2;
	    }
	    sf_warning("iter2=%d norm2=%g",iter2+1,norm);
	}

	/* F[n] = iter( r = S[B n + (I-BD) r] ) */
	for (iter=0; iter < niter; iter++) {
	    for (id=0; id < nd; id++) {
		rat[id] = (den[id]*num2[id]+eps*rat[id])/(den[id]*den[id]+eps);
	    }
	    sf_sharpen(rat);

	    sf_weight_lop(false,false,nd,nd,rat,tmp);
	    sf_weight_lop(true, false,nd,nd,rat,tmp);

	    if (verb) {
		if (iter==0) {
		    norm0=cblas_dsdot(nd,rat,1,rat,1);
		    norm=1.0;
		} else {
		    norm=cblas_dsdot(nd,rat,1,rat,1)/norm0;
		}
		sf_warning("\titer=%d norm=%g",iter+1,norm);
	    }
	}
    }

    sf_floatwrite(rat,nd,frat);

    exit(0);
}
