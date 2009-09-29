/* 1-D data regularization using freqlet transform and PEF frequencies*/
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

#include <rsf.h>

#include "freqintpef.h"
#include "freqlets.h"
#include "cweight.h"

int main(int argc, char *argv[])
{
    int nd, i1, n1, i2, n2, nw, n1w, niter, i, ncycle;
    float *w0, d1, o1, *ww, *crd;
    char *type;
    sf_complex *pp, *qq, *mm, *z0;
    sf_file in, out, w, coord;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("freq");
    coord = sf_input("coord");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    /* get data size */
    if (!sf_histint(in,"n1",&nd)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    /* specify model size */
    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* output samples */
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");   
    /* output sampling */
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");   
    /* output origin */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);
    
    if (!sf_histint(w,"n1",&nw)) sf_error("No n1= in freq");

    if (SF_FLOAT == sf_gettype(w)) {
	w0 = sf_floatalloc(nw);
	z0 = NULL;
    } else if (SF_COMPLEX == sf_gettype(w)) {
	w0 = NULL;
	z0 = sf_complexalloc(nw);
    } else {
	sf_error("Need float or complex type in freq");
	w0 = NULL;
	z0 = NULL;
    }
    n1w = n1*nw;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations for inversion */

    if (!sf_getint("ncycle",&ncycle)) ncycle=1;
    /* number of IRLS iterations */
    
    pp = sf_complexalloc(nd);
    crd = sf_floatalloc(nd);
    qq = sf_complexalloc(n1w);
    mm = sf_complexalloc(n1);

    if (ncycle > 0) {
	ww = sf_floatalloc(n1w);
	cweight_init(ww);
    } else {
	ww = NULL;
    }

    if (NULL == (type=sf_getstring("type"))) type="bior";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */
    
    sf_floatread(crd,nd,coord);

    freqint_init(nd,crd,n1,d1,o1,sf_lin_int,2,true,true,type[0],nw,w0,z0);
    
    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	/* read frequencies */
	if (NULL != w0) {
	    sf_floatread(w0,nw,w);
	} else {
	    sf_complexread(z0,nw,w);
	}

	/* read data */
	sf_complexread(pp,nd,in);

	/* apply adjoint */
	freqint_lop(true,false,n1w,nd,qq,pp);

	/* do inversion if ncycle > 0 */
	for (i=0; i < ncycle; i++) {	    
	    /* weight by absolute value */
	    for (i1=0; i1 < n1w; i1++) {
		ww[i1] = cabsf(qq[i1]);
	    }

	    /* sparse inverse */
	    sf_csolver_prec (freqint_lop,sf_ccgstep,cweight_lop,n1w,
			     n1w,nd,qq,pp,niter,0.,"verb",true,"end");
	    sf_ccgstep_close();
	} 

	/* reconstruct regular data */
	freqlets_lop(false,false,n1w,n1,qq,mm);
	sf_complexwrite(mm,n1,out);
    }

    exit(0);
}
