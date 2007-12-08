/* 1-D freqlet transform */
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

#include "freqlet.h"

int main(int argc, char *argv[])
{
    int i1, n1, i2, n2, iw, nw, niter;
    bool inv, adj, unit;
    float *w0, d1;
    char *type;
    sf_complex *pp, *qq;
    sf_file in, out, w;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("freq");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (SF_FLOAT != sf_gettype(w)) sf_error("Need float freq");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    
    if (!sf_histint(w,"n1",&nw)) sf_error("No n1= in freq");
    w0 = sf_floatalloc(nw);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, do adjoint transform */

    if (!sf_getbool("unit",&unit)) unit=false;
    /* if y, use unitary scaling */

    if (!sf_getbool("niter",&niter)) niter=0;
    /* number of iterations for inversion */

    if (adj) {
	n2 = sf_leftsize(in,1);
	sf_putint(out,"n2",nw);
	(void) sf_shiftdim(in, out, 2);
    } else {
	n2 = sf_leftsize(in,2);
	sf_unshiftdim(in, out, 2);
    } 

    pp = sf_complexalloc(n1);
    qq = sf_complexalloc(n1);
    
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* sampling in the input file */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    freqlet_init(n1,inv,unit,type[0]);
    
    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(w0,nw,w);

	if (adj) {
	    sf_complexread(pp,n1,in);
	} else {
	    for (i1=0; i1 < n1; i1++) {
		pp[i1] = sf_cmplx(0.,0.);
	    }
	} 

	/* loop over frequencies */
	for (iw=0; iw < nw; iw++) {
	    freqlet_set(w0[iw]* 2*SF_PI*d1);
	    
	    if (adj) {
		freqlet_lop(true,false,n1,n1,qq,pp);
		sf_complexwrite(qq,n1,out);
	    } else {
		sf_complexread(qq,n1,in);		
		freqlet_lop(false,true,n1,n1,qq,pp);
	    } 
	}
	
	if (!adj) {
	    if (inv) {
		for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    pp[i1] /= nw;
#else
		    pp[i1] = sf_crmul(pp[i1],1.0f/nw);
#endif
		}
	    } 
	    sf_complexwrite(pp,n1,out);
	}
    }

    exit(0);
}
