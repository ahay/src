/* 1-D seislet frame */
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

#include "freqlets.h"

int main(int argc, char *argv[])
{
    int i1, n1, i2, n2, nw, n1w, ncycle;
    bool inv, verb;
    float *w0, d1, perc;
    char *type;
    sf_complex *pp, *qq, *z0;
    sf_file in, out, w;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("freq");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

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

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getint("ncycle",&ncycle)) ncycle=0;
    /* number of iterations */
    
    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (inv) {
	n2 = sf_leftsize(in,2);
	sf_unshiftdim(in, out, 2);
	sf_putint(out,"n3",1);
    } else {
	n2 = sf_leftsize(in,1);
	sf_putint(out,"n2",nw);
	(void) sf_shiftdim(in, out, 2);
    }

    pp = sf_complexalloc(n1);
    qq = sf_complexalloc(n1w);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* sampling in the input file */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    freqlets_init(n1,d1,true,true,type[0],nw,w0,z0);
    
    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	if (NULL != w0) {
	    sf_floatread(w0,nw,w);
	} else {
	    sf_complexread(z0,nw,w);
	}

	if (inv) {
	    sf_complexread(qq,n1w,in);
	} else {
	    sf_complexread(pp,n1,in);
	} 

	freqlets_lop(!inv,false,n1w,n1,qq,pp);
	if (!inv && ncycle > 0) 
	    sf_csharpinv(freqlets_lop,1./nw,ncycle,perc,verb,n1w,n1,qq,pp); 

	if (inv) {
	    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i1] /= nw;
#else
		pp[i1] = sf_crmul(pp[i1],1.0f/nw);
#endif
	    }
	    sf_complexwrite(pp,n1,out);
	} else {
	    sf_complexwrite(qq,n1w,out);
	} 
    }
		
    exit(0);
}
