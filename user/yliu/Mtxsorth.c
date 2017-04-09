/* Streaming orthogonalize signal and noise in t-x domain. */
/*
  Copyright (C) 2017 Jilin University
  
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
#include <math.h>

int main(int argc, char* argv[])
{
    int m[SF_MAX_DIM],i,it,ix,n1,n2,mdim,nd;
    float eps1, eps2, *w, *s, *n, remove,*sig2,*noi2;
    sf_file in,noise,fnoi2,fsig2;
    bool verb;

    sf_init(argc,argv);

    in = sf_input("in");
    noise = sf_input("noise");

    fnoi2 = sf_output("out");
    fsig2 = sf_output("sig2"); 

    if (!sf_getfloat("eps1",&eps1)) sf_error("Need eps1=");
    /* regularization 1*/
    eps1 *= eps1;
    if (!sf_getfloat("eps2",&eps2)) sf_error("Need eps2=");
    /* regularization 2*/
    eps2 *= eps2;

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    mdim = sf_filedims(in,m);
    nd = 1;
    for (i=0; i < mdim; i++) {
	nd *= m[i];
    }
    n1=m[0];
    n2=m[1];

    s = sf_floatalloc(nd);
    n = sf_floatalloc(nd);
    w = sf_floatalloc(nd);
    noi2 = sf_floatalloc(nd);
    sig2 = sf_floatalloc(nd);

    sf_floatread(s,nd,in);
    sf_floatread(n,nd,noise);

    for (i=0; i < nd; i++) {
	noi2[i] = n[i];
	sig2[i] = s[i];
    }

    w[0]=0.0f;
    for(it=0;it < n1; it++){
	for(ix=0;ix < n2; ix++){
	    if(ix-1<0 || it-1<0) {
		w[ix*n1+it]=0.;
	    } else {
		w[ix*n1+it]=(s[ix*n1+it]*n[ix*n1+it]+
			     eps1*w[(ix-1)*n1+it]+eps2*w[ix*n1+it-1])/
		    (s[ix*n1+it]*s[ix*n1+it]+eps1+eps2);
	    }
	}
    }
    
    for (i=0; i < nd; i++) {
	remove = w[i]*sig2[i];
	noi2[i] -= remove;
	sig2[i] += remove;
    }
    
    sf_floatwrite(noi2,nd,fnoi2);
    sf_floatwrite(sig2,nd,fsig2);

    exit(0);
}
/* 	$Id$	 */

