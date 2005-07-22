/* Riemannian Wavefield Extrapolation: modeling */
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
#include "rweone.h"

int main(int argc, char* argv[])
{
    sf_file Fd, Fi, Fm, Fr;
    axa ag,at,aw,ar;
    int ig,it,   ir;

    int method;
    bool verb;
    bool mod;

    complex float **dat;
    float         **img;
    float         **aa,**bb,**mm;

    complex float **ab;
    float         **a0,**b0;

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getint("method",&method)) method=0; /* extrapolation method */
    if(! sf_getbool("mod",&mod)) mod=false;     /* mod=y modeling
						   mod=n migration */

    Fd = sf_input("in");
    Fi = sf_output("out");

    Fm = sf_input("abm");
    Fr = sf_input("abr");

    iaxa(Fd,&ag,1);       /* x='position' (can be angle) */
    iaxa(Fd,&aw,2);       /* w=frequency */

    iaxa(Fm,&at,1);       /* z=extrapolation (can be time) */
    iaxa(Fr,&ar,2);       /* a,b reference */
    if(method==0) ar.n=1; /* pure F-D */

    if(verb) {
	raxa(ag);
	raxa(at);
	raxa(aw);
	raxa(ar);
    }

    oaxa(Fi,&at,1);
    oaxa(Fi,&ag,2);
    sf_settype(Fi,SF_FLOAT);

    /* read data */
    dat = sf_complexalloc2(ag.n,aw.n);
    sf_complexread(dat[0],ag.n*aw.n,Fd);

    /* read ABM */
    aa = sf_floatalloc2  (at.n,ag.n);
    bb = sf_floatalloc2  (at.n,ag.n);
    mm = sf_floatalloc2  (at.n,ag.n);

    sf_floatread(aa[0],at.n*ag.n,Fm); /* a coef */
    sf_floatread(bb[0],at.n*ag.n,Fm); /* b coef */
    sf_floatread(mm[0],at.n*ag.n,Fm); /* mask */

    /* read ABr */
    ab = sf_complexalloc2(at.n,ar.n);
    a0 = sf_floatalloc2  (at.n,ar.n);
    b0 = sf_floatalloc2  (at.n,ar.n);

    sf_complexread(ab[0],at.n*ar.n,Fr);
    for(ir=0;ir<ar.n;ir++) {
	for(it=0;it<at.n;it++) {
	    a0[ir][it] = crealf(ab[ir][it]);
	    b0[ir][it] = cimagf(ab[ir][it]);
	}
    }

    /* allocate image */
    img = sf_floatalloc2  (at.n,ag.n);
    for(ig=0;ig<ag.n;ig++) {
	for(it=0;it<at.n;it++) {
	    img[ig][it] = 0.;
	}
    }

    /* model */
    rweone_init(ag,at,aw,ar,method);
    rweone_main(dat,img,aa,bb,mm,a0,b0);

    /* write image */
    sf_floatwrite(img[0],ag.n*at.n,Fi);
}
