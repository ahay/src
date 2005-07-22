/* Riemannian Wavefield Extrapolation: 
   zero-offset modeling/migration 
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

#include <rsf.h>
#include "rweone.h"

int main(int argc, char* argv[])
{
    sf_file Fd, Fi, Fm, Fr;
    axa ag,at,aw,ar;
    int ig,it,iw,ir;

    int method;
    bool verb;
    bool inv;

    complex float **dat;
    float         **img;
    float         **aa,**bb,**mm;

    complex float **ab;
    float         **a0,**b0;

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getint("method",&method)) method=0; /* extrapolation method */
    if(! sf_getbool("inv",&inv)) inv=false;     /* y=modeling; n=migration */
						
    Fm = sf_input("abm");
    Fr = sf_input("abr");

    iaxa(Fm,&at,1);       /* 'extrapolation axis' (can be time) */
    iaxa(Fr,&ar,2);       /* a,b reference */
    if(method==0) ar.n=1; /* pure F-D */

    if(inv) {  /* modeling */
	Fi = sf_input ( "in");
	Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
	if (SF_FLOAT !=sf_gettype(Fi)) sf_error("Need float image");

	if (!sf_getint  ("nw",&aw.n)) sf_error ("Need nw=");
	if (!sf_getfloat("dw",&aw.d)) sf_error ("Need dw=");
	if (!sf_getfloat("ow",&aw.o)) aw.o=0.;

	iaxa(Fi,&ag,1);
	iaxa(Fi,&at,2);

	oaxa(Fd,&ag,1);
	oaxa(Fd,&aw,2);

	dat = sf_complexalloc2(ag.n,aw.n);
	img = sf_floatalloc2  (ag.n,at.n);

    } else {   /* migration */
	Fd = sf_input("in");
	Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
	if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");

	iaxa(Fd,&ag,1);       /* 'position axis' (can be angle) */
	iaxa(Fd,&aw,2);       /* frequency */

	oaxa(Fi,&ag,1);
	oaxa(Fi,&at,2);

	dat = sf_complexalloc2(ag.n,aw.n);
	img = sf_floatalloc2  (ag.n,at.n);
    }

    if(verb) {
	raxa(ag);
	raxa(at);
	raxa(aw);
	raxa(ar);
    }

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

    if(inv) { /* modeling */
	sf_floatread  (img[0],ag.n*at.n,Fi);

	for(iw=0;iw<aw.n;iw++) {
	    for(ig=0;ig<ag.n;ig++) {
		dat[iw][ig] = 0.;
	    }
	}

    } else {  /* migration */
	sf_complexread(dat[0],ag.n*aw.n,Fd);
	
	for(it=0;it<at.n;it++) {
	    for(ig=0;ig<ag.n;ig++) {
		img[it][ig] = 0.;
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* execute */
    rweone_init(ag,at,aw,ar,method);

    rweone_main(inv,dat,img,aa,bb,mm,a0,b0);
    /* execute */
    /*------------------------------------------------------------*/

    if(inv) { /* modeling */
	sf_complexwrite(dat[0],ag.n*aw.n,Fd);
    } else {  /* migration */
	sf_floatwrite  (img[0],ag.n*at.n,Fi);
    }
}
