/* Riemannian Wavefield Extrapolation: zero-offset modeling/migration */
/*
  Copyright (C) 2006 Colorado School of Mines
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
    sf_file Fd=NULL, Fi=NULL, Fm=NULL, Fr=NULL;
    sf_axis ag,at,aw,ar;
    int ig,it,iw,ir;
    int ng,nt,nw,nr;
    float ow,dw;

    int method;
    bool  verb;
    bool   adj;

    sf_complex **dat=NULL;
    float         **img=NULL;
    sf_complex  *wfl=NULL;
    float         **aa=NULL,**bb=NULL,**mm=NULL;

    sf_complex **ab=NULL;
    float         **a0=NULL,**b0=NULL;

    float w;
    char *met="";

    sf_init(argc,argv);

    if(! sf_getbool("verb", &verb))     verb=false;
    if(! sf_getint("method",&method)) method=0;    /* extrapolation method */
    if(! sf_getbool("adj",  &adj))       adj=false;/* y=modeling; n=migration */
    Fm = sf_input("abm");
    Fr = sf_input("abr");

    at = sf_iaxa(Fm,2); sf_setlabel(at,"t"); 
    ar = sf_iaxa(Fr,1); sf_setlabel(ar,"r"); /* a,b reference */
    if(method==0) sf_setn(ar,1); /* pure F-D */
    nr = sf_n(ar);

    if(adj) {  /* modeling */
	Fi = sf_input ( "in");
	Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
	if (SF_FLOAT !=sf_gettype(Fi)) sf_error("Need float image");

	if (!sf_getint  ("nw",&nw)) sf_error ("Need nw=");
	if (!sf_getfloat("dw",&dw)) sf_error ("Need dw=");
	if (!sf_getfloat("ow",&ow)) ow=0.;
	aw = sf_maxa(nw,ow,dw); sf_setlabel(aw,"w");

	ag = sf_iaxa(Fi,1); sf_setlabel(ag,"g");
	at = sf_iaxa(Fi,2); sf_setlabel(at,"t");

	sf_oaxa(Fd,ag,1);
	sf_oaxa(Fd,at,2);
	sf_oaxa(Fd,aw,3);
    } else {   /* migration */
	Fd = sf_input("in");
	Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
	if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");

        /* 'position axis' (can be angle) */
	ag = sf_iaxa(Fd,1); sf_setlabel(ag,"g"); 
        /* 'extrapolation axis' (can be time) */
	at = sf_iaxa(Fd,2); sf_setlabel(at,"t");
	/* frequency axis */
	aw = sf_iaxa(Fd,3); sf_setlabel(aw,"w");
	nw = sf_n(aw); ow = sf_o(aw); dw = sf_d(aw);

	sf_oaxa(Fi,ag,1);
	sf_oaxa(Fi,at,2);
	sf_putint(Fi,"n3",1);
    }
    ng = sf_n(ag);
    nt = sf_n(at);

    img = sf_floatalloc2  (ng,nt);
    dat = sf_complexalloc2(ng,nt);
    wfl = sf_complexalloc (ng);

    if(verb) {
	sf_raxa(ag);
	sf_raxa(at);
	sf_raxa(aw);
	sf_raxa(ar);
    }

    /* read ABM */
    aa = sf_floatalloc2  (ng,nt);
    bb = sf_floatalloc2  (ng,nt);
    mm = sf_floatalloc2  (ng,nt);

    sf_floatread(aa[0],ng*nt,Fm); /* a coef */
    sf_floatread(bb[0],ng*nt,Fm); /* b coef */
    sf_floatread(mm[0],ng*nt,Fm); /* mask */

    /* read ABr */
    ab = sf_complexalloc2(nr,nt);
    a0 = sf_floatalloc2  (nr,nt);
    b0 = sf_floatalloc2  (nr,nt);

    sf_complexread(ab[0],nr*nt,Fr);
    for(it=0;it<nt;it++) {
	for(ir=0;ir<nr;ir++) {
	    a0[it][ir] = crealf(ab[it][ir]);
	    b0[it][ir] = cimagf(ab[it][ir]);
	}
    }

/*------------------------------------------------------------*/
    switch (method) {
	case 3: met="PSC"; break;
	case 2: met="FFD"; break;
	case 1: met="SSF"; break;
	case 0: met="XFD"; break;
    }

    /* from hertz to radian */
    dw *= 2.*SF_PI;
    ow *= 2.*SF_PI;

    rweone_init(ag,at,aw,ar,method,verb);
    switch(method) {
	case 3: rweone_psc_coef(aa,bb,a0,b0); break;
	case 2: rweone_ffd_coef(aa,bb,a0,b0); break;
	case 1: ; /* SSF */                   break;
	case 0: rweone_xfd_coef(aa,bb);       break;
    }

    if(adj) { /* modeling */
	for(it=0;it<nt;it++) {
	    for(ig=0;ig<ng;ig++) { 
		dat[it][ig] = sf_cmplx(0.,0.);
	    }
	}
    } else {  /* migration */
	for(it=0;it<nt;it++) {
	    for(ig=0;ig<ng;ig++) {
		img[it][ig] = 0.;
	    }
	}
    }

/*------------------------------------------------------------*/
    if( adj) sf_floatread  (img[0],ng*nt,Fi);

    for(iw=0;iw<nw;iw++) {
	w=ow+iw*dw;
	sf_warning("%s %d %d;",met,iw,nw);
	
	if(adj) {  /* modeling */
	    w*=-2; /*      causal, two-way time */

	    for(ig=0;ig<ng;ig++) {
		wfl[ig] = sf_cmplx(0.,0.);
	    }

	    it=nt-1; rweone_zoi(adj,wfl,img[it]);
	    for(it=nt-2;it>=0;it--) {
		if(method!=0) rweone_fk(w,wfl,aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(w,wfl,aa[it],it);
		rweone_zoi(adj,wfl,img[it]);

		for(ig=0;ig<ng;ig++) {
		    dat[it][ig] = wfl[ig];
		}
	    }

	    sf_complexwrite(dat[0],ng*nt,Fd);

	} else {   /* migration */
	    w*=+2; /* anti-causal, two-way time */

	    sf_complexread(dat[0],ng*nt,Fd);

	    for(ig=0;ig<ng;ig++) {
		wfl[ig] = sf_cmplx(0.,0.);
	    }

	    for(it=0;it<=nt-2;it++) {
		for(ig=0;ig<ng;ig++) {
#ifdef SF_HAS_COMPLEX_H
		    wfl[ig] += dat[it][ig];
#else
		    wfl[ig] = sf_cadd(wfl[ig],dat[it][ig]);
#endif
		}

		rweone_zoi(adj,wfl,img[it]);		
		if(method!=0) rweone_fk(w,wfl,aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(w,wfl,aa[it],it);
	    }
	    it=nt-1; rweone_zoi(adj,wfl,img[it]);
	}
    }
    sf_warning(".");

    if(!adj) sf_floatwrite  (img[0],ng*nt,Fi);
/*------------------------------------------------------------*/

    exit(0);
}
