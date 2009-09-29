/* Riemannian Wavefield Extrapolation: shot-record migration. */

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
    sf_file Fw_s=NULL, Fw_r=NULL;
    sf_file Fi=NULL, Fm=NULL, Fr=NULL;
    sf_axis ag,at,aw,ar;
    int ig,it,iw,ir;
    int ng,nt,nw,nr;

    int method;
    bool  verb;
    bool   adj;

    sf_complex **dat_s=NULL, *wfl_s=NULL;
    sf_complex **dat_r=NULL, *wfl_r=NULL;
    float      **img=NULL;
    float      **aa=NULL,**bb=NULL,**mm=NULL;

    sf_complex **ab=NULL;
    float      **a0=NULL,**b0=NULL;

    float w,ws,wr,w0,dw;
    char *met="";

    sf_init(argc,argv);

    if(! sf_getbool("verb", &verb))     verb=false;
    if(! sf_getint("method",&method)) method=0;    /* extrapolation method */
    if(! sf_getbool("adj",  &adj))       adj=false;/* y=modeling; n=migration */
    Fm = sf_input("abm");
    Fr = sf_input("abr");

    at=sf_iaxa(Fm,2); sf_setlabel(at,"t"); /* 'extrapolation axis' */
    ar=sf_iaxa(Fr,1); sf_setlabel(ar,"r"); /* a,b reference */
    if(method==0) sf_setn(ar,1); /* pure F-D */
    nr=sf_n(ar); 

    Fw_s = sf_input ( "in");
    if (SF_COMPLEX !=sf_gettype(Fw_s)) sf_error("Need complex source");

    /* 'position axis' (could be angle) */
    ag = sf_iaxa(Fw_s,1); ng=sf_n(ag); sf_setlabel(ag,"g"); 

    /* 'extrapolation axis' (could be time) */
    at = sf_iaxa(Fw_s,2); nt=sf_n(at); sf_setlabel(at,"t");
    aw = sf_iaxa(Fw_s,3); sf_setlabel(aw,"w"); /* frequency */

    if(adj) { /* modeling */
	Fw_r = sf_output("out"); 
	sf_settype(Fw_r,SF_COMPLEX);

	Fi   = sf_input ("img");
	if (SF_FLOAT !=sf_gettype(Fi)) sf_error("Need float image");

	sf_oaxa(Fw_r,ag,1);
	sf_oaxa(Fw_r,at,2);
	sf_oaxa(Fw_r,aw,3);
    } else {  /* migration */
	Fw_r = sf_input ("rwf");
	if (SF_COMPLEX !=sf_gettype(Fw_r)) sf_error("Need complex data");

	Fi   = sf_output("out"); 
	sf_settype(Fi,SF_FLOAT);
		
	sf_oaxa(Fi,ag,1);
	sf_oaxa(Fi,at,2);
	sf_putint(Fi,"n3",1);
    }

    img   = sf_floatalloc2  (ng,nt);
    dat_s = sf_complexalloc2(ng,nt);
    dat_r = sf_complexalloc2(ng,nt);
    wfl_s = sf_complexalloc (ng);
    wfl_r = sf_complexalloc (ng);

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
    nw = sf_n(aw);
    dw = sf_d(aw) * 2.*SF_PI; 
    w0 = sf_o(aw) * 2.*SF_PI;

    rweone_init(ag,at,aw,ar,method,verb);
    switch(method) {
	case 3: rweone_psc_coef(aa,bb,a0,b0); break;
	case 2: rweone_ffd_coef(aa,bb,a0,b0); break;
	case 1: ;/* SSF */                    break;
	case 0: rweone_xfd_coef(aa,bb);       break;
    }

    if(adj) { /* modeling */
	
    } else {  /* migration */
	for(it=0;it<nt;it++) {
	    for(ig=0;ig<ng;ig++) {
		img[it][ig] = 0.;
	    }
	}
    }

/*------------------------------------------------------------*/
    if( adj)  sf_floatread  (img[0],ng*nt,Fi);

    for(iw=0;iw<nw;iw++) {
	w=w0+iw*dw;
	sf_warning("%s %d %d",met,iw,nw);

	if(adj) {
	    sf_complexread(dat_s[0],ng*nt,Fw_s);

	    ws = -w; /*      causal */
	    for(ig=0;ig<ng;ig++) {
		wfl_s[ig] = sf_cmplx(0.,0.);
	    }
	    for(it=0;it<nt;it++) {
		for(ig=0;ig<ng;ig++) {
#ifdef SF_HAS_COMPLEX_H
		    wfl_s[ig] += dat_s[it][ig];
		    dat_r[it][ig] = wfl_s[ig]*img[it][ig];
#else
		    wfl_s[ig] = sf_cadd(wfl_s[ig],dat_s[it][ig]);
		    dat_r[it][ig] = sf_crmul(wfl_s[ig],img[it][ig]);
#endif
		}
		
		if(method!=0) rweone_fk(ws,wfl_s,aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(ws,wfl_s,aa[it],it);
		rweone_tap(wfl_s);
	    }

	    for(it=0;it<nt;it++) {
		for(ig=0;ig<ng;ig++) {
		    dat_s[it][ig] = dat_r[it][ig];
		}
	    }

	    wr = -w; /*      causal */
	    for(ig=0;ig<ng;ig++) {
		wfl_r[ig] = sf_cmplx(0.,0.);
	    }
	    for(it=nt-1;it>=0;it--) {
		for(ig=0;ig<ng;ig++) {
#ifdef SF_HAS_COMPLEX_H
		    wfl_r[ig] += dat_s[it][ig];
#else
		    wfl_r[ig] = sf_cadd(wfl_r[ig],dat_s[it][ig]);
#endif
		    dat_r[it][ig] = wfl_r[ig];
		}
		
		if(method!=0) rweone_fk(wr,wfl_r,aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(wr,wfl_r,aa[it],it);
		rweone_tap(wfl_r);
	    }

	    sf_complexwrite(dat_r[0],ng*nt,Fw_r);

	} else {
	    ws = -w; /*      causal */
	    wr = +w; /* anti-causal */

	    sf_complexread(dat_s[0],ng*nt,Fw_s);
	    sf_complexread(dat_r[0],ng*nt,Fw_r);

	    for(ig=0;ig<ng;ig++) {
		wfl_s[ig] = sf_cmplx(0,0);
		wfl_r[ig] = sf_cmplx(0,0);
	    }

	    for(it=0;it<=nt-2;it++) {
		for(ig=0;ig<ng;ig++) {
#ifdef SF_HAS_COMPLEX_H
		    wfl_s[ig] += dat_s[it][ig];
		    wfl_r[ig] += dat_r[it][ig];
#else
		    wfl_s[ig] = sf_cadd(wfl_s[ig],dat_s[it][ig]);
		    wfl_r[ig] = sf_cadd(wfl_r[ig],dat_r[it][ig]);
#endif
		}

		rweone_spi(wfl_s,wfl_r,img[it]);
		
		if(method!=0) {
		    rweone_fk(ws,wfl_s,aa[it],a0[it],b0[it],mm[it],it);
		    rweone_fk(wr,wfl_r,aa[it],a0[it],b0[it],mm[it],it);
		} else {
		    rweone_fx(ws,wfl_s,aa[it],it);
		    rweone_fx(wr,wfl_r,aa[it],it);
		}
	    }
	    it=nt-1; rweone_spi(wfl_s,wfl_r,img[it]);
	}
    }

    if(!adj)  sf_floatwrite  (img[0],ng*nt,Fi);
/*------------------------------------------------------------*/

    exit(0);
}
