/* 
 * Riemannian Wavefield Extrapolation: 
 * shot-record migration 
 * pcs 2005
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
    sf_file Fw_s=NULL, Fw_r=NULL;
    sf_file Fi=NULL, Fm=NULL, Fr=NULL;
    axa ag,at,aw,ar,aj;
    int ig,it,iw,ir;

    int method;
    bool verb;

    complex float **dat_s, *wfl_s;
    complex float **dat_r, *wfl_r;
    float         **img;
    float         **aa,**bb,**mm;

    complex float **ab;
    float         **a0,**b0;

    float w,ws,wr;
    char *met="";
    
    sf_init(argc,argv);

    if(! sf_getbool("verb", &verb))     verb=false;
    if(! sf_getint("method",&method)) method=0;    /* extrapolation method */
						
    Fm = sf_input("abm");
    Fr = sf_input("abr");

    iaxa(Fm,&at,2); at.l="t"; /* 'extrapolation axis' (can be time) */
    iaxa(Fr,&ar,1); ar.l="r"; /* a,b reference */
    if(method==0) ar.n=1; /* pure F-D */

    aj.n=1;
    aj.o=0.;
    aj.d=1.;
    aj.l="";

    Fw_s = sf_input ( "in");
    Fw_r = sf_input ("rwf");
    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    if (SF_COMPLEX !=sf_gettype(Fw_s)) sf_error("Need complex source");
    if (SF_COMPLEX !=sf_gettype(Fw_r)) sf_error("Need complex data");
    
    iaxa(Fw_s,&ag,1); ag.l="g"; /* 'position axis' (can be angle) */
    iaxa(Fw_s,&at,2); at.l="t";
    iaxa(Fw_s,&aw,3); aw.l="w"; /* frequency */
    
    oaxa(Fi,&ag,1);
    oaxa(Fi,&at,2);
    oaxa(Fi,&aj,3);
    
    dat_s = sf_complexalloc2(ag.n,at.n);
    dat_r = sf_complexalloc2(ag.n,at.n);
    wfl_s = sf_complexalloc (ag.n);
    wfl_r = sf_complexalloc (ag.n);
    img   = sf_floatalloc2  (ag.n,at.n);
    
    if(verb) {
	raxa(ag);
	raxa(at);
	raxa(aw);
	raxa(ar);
    }

    /* read ABM */
    aa = sf_floatalloc2  (ag.n,at.n);
    bb = sf_floatalloc2  (ag.n,at.n);
    mm = sf_floatalloc2  (ag.n,at.n);

    sf_floatread(aa[0],ag.n*at.n,Fm); /* a coef */
    sf_floatread(bb[0],ag.n*at.n,Fm); /* b coef */
    sf_floatread(mm[0],ag.n*at.n,Fm); /* mask */

    /* read ABr */
    ab = sf_complexalloc2(ar.n,at.n);
    a0 = sf_floatalloc2  (ar.n,at.n);
    b0 = sf_floatalloc2  (ar.n,at.n);

    sf_complexread(ab[0],ar.n*at.n,Fr);
    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
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
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    rweone_init(ag,at,aw,ar,method,verb);
    switch(method) {
	case 3: rweone_psc_coef(aa,bb,a0,b0); break;
	case 2: rweone_ffd_coef(aa,bb,a0,b0); break;
	case 1: ;/* SSF */                    break;
	case 0: rweone_xfd_coef(aa,bb);       break;
    }

    for(it=0;it<at.n;it++) {
	for(ig=0;ig<ag.n;ig++) {
	    img[it][ig] = 0.;
	}
    }

/*------------------------------------------------------------*/        
    for(iw=0;iw<aw.n;iw++) {
	w=aw.o+iw*aw.d;
	ws = -w; /*      causal */
	wr = +w; /* anti-causal */

	sf_warning("%s %d %d",met,iw,aw.n);
	
	sf_complexread(dat_s[0],ag.n*at.n,Fw_s);
	sf_complexread(dat_r[0],ag.n*at.n,Fw_r);
    
	for(ig=0;ig<ag.n;ig++) {
	    wfl_s[ig] = 0;
	    wfl_r[ig] = 0;
	}
	
	for(it=0;it<=at.n-2;it++) {
	    for(ig=0;ig<ag.n;ig++) {
		wfl_s[ig] += dat_s[it][ig];
		wfl_r[ig] += dat_r[it][ig];
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
	it=at.n-1; rweone_spi(wfl_s,wfl_r,img[it]);
    }
/*------------------------------------------------------------*/
    
    sf_floatwrite  (img[0],ag.n*at.n,Fi);
}
