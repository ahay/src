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
#include "rwespm.h"

int main(int argc, char* argv[])
{
    sf_file Fw_s=NULL, Fw_r=NULL;
    sf_file Fi=NULL, Fm=NULL, Fr=NULL;
    axa ag,at,aw,ar;
    int ig,it,   ir;

    int method;
    bool verb;

    complex float **wfl_s;
    complex float **wfl_r;
    float         **img;
    float         **aa,**bb,**mm;

    complex float **ab;
    float         **a0,**b0;

    sf_init(argc,argv);

    if(! sf_getbool("verb", &verb))     verb=false;
    if(! sf_getint("method",&method)) method=0;    /* extrapolation method */
						
    Fm = sf_input("abm");
    Fr = sf_input("abr");

    iaxa(Fm,&at,2); at.l="t"; /* 'extrapolation axis' (can be time) */
    iaxa(Fr,&ar,1); ar.l="r"; /* a,b reference */
    if(method==0) ar.n=1; /* pure F-D */

    Fw_s = sf_input ( "in");
    Fw_r = sf_input ("rwf");
    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    if (SF_COMPLEX !=sf_gettype(Fw_s)) sf_error("Need complex source");
    if (SF_COMPLEX !=sf_gettype(Fw_r)) sf_error("Need complex data");
    
    iaxa(Fw_s,&ag,1); ag.l="g"; /* 'position axis' (can be angle) */
    iaxa(Fw_s,&aw,2); aw.l="w"; /* frequency */
    
    oaxa(Fi,&ag,1);
    oaxa(Fi,&at,2);
    
    wfl_s = sf_complexalloc2(ag.n,aw.n);
    wfl_r = sf_complexalloc2(ag.n,aw.n);
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

    sf_complexread(wfl_s[0],ag.n*aw.n,Fw_s);
    sf_complexread(wfl_r[0],ag.n*aw.n,Fw_r);
    
    for(it=0;it<at.n;it++) {
	for(ig=0;ig<ag.n;ig++) {
	    img[it][ig] = 0.;
	}
    }
    
    /*------------------------------------------------------------*/
    /* execute */
    rwespm_init(ag,at,aw,ar,method);
    
    rwespm_main(wfl_s,wfl_r,img,aa,bb,mm,a0,b0);
    /* execute */
    /*------------------------------------------------------------*/
    
    sf_floatwrite  (img[0],ag.n*at.n,Fi);
}
