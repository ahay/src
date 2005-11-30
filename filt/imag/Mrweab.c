/* 
 * Riemannian Wavefield Extrapolation:
 * a,b coefficients
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

#define LOOPRC(a) for(it=0;it<at.n;it++){ for(ig=0;ig<ag.n;ig++){ {a} }}

int main(int argc, char* argv[])
{

    sf_file Fi=NULL, Fo=NULL, Fr=NULL, Fs=NULL; /* I/O files */

    complex float **rays;
    complex float **ab;

    float **x, **z;
    float **h1, **h2;
    float  *qq, *gg;
    float **ss, **aa, **bb;
    float **ma, **mb, **mm;

    axa at,ag,ar;
    int it,ig;
    float gx,gz;
    float eps,peps;

    bool verb;
    int naref, nbref;
    float mina, maxa, dela;
    float minb, maxb, delb;
    float a0,b0;
    int   ia,ib,ii;

    float tiny;

    sf_init(argc,argv);

    Fi = sf_input ( "in");
    Fs = sf_input ("slo");
    Fo = sf_output("out");
    Fr = sf_output("abr");

    if(! sf_getbool("verb",&verb)) verb=false;

    if(! sf_getint("naref",&naref)) naref=1;
    if(! sf_getint("nbref",&nbref)) nbref=1;

    if(verb) sf_warning("naref=%d nbref=%d",naref,nbref);

    iaxa(Fi,&ag,1); if(verb) raxa(ag);
    iaxa(Fi,&at,2); if(verb) raxa(at);

    ag.l="g"; oaxa(Fo,&ag,1);
    at.l="t"; oaxa(Fo,&at,2);

    sf_settype(Fo,SF_FLOAT);
    sf_putint(Fo,"n3",5);

    if(! sf_getfloat("peps",&peps)) peps=0.01;

    ar.n=naref*nbref;
    ar.o=0.;
    ar.d=1.;
    ar.l="r";
    oaxa(Fr,&ar,1);
    oaxa(Fr,&at,2);
    sf_settype(Fr,SF_COMPLEX);
    sf_putfloat(Fr,"esize",8);

    rays = sf_complexalloc2(ag.n,at.n);
    ab   = sf_complexalloc2(naref*nbref,at.n);

    x    = sf_floatalloc2(ag.n,at.n);
    z    = sf_floatalloc2(ag.n,at.n);
    h1   = sf_floatalloc2(ag.n,at.n);
    h2   = sf_floatalloc2(ag.n,at.n);

    ss   = sf_floatalloc2(ag.n,at.n);
    aa   = sf_floatalloc2(ag.n,at.n);
    bb   = sf_floatalloc2(ag.n,at.n);
    ma   = sf_floatalloc2(ag.n,at.n);
    mb   = sf_floatalloc2(ag.n,at.n);
    mm   = sf_floatalloc2(ag.n,at.n);

    qq   = sf_floatalloc (ag.n*at.n);
    gg   = sf_floatalloc (ag.n);

    sf_floatread  (ss[0],  ag.n*at.n,Fs);  /* slowness */
    sf_complexread(rays[0],ag.n*at.n,Fi);  /* rays */
    LOOPRC( z[it][ig] = cimagf(rays[it][ig]);
	    x[it][ig] = crealf(rays[it][ig]); );
    
    /* h1=alpha */
    for(it=0;it<at.n-1;it++) {
	for(ig=0;ig<ag.n;ig++) {
	    gx = (x[it+1][ig] - x[it][ig]) / at.d;
	    gz = (z[it+1][ig] - z[it][ig]) / at.d;
	    h1[it][ig] = sqrtf(gx*gx+gz*gz);
	}
    }
    for(ig=0;ig<ag.n;ig++) {
	h1[at.n-1][ig] = h1[at.n-2][ig];
    }

    /* h2=J */
    for(it=0;it<at.n;it++) {
	for(ig=0;ig<ag.n-1;ig++) {
	    gx = (x[it][ig+1] - x[it][ig]) / (ag.d);
	    gz = (z[it][ig+1] - z[it][ig]) / (ag.d);
	    h2[it][ig] = sqrtf(gx*gx+gz*gz);
	}
    }
    for(it=0;it<at.n;it++) {
	h2[it][ag.n-1] = h2[it][ag.n-2];
    }

    /* avoid small h2=J */
    ii=0; LOOPRC( qq[ii] = SF_ABS(h2[it][ig]); ii++; );
/*    eps = peps * sf_quantile(at.n*ag.n/2,at.n*ag.n,qq);*/
    eps = peps * (
	sf_quantile(          0,at.n*ag.n,qq) + 
	sf_quantile(at.n*ag.n-1,at.n*ag.n,qq)
	) /2.;
    LOOPRC( if(SF_ABS(h2[it][ig]) < eps) h2[it][ig]=eps; );

    LOOPRC( 
	aa[it][ig] = ss[it][ig] * h1[it][ig];
	bb[it][ig] = h1[it][ig] / h2[it][ig];
	mm[it][ig] = 1.;
	ma[it][ig] = 1.;
	mb[it][ig] = 1.;
	);

    tiny=0.1;

    /* compute reference a,b */
    for(it=0;it<at.n;it++) {

	for(ig=0;ig<ag.n;ig++){
	    gg[ig] = aa[it][ig];
	}
	mina = sf_quantile(     0,ag.n,gg);
	maxa = sf_quantile(ag.n-1,ag.n,gg);
	dela = (maxa-mina)/naref;

	for(ig=0;ig<ag.n;ig++){
	    gg[ig] = bb[it][ig];
	}
	minb = sf_quantile(     0,ag.n,gg);
	maxb = sf_quantile(ag.n-1,ag.n,gg);
	delb = (maxb-minb)/nbref;

	/* reference a and b */
	ii=0;
	for(ia=0;ia<naref;ia++) {
	    a0 = mina + 0.5*dela + ia*dela;
	    for(ib=0;ib<nbref;ib++) {
		b0 = minb + 0.5*delb + ib*delb;
		ab[it][ii] = a0+I*b0;
		ii++;
	    }
	}
	
	/* mask for a */
	for(ia=0;ia<naref;ia++) {
	    a0 = mina + 0.5*dela + ia*dela;
	    for(ig=0;ig<ag.n;ig++) {
		if( a0-(0.5+tiny)*dela <= aa[it][ig] &&
		    a0+(0.5+tiny)*dela >= aa[it][ig]  ) ma[it][ig] = ia;
	    }
	}
	/* mask for b */
	for(ib=0;ib<nbref;ib++) {
	    b0 = minb + 0.5*delb + ib*delb;
	    for(ig=0;ig<ag.n;ig++) {
		if( b0-(0.5+tiny)*delb <= bb[it][ig] &&
		    b0+(0.5+tiny)*delb >= bb[it][ig]  ) mb[it][ig] = ib;
	    }
	}
    }

    /* reference medium mask */
    ii=0;
    for(ia=0;ia<naref;ia++) {
	for(ib=0;ib<nbref;ib++) {
	    
	    LOOPRC(
		if( ma[it][ig]==ia &&
		    mb[it][ig]==ib  ) mm[it][ig]=ii;
		);
	    ii++;
	}
    }

    /* for compatibility with older RWE programs */
    LOOPRC(
	ma[it][ig] +=1;
	mb[it][ig] +=1;
	mm[it][ig] +=1;
	);
    
    sf_floatwrite(aa[0],ag.n*at.n,Fo);
    sf_floatwrite(bb[0],ag.n*at.n,Fo);
    sf_floatwrite(mm[0],ag.n*at.n,Fo);    
    sf_floatwrite(ma[0],ag.n*at.n,Fo);
    sf_floatwrite(mb[0],ag.n*at.n,Fo);
    
    sf_complexwrite(ab[0],naref*nbref*at.n,Fr);
}
