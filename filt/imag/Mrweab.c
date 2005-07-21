/* Riemannian Wavefield Extrapolation a,b coefficients */
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

#define LOOPRC(a) for(ig=0;ig<ag.n;ig++){ for(it=0;it<at.n;it++){ {a} }}

int main(int argc, char* argv[])
{

    sf_file Fi, Fo, Fr, Fs; /* I/O files */

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

    Fi = sf_input("in");
    Fs = sf_input("slo");
    Fo = sf_output("out");
    Fr = sf_output("abr");

    if(! sf_getbool("verb",&verb)) verb=false;

    if(! sf_getint("naref",&naref)) naref=1;
    if(! sf_getint("nbref",&nbref)) nbref=1;

    if(verb) sf_warning("naref=%d nbref=%d",naref,nbref);

    iaxa(Fi,&at,1); if(verb) raxa(at);
    iaxa(Fi,&ag,2); if(verb) raxa(ag);

    sf_settype(Fo,SF_FLOAT);
    sf_putint(Fo,"n3",5);

    if(! sf_getfloat("peps",&peps)) peps=0.01;

    ar.n=naref*nbref;
    ar.o=0.;
    ar.d=1.;
    oaxa(Fr,&at,1);
    oaxa(Fr,&ar,2);
    sf_settype(Fr,SF_COMPLEX);
    sf_putfloat(Fr,"esize",8);

    rays = sf_complexalloc2(at.n,ag.n);
    ab   = sf_complexalloc2(at.n,naref*nbref);

    x    = sf_floatalloc2(at.n,ag.n);
    z    = sf_floatalloc2(at.n,ag.n);
    h1   = sf_floatalloc2(at.n,ag.n);
    h2   = sf_floatalloc2(at.n,ag.n);

    ss   = sf_floatalloc2(at.n,ag.n);
    aa   = sf_floatalloc2(at.n,ag.n);
    bb   = sf_floatalloc2(at.n,ag.n);
    ma   = sf_floatalloc2(at.n,ag.n);
    mb   = sf_floatalloc2(at.n,ag.n);
    mm   = sf_floatalloc2(at.n,ag.n);

    qq   = sf_floatalloc (at.n*ag.n);
    gg   = sf_floatalloc (     ag.n);

    sf_floatread  (ss[0],  at.n*ag.n,Fs);  /* slowness */
    sf_complexread(rays[0],at.n*ag.n,Fi);  /* rays */
    LOOPRC( z[ig][it] = cimagf(rays[ig][it]);
	    x[ig][it] = crealf(rays[ig][it]); );
    
    /* h1=alpha */
    for(ig=0;ig<ag.n;ig++) {
	for(it=0;it<at.n-1;it++) {
	    gx = (x[ig][it+1] - x[ig][it]) / at.d;
	    gz = (z[ig][it+1] - z[ig][it]) / at.d;
	    h1[ig][it] = sqrtf(gx*gx+gz*gz);
	}
    }
    for(ig=0;ig<ag.n;ig++) {
	h1[ig][at.n-1] = h1[ig][at.n-2];
    }

    /* h2=J */
    for(ig=1;ig<ag.n-1;ig++) {
	for(it=0;it<at.n;it++) {
	    gx = (x[ig+1][it] - x[ig-1][it]) / (2*ag.d);
	    gz = (z[ig+1][it] - z[ig-1][it]) / (2*ag.d);
	    h2[ig][it] = sqrtf(gx*gx+gz*gz);
	}
    }
    for(it=0;it<at.n;it++) {
	h2[     0][it] = h2[     1][it];
	h2[ag.n-1][it] = h2[ag.n-2][it];
    }

    /* avoid small h2=J */
    ii=0; LOOPRC( qq[ii] = SF_ABS(h2[ig][it]); ii++; );
/*    eps = peps * sf_quantile(at.n*ag.n/2,at.n*ag.n,qq);*/
    eps = peps * (
	sf_quantile(          0,at.n*ag.n,qq) + 
	sf_quantile(at.n*ag.n-1,at.n*ag.n,qq)
	) /2.;
    LOOPRC( if(SF_ABS(h2[ig][it]) < eps) h2[ig][it]=eps; );

    LOOPRC( 
	aa[ig][it] = ss[ig][it] * h1[ig][it];
	bb[ig][it] = h1[ig][it] / h2[ig][it];
	mm[ig][it] = 1.;
	ma[ig][it] = 1.;
	mb[ig][it] = 1.;
	);

    tiny=0.1;

    /* compute reference a,b */
    for(it=0;it<at.n;it++) {

	for(ig=0;ig<ag.n;ig++){
	    gg[ig] = aa[ig][it];
	}
	mina = sf_quantile(     0,ag.n,gg);
	maxa = sf_quantile(ag.n-1,ag.n,gg);
	dela = (maxa-mina)/naref;

	for(ig=0;ig<ag.n;ig++){
	    gg[ig] = bb[ig][it];
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
		ab[ii][it] = a0+I*b0;
		ii++;
	    }
	}
	
	/* mask for a */
	for(ia=0;ia<naref;ia++) {
	    a0 = mina + 0.5*dela + ia*dela;
	    for(ig=0;ig<ag.n;ig++) {
		if( a0-(0.5+tiny)*dela <= aa[ig][it] &&
		    a0+(0.5+tiny)*dela >= aa[ig][it]  ) ma[ig][it] = ia;
	    }
	}
	/* mask for b */
	for(ib=0;ib<nbref;ib++) {
	    b0 = minb + 0.5*delb + ib*delb;
	    for(ig=0;ig<ag.n;ig++) {
		if( b0-(0.5+tiny)*delb <= bb[ig][it] &&
		    b0+(0.5+tiny)*delb >= bb[ig][it]  ) mb[ig][it] = ib;
	    }
	}
    }

    /* reference medium mask */
    ii=0;
    for(ia=0;ia<naref;ia++) {
	for(ib=0;ib<nbref;ib++) {
	    
	    LOOPRC(
		if( ma[ig][it]==ia &&
		    mb[ig][it]==ib  ) mm[ig][it]=ii;
		);
	    ii++;
	}
    }

    sf_floatwrite(aa[0],at.n*ag.n,Fo);
    sf_floatwrite(bb[0],at.n*ag.n,Fo);
    sf_floatwrite(mm[0],at.n*ag.n,Fo);    
    sf_floatwrite(ma[0],at.n*ag.n,Fo);
    sf_floatwrite(mb[0],at.n*ag.n,Fo);

    sf_complexwrite(ab[0],at.n*naref*nbref,Fr);
}
