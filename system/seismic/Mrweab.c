/* Riemannian Wavefield Extrapolation: a,b coefficients */
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

#define LOOPRC(a) for(it=0;it<nt;it++){ for(ig=0;ig<ng;ig++){ {a} }}

int main(int argc, char* argv[])
{

    sf_file Fi=NULL, Fo=NULL, Fr=NULL, Fs=NULL; /* I/O files */

    sf_complex **rays=NULL;
    sf_complex **ab=NULL;

    float **x=NULL, **z=NULL;
    float **h1=NULL, **h2=NULL;
    float  *qq=NULL, *gg=NULL;
    float **ss=NULL, **aa=NULL, **bb=NULL;
    float **ma=NULL, **mb=NULL, **mm=NULL;

    sf_axis at,ag,ar;
    int it,ig;
    int nt,ng;
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

    ag = sf_iaxa(Fi,1); sf_setlabel(ag,"g"); ng=sf_n(ag); if(verb) sf_raxa(ag);
    at = sf_iaxa(Fi,2); sf_setlabel(at,"t"); nt=sf_n(at); if(verb) sf_raxa(at);

    sf_oaxa(Fo,ag,1);
    sf_oaxa(Fo,at,2);

    sf_settype(Fo,SF_FLOAT);
    sf_putint(Fo,"n3",5);

    if(! sf_getfloat("peps",&peps)) peps=0.01;

    ar=sf_maxa(naref*nbref,0.,1.);
    sf_setlabel(ar,"r");
    sf_oaxa(Fr,ar,1);
    sf_oaxa(Fr,at,2);
    sf_settype(Fr,SF_COMPLEX);
    sf_putfloat(Fr,"esize",8);

    rays = sf_complexalloc2(ng,nt);
    ab   = sf_complexalloc2(naref*nbref,nt);

    x    = sf_floatalloc2(ng,nt);
    z    = sf_floatalloc2(ng,nt);
    h1   = sf_floatalloc2(ng,nt);
    h2   = sf_floatalloc2(ng,nt);

    ss   = sf_floatalloc2(ng,nt);
    aa   = sf_floatalloc2(ng,nt);
    bb   = sf_floatalloc2(ng,nt);
    ma   = sf_floatalloc2(ng,nt);
    mb   = sf_floatalloc2(ng,nt);
    mm   = sf_floatalloc2(ng,nt);

    qq   = sf_floatalloc (ng*nt);
    gg   = sf_floatalloc (ng);

    sf_floatread  (ss[0],  ng*nt,Fs);  /* slowness */
    sf_complexread(rays[0],ng*nt,Fi);  /* rays */
    LOOPRC( z[it][ig] = cimagf(rays[it][ig]);
	    x[it][ig] = crealf(rays[it][ig]); );

    /* h1=alpha */
    for(it=0;it<nt-1;it++) {
	for(ig=0;ig<ng;ig++) {
	    gx = (x[it+1][ig] - x[it][ig]) / sf_d(at);
	    gz = (z[it+1][ig] - z[it][ig]) / sf_d(at);
	    h1[it][ig] = hypotf(gx,gz);
	}
    }
    for(ig=0;ig<ng;ig++) {
	h1[nt-1][ig] = h1[nt-2][ig];
    }

    /* h2=J */
    for(it=0;it<nt;it++) {
	for(ig=0;ig<ng-1;ig++) {
	    gx = (x[it][ig+1] - x[it][ig]) / sf_d(ag);
	    gz = (z[it][ig+1] - z[it][ig]) / sf_d(ag);
	    h2[it][ig] = hypotf(gx,gz);
	}
    }
    for(it=0;it<nt;it++) {
	h2[it][ng-1] = h2[it][ng-2];
    }

    /* avoid small h2=J */
    ii=0; LOOPRC( qq[ii] = SF_ABS(h2[it][ig]); ii++; );
/*    eps = peps * sf_quantile(nt*ng/2,nt*ng,qq);*/
    eps = peps * (
	sf_quantile(      0,nt*ng,qq) +
	sf_quantile(nt*ng-1,nt*ng,qq)
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
    for(it=0;it<nt;it++) {

	for(ig=0;ig<ng;ig++){
	    gg[ig] = aa[it][ig];
	}
	mina = sf_quantile(     0,ng,gg);
	maxa = sf_quantile(ng-1,ng,gg);
	dela = (maxa-mina)/naref;

	for(ig=0;ig<ng;ig++){
	    gg[ig] = bb[it][ig];
	}
	minb = sf_quantile(     0,ng,gg);
	maxb = sf_quantile(ng-1,ng,gg);
	delb = (maxb-minb)/nbref;

	/* reference a and b */
	ii=0;
	for(ia=0;ia<naref;ia++) {
	    a0 = mina + 0.5*dela + ia*dela;
	    for(ib=0;ib<nbref;ib++) {
		b0 = minb + 0.5*delb + ib*delb;
		ab[it][ii] = sf_cmplx(a0,b0);
		ii++;
	    }
	}

	/* mask for a */
	for(ia=0;ia<naref;ia++) {
	    a0 = mina + 0.5*dela + ia*dela;
	    for(ig=0;ig<ng;ig++) {
		if( a0-(0.5+tiny)*dela <= aa[it][ig] &&
		    a0+(0.5+tiny)*dela >= aa[it][ig]  ) ma[it][ig] = ia;
	    }
	}
	/* mask for b */
	for(ib=0;ib<nbref;ib++) {
	    b0 = minb + 0.5*delb + ib*delb;
	    for(ig=0;ig<ng;ig++) {
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

    sf_floatwrite(aa[0],ng*nt,Fo);
    sf_floatwrite(bb[0],ng*nt,Fo);
    sf_floatwrite(mm[0],ng*nt,Fo);
    sf_floatwrite(ma[0],ng*nt,Fo);
    sf_floatwrite(mb[0],ng*nt,Fo);

    sf_complexwrite(ab[0],naref*nbref*nt,Fr);


    exit(0);
}
