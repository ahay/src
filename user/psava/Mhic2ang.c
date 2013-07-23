/* angle decomposition of CIPs */

/*
  Copyright (C) 2010 Colorado School of Mines

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "aniang.h"

int main(int argc, char* argv[])
{
    bool verb;
    bool  adj;
    bool anis;
    sf_file Fcip=NULL;	/*   lag-domain CIPs */
    sf_file Fang=NULL;	/* angle-domain CIPs */
    sf_file Fvel=NULL;  /*   velocity @ CIPs */
    sf_file Fnor=NULL;	/*    normals @ CIPs */
    sf_file Ftlt=NULL;	/*       tilt @ CIPs */
    sf_file Fani=NULL;  /* anisotropy @ CIPs */

    sf_axis ahx,ahy,ahz,aht,ac,ath,aph,aps,aj;
    int     ihx,ihy,ihz,iht,ic,ith,iph;

    
    /* angle parameters */
    int   nth,nph,nps,nhx,nhy,nhz,nht;
    float oth,oph,ops,ohx,ohy,ohz,oht;
    float dth,dph,dps,dhx,dhy,dhz,dht;
    float phi;
    float tht;
    float psi;
    float v_s,v_r;
    float cosum,codif,sitovel;

    /*  arrays                  1   2   3   4  */
    float     ****cip;      /* nhx-nhy-nhz-nht */
    float       **ang;      /* nph-nth         */
    float        *vep;      /* nc              */
    float        *ves;      /* nc              */
    float        *eps=NULL; /* nc              */
    float        *dlt=NULL; /* nc              */

    vc3d           vv;      /* azimuth reference vector */
    vc3d          *nn;      /* normal vectors  */
    vc3d          *tt=NULL; /*   tilt vectors  */
    vc3d          *aa;      /* in-plane reference vector */
    vc3d           qq;
    vc3d           jk;      /* temp vector */
    float    hx,hy,hz;

    float tau; /* time lag */
    int   jht; /* tau axis index */
    float fht; /* tau axis weight */

    float ssn; /* slant-stack normalization */

    float *ttipar;
    /*-----------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);    
    
#ifdef _OPENMP
    omp_init(); /* OMP parameters */
#endif

    if(! sf_getbool("verb",&verb)) verb=false;	/* verbosity flag */
    if(! sf_getbool("anis",&anis)) anis=false;	/* anisotropy flag */
    if(! sf_getbool("adj", &adj))   adj=true;	/* adj flag */    
    /* 
     * ADJ: cip to ang
     * FOR: ang to cip
     */

    sf_warning("verb=%d",verb);
    sf_warning("anis=%d",anis);

    /* select anisotropy model */
    if(anis) sf_warning("ANI model");
    else     sf_warning("ISO model");

    if(adj) {
	Fcip=sf_input ( "in"); /* CIP file */
	Fang=sf_output("out"); /* ANG file */
    } else {
	Fcip=sf_output("out"); /* CIP file */
	Fang=sf_input ("in");  /* ANG file */
    }
    Fvel=sf_input ("vel");     /* velocity file  */
    Fnor=sf_input ("nor");     /* normal vectors */
    if(anis) {
	Ftlt=sf_input ("tlt"); /*   tilt vectors */
	Fani=sf_input ("ani"); /*     anisotropy */
    }

    aj  = sf_maxa(1,0,1);

    if(adj) {
	/* input axes */
	ahx = sf_iaxa(Fcip,1); sf_setlabel(ahx,"hx");
	ahy = sf_iaxa(Fcip,2); sf_setlabel(ahy,"hy");
	ahz = sf_iaxa(Fcip,3); sf_setlabel(ahz,"hz");
	aht = sf_iaxa(Fcip,4); sf_setlabel(aht,"ht");

	/* CIP axis */
	ac  = sf_iaxa(Fcip,5); sf_setlabel(ac ,"c ");
		
	/* reflection angle */
	if(! sf_getint  ("nth",&nth)) nth=90;
	if(! sf_getfloat("oth",&oth)) oth=0;
	if(! sf_getfloat("dth",&dth)) dth=1.;
	ath = sf_maxa(nth,oth,dth);
	sf_setlabel(ath,"th");
	sf_setunit (ath,"deg");
	
	/* azimuth angle */
	if(! sf_getint  ("nph",&nph)) nph=360;
	if(! sf_getfloat("oph",&oph)) oph=-180;
	if(! sf_getfloat("dph",&dph)) dph=1.;
	aph = sf_maxa(nph,oph,dph);
	sf_setlabel(aph,"ph");
	sf_setunit (aph,"deg");

	/* output axes */
	sf_oaxa(Fang,ath,1);
	sf_oaxa(Fang,aph,2);
	sf_oaxa(Fang,ac ,3);
	sf_oaxa(Fang,aj ,4);
	sf_oaxa(Fang,aj ,5);

    } else {

	/* lag in x */
	if(! sf_getint  ("nhx",&nhx)) nhx=1;
	if(! sf_getfloat("ohx",&ohx)) ohx=0;
	if(! sf_getfloat("dhx",&dhx)) dhx=1.;
	ahx = sf_maxa(nhx,ohx,dhx);
	sf_setlabel(ahx,"hx");
	sf_setunit (ahx,"");

	/* lag in y */
	if(! sf_getint  ("nhy",&nhy)) nhy=1;
	if(! sf_getfloat("ohy",&ohy)) ohy=0;
	if(! sf_getfloat("dhy",&dhy)) dhy=1.;
	ahy = sf_maxa(nhy,ohy,dhy);
	sf_setlabel(ahy,"hy");
	sf_setunit (ahy,"");

	/* lag in z */
	nhz=1;
	ohz=0.;
	dhz=1.;
	ahz = sf_maxa(nhz,ohz,dhz);
	sf_setlabel(ahz,"hz");
	sf_setunit (ahz,"");

	/* lag in t */
	if(! sf_getint  ("nht",&nht)) nht=1;
	if(! sf_getfloat("oht",&oht)) oht=0.;
	if(! sf_getfloat("dht",&dht)) dht=1.;
	aht = sf_maxa(nht,oht,dht);
	sf_setlabel(aht,"ht");
	sf_setunit (aht,"");

	/* reflection angle */
	ath = sf_iaxa(Fang,1); sf_setlabel(ath,"th");
	/* azimuth angle */
	aph = sf_iaxa(Fang,2); sf_setlabel(aph,"ph");
	/* CIP axis */
	ac  = sf_iaxa(Fang,3); sf_setlabel(ac ,"c ");

	/* output axes */
	sf_oaxa(Fcip,ahx,1);
	sf_oaxa(Fcip,ahy,2);
	sf_oaxa(Fcip,ahz,3);
	sf_oaxa(Fcip,aht,4);
	sf_oaxa(Fcip,ac ,5);
    }

    if (verb){
	sf_raxa(ahx);
	sf_raxa(ahy);
	sf_raxa(ahz);
	sf_raxa(aht);
	sf_raxa(ac);
	sf_raxa(ath);
	sf_raxa(aph);
    }

    if(anis) {
	/* deviation angle */
	if(! sf_getint  ("nps",&nps)) nps=251;
	if(! sf_getfloat("ops",&ops)) ops=-25;
	if(! sf_getfloat("dps",&dps)) dps=0.2;
	aps = sf_maxa(nps,ops,dps);
	sf_setlabel(aps,"ps");
	sf_setunit (aps,"deg");

	if(verb) sf_raxa(aps);
    } else {
	aps = NULL;
    }

    /*------------------------------------------------------------*/
    /* allocate arrays */
    cip = sf_floatalloc4  (sf_n(ahx),sf_n(ahy),sf_n(ahz),sf_n(aht));
    ang = sf_floatalloc2  (sf_n(ath),sf_n(aph));

    /* read velocity */
    vep = sf_floatalloc  (sf_n(ac));    
    sf_floatread(vep,sf_n(ac),Fvel);

    ves = sf_floatalloc  (sf_n(ac));    
    sf_floatread(ves,sf_n(ac),Fvel);
	
    /*------------------------------------------------------------*/
    /* read normals */
    nn  = (vc3d*) sf_alloc(sf_n(ac),sizeof(*nn)); /* normals  */
    vc3dread1(Fnor,nn,sf_n(ac));

    if(anis) {
	/* read anisotropy */
	eps = sf_floatalloc   (sf_n(ac));
	sf_floatread(eps,sf_n(ac),Fani);

	dlt = sf_floatalloc   (sf_n(ac));
	sf_floatread(dlt,sf_n(ac),Fani);

	/* read tilts */	
	tt  = (vc3d*) sf_alloc(sf_n(ac),sizeof(*tt));
	vc3dread1(Ftlt,tt,sf_n(ac));
    }

    /*------------------------------------------------------------*/
    /* in-plane azimuth reference */
    vv.dx=1;
    vv.dy=0;
    vv.dz=0;

    aa  = (vc3d*) sf_alloc(sf_n(ac),sizeof(*aa));
    for(ic=0;ic<sf_n(ac);ic++) {
	jk    =vcp3d(&nn[ic],&vv);
	aa[ic]=vcp3d(&jk,&nn[ic]);
    }

    /*------------------------------------------------------------*/
    ssn = 1./sqrt(sf_n(ahx)*sf_n(ahy)*sf_n(ahz));

    /*------------------------------------------------------------*/
    /* loop over CIPs */
/*    if(verb) fprintf(stderr,"ic\n");*/
    for(ic=0;ic<sf_n(ac);ic++) {
/*	if(verb) fprintf(stderr,"\b\b\b\b\b%d",ic);*/

	if(adj) {

	    /* read CIP */
	    sf_floatread(cip[0][0][0],sf_n(ahx)*sf_n(ahy)*sf_n(ahz)*sf_n(aht),Fcip);
	    
	    /* init ANG */
	    for    (iph=0;iph<sf_n(aph);iph++) {
		for(ith=0;ith<sf_n(ath);ith++) {
		    ang[iph][ith]=0;
		}
	    }
	} else {
	    
	    /* init CIP */
	    for            (iht=0;iht<sf_n(aht);iht++) {
		for        (ihz=0;ihz<sf_n(ahz);ihz++) {
		    for    (ihy=0;ihy<sf_n(ahy);ihy++) {
			for(ihx=0;ihx<sf_n(ahx);ihx++) {
			    cip[iht][ihz][ihy][ihx]=0;
			}
		    }
		}
	    }

	    /* read ANG */
	    sf_floatread(ang[0],sf_n(ath)*sf_n(aph),Fang);
	}

	/* phi loop */

	nph = sf_n(aph);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)				\
    private(iph,phi,jk,qq,						\
	    ith,tht,							\
	    ihy,ihx,hy,hx,hz,						\
	    tau,jht,fht,cosum,codif,v_s,v_r,psi,sitovel)	\
    shared( nph,aph,ath,aps,ahy,ahx,aht,cip,ang,vep,ves,eps,dlt)
#endif
	for(iph=0;iph<nph;iph++) {
	    phi=(180+sf_o(aph)+iph*sf_d(aph))/180.*SF_PI;
	    /* use '180' to reverse illumination direction: */
	    /* at a CIP, look toward the source */

	    /* reflection azimuth vector */
	    jk = rot3d(nn,aa,phi);
	    qq = nor3d(&jk);
	    
	    /* theta loop */
	    for(ith=0;ith<sf_n(ath);ith++) {
		tht=(sf_o(ath)+ith*sf_d(ath))/180.*SF_PI;
		
		if(anis) {

		    ttipar = psitti(nn,&qq,tt,aa,
				    tht,phi,aps,
				    vep[ic],ves[ic],eps[ic],dlt[ic]);
		    psi = ttipar[0];
		    v_s = ttipar[1];
		    v_r = ttipar[2];

		    psi *= SF_PI/180.;
                    cosum = cosf(tht+psi);
                    codif = cosf(tht-psi);

                    sitovel = sinf(2*tht)/(v_s*cosum + v_r*codif);
		} else {
		    sitovel = sinf(tht)/vep[ic];
		}

		/* lag loops */
		if(adj) {
		    for    (ihy=0;ihy<sf_n(ahy);ihy++) { hy=sf_o(ahy)+ihy*sf_d(ahy);
			for(ihx=0;ihx<sf_n(ahx);ihx++) { hx=sf_o(ahx)+ihx*sf_d(ahx);
 
			    hz = -(hx*(nn[ic].dx)+hy*(nn[ic].dy))/(nn[ic].dz);
			    tau = -((qq.dx)*hx+(qq.dy)*hy+(qq.dz)*hz)*sitovel;			    
			    jht=0.5+(tau-sf_o(aht))/sf_d(aht);

			    if(jht>=0 && jht<sf_n(aht)-1) {
				fht= (tau-sf_o(aht))/sf_d(aht)-jht;
				ang[iph][ith] += (1-fht)*ssn*cip[jht  ][0][ihy][ihx]
				    +               fht *ssn*cip[jht+1][0][ihy][ihx]; 
			    }
			    
			} /* hx */
		    } /* hy */
		} else {
		    for    (ihy=0;ihy<sf_n(ahy);ihy++) { hy=sf_o(ahy)+ihy*sf_d(ahy);
			for(ihx=0;ihx<sf_n(ahx);ihx++) { hx=sf_o(ahx)+ihx*sf_d(ahx);
		    
			    hz = -(hx*(nn[ic].dx)+hy*(nn[ic].dx))/(nn[ic].dz);
			    tau = -((qq.dx)*hx+(qq.dy)*hy+(qq.dz)*hz)*sitovel; 		    
			    jht=0.5+(tau-sf_o(aht))/sf_d(aht);

			    if(jht>=0 && jht<sf_n(aht)-1) {
				fht= (tau-sf_o(aht))/sf_d(aht)-jht;
				cip[jht  ][0][ihy][ihx] += (1-fht)*ssn*ang[iph][ith];
				cip[jht+1][0][ihy][ihx] +=    fht *ssn*ang[iph][ith];
			    }
			    
			} /* hx */
		    } /* hy */
		}

	    } /* th */
	} /* ph */
	
	if(adj) {
	    /* write ANG */
	    sf_floatwrite(ang[0],sf_n(ath)*sf_n(aph),Fang);
	} else {
	    /* write CIP */
	    sf_floatwrite(cip[0][0][0],sf_n(ahx)*sf_n(ahy)*sf_n(ahz)*sf_n(aht),Fcip);
	}

    }
    if(verb) fprintf(stderr,"\n");
    /*------------------------------------------------------------*/
 
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"free memory...");
    free(***cip);free(**cip);free(*cip);free(cip);
    ;                        free(*ang);free(ang);
    ;                                   free(vep);
    ;                                   free (nn);
    ;                                   free (aa);
    if(anis) {
	free(ves);
	free(eps);
	free(dlt);
	free(tt);
    }
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    exit(0);
}		
