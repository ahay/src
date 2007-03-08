/* Elastic time-domain FD modeling */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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
#include "fdutil.h"

#define NOP 4 /* derivative operator half-size */

/* Muir's derivative operator */
#define C1 +0.598144 // 1225/ 1024      /2
#define C2 -0.039876 //-1225/(1024*  15)/2
#define C3 +0.004785 // 1225/(1024* 125)/2
#define C4 -0.000348 //-1225/(1024*1715)/2

/*  forward FD derivative stencils */
#define F1(a,i2,i1,s) (C4*(a[i2  ][i1+4] - a[i2  ][i1-3]) +	\
		       C3*(a[i2  ][i1+3] - a[i2  ][i1-2]) +	\
		       C2*(a[i2  ][i1+2] - a[i2  ][i1-1]) +	\
                       C1*(a[i2  ][i1+1] - a[i2  ][i1  ])  )*s
#define F2(a,i2,i1,s) (C4*(a[i2+4][i1  ] - a[i2-3][i1  ]) +	\
		       C3*(a[i2+3][i1  ] - a[i2-2][i1  ]) +	\
		       C2*(a[i2+2][i1  ] - a[i2-1][i1  ]) +	\
                       C1*(a[i2+1][i1  ] - a[i2  ][i1  ])  )*s

/* backward FD derivative stencils */
#define B1(a,i2,i1,s) (C4*(a[i2  ][i1+3] - a[i2  ][i1-4]) +	\
		       C3*(a[i2  ][i1+2] - a[i2  ][i1-3]) +	\
		       C2*(a[i2  ][i1+1] - a[i2  ][i1-2]) +	\
                       C1*(a[i2  ][i1  ] - a[i2  ][i1-1])  )*s
#define B2(a,i2,i1,s) (C4*(a[i2+3][i1  ] - a[i2-4][i1  ]) +	\
		       C3*(a[i2+2][i1  ] - a[i2-3][i1  ]) +	\
		       C2*(a[i2+1][i1  ] - a[i2-2][i1  ]) +	\
                       C1*(a[i2  ][i1  ] - a[i2-1][i1  ])  )*s

int main(int argc, char* argv[])
{
    bool verb,free,snap,ssou,opot;
    int  jsnap;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    float **dd=NULL;           /* data      */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */

    float **roin=NULL;         /* density   */
    float **ro=NULL;           /* density  in expanded domain */

    float **c11in=NULL;        /* stiffness */
    float **c13in=NULL;
    float **c33in=NULL;
    float **c44in=NULL;

    float **c11=NULL;          /* stiffness in expanded domain */
    float **c13=NULL;
    float **c33=NULL;
    float **c44=NULL;

    float **um1,**uo1,**up1,**ua1,**ut1; /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **um2,**uo2,**up2,**ua2,**ut2;

    float **t11,**t12,**t22;       /* stress/strain tensor */ 
    float   s11,  s12,  s22;

    float **qp,**qs;               /* potential (P waves, S waves)  */
    float **vp,**vs;               /* potential (P waves, S waves)  */

    /* cube axes */
    sf_axis at,a1,a2,as,ar,ac;
    int     nt,n1,n2,ns,nr,nc,nb;
    int     it,i1,i2;
    float   dt,d1,d2,id1,id2,dt2;

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    fdm2d    fdm;
    abcone2d abcp,abcs;     /* abc */
    sponge2d spo;

    int ompchunk;
    int nbell;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&free)) free=false; /* free surface flag */
    if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
    if(! sf_getbool("opot",&opot)) opot=false; /* output potential */

    if(! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
    sf_warning("nbell=%d",nbell);

    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    Fden = sf_input ("den"); /* density   */

    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    a1 = sf_iaxa(Fccc,1); sf_setlabel(a1,"z"); if(verb) sf_raxa(a1); /* depth */
    a2 = sf_iaxa(Fccc,2); sf_setlabel(a2,"x"); if(verb) sf_raxa(a2); /* space */

    /* 2D vector components */
    nc=2;
    ac=sf_maxa(nc,0,1);

    nt = sf_n(at); dt = sf_d(at);
    ns = sf_n(as);
    nr = sf_n(ar);
    n1 = sf_n(a1); d1 = sf_d(a1);
    n2 = sf_n(a2); d2 = sf_d(a2);

    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,free,a1,a2,nb,ompchunk);
    fdbell_init(nbell);

    sf_setn(a1,fdm->n1pad); sf_seto(a1,fdm->o1pad); if(verb) sf_raxa(a1);
    sf_setn(a2,fdm->n2pad); sf_seto(a2,fdm->o2pad); if(verb) sf_raxa(a2);
    /*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header */
    if(snap) {
	sf_setn(at,nt/jsnap);
	sf_setd(at,dt*jsnap);

	sf_oaxa(Fwfl,a1,1);
	sf_oaxa(Fwfl,a2,2);
	sf_oaxa(Fwfl,ac,3);
	sf_oaxa(Fwfl,at,4);
    }

    /* source array */
    ww=sf_floatalloc3(ns,nc,nt); 
    sf_floatread(ww[0][0],nt*nc*ns,Fwav);

    /* data array */
    dd=sf_floatalloc2(nr,nc);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    dt2 = dt*dt;
    id1 = 2/d1;
    id2 = 2/d2;

    /*------------------------------------------------------------*/ 
    /* setup model */
    roin=sf_floatalloc2(n1   ,n2   ); 
    ro  =sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    sf_floatread(roin[0],n1*n2,Fden); 
    expand(roin,ro,fdm);

    c11in=sf_floatalloc2(n1   ,n2   ); 
    c13in=sf_floatalloc2(n1   ,n2   ); 
    c33in=sf_floatalloc2(n1   ,n2   ); 
    c44in=sf_floatalloc2(n1   ,n2   ); 

    c11  =sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    c13  =sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    c33  =sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    c44  =sf_floatalloc2(fdm->n1pad,fdm->n2pad); 

    sf_floatread(c11in[0],n1*n2,Fccc );
    sf_floatread(c13in[0],n1*n2,Fccc );
    sf_floatread(c33in[0],n1*n2,Fccc );
    sf_floatread(c44in[0],n1*n2,Fccc );

    expand(c11in,c11,fdm);
    expand(c13in,c13,fdm);
    expand(c33in,c33,fdm);
    expand(c44in,c44,fdm);

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    um1=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    uo1=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    up1=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    ua1=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    um2=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    uo2=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    up2=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    ua2=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    qp=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    qs=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    t11=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    t12=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    t22=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    for    (i2=0; i2<fdm->n2pad; i2++) {
	for(i1=0; i1<fdm->n1pad; i1++) {
	    um1[i2][i1]=0; um2[i2][i1]=0;
	    uo1[i2][i1]=0; uo2[i2][i1]=0;
	    up1[i2][i1]=0; up2[i2][i1]=0;
	    ua1[i2][i1]=0; ua2[i2][i1]=0;
	}
    }

    /* one-way abc setup   */
    vp = sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    vs = sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    for    (i2=0; i2<fdm->n2pad; i2++) {
	for(i1=0; i1<fdm->n1pad; i1++) {
	    vp[i2][i1] = sqrt( c11[i2][i1]/ro[i2][i1] );
	    vs[i2][i1] = sqrt( c13[i2][i1]/ro[i2][i1] );
	}
    }
    abcp = abcone2d_make(NOP,dt,vp,free,fdm);
    abcs = abcone2d_make(NOP,dt,vs,free,fdm);

    /* sponge abc setup */
    spo = sponge2d_make(fdm);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

	/*------------------------------------------------------------*/
	/* from displacement to strain                                */
	/*------------------------------------------------------------*/
	/* ezz(x     ,z     ) = Fz( uz(x     ,z-dz/2) )
	   exx(x     ,z     ) = Fx( ux(x-dx/2,z     ) ) 
	   exz(x-dx/2,z-dz/2) = Bx( uz(x     ,z-dz/2) ) + */
	/*                      Bz( ux(x-dx/2,z     ) )   */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i1,i2) shared(fdm,t11,t12,t22,uo1,uo2,id1,id2)
#endif
	for    (i2=NOP; i2<fdm->n2pad-NOP; i2++) {
	    for(i1=NOP; i1<fdm->n1pad-NOP; i1++) {
		
		t11[i2][i1] = F1(uo1,i2,i1,id1);
		
		t22[i2][i1] = F2(uo2,i2,i1,id2);

	        t12[i2][i1] = B2(uo1,i2,i1,id2) 
		    +         B1(uo2,i2,i1,id1);
	    }
	}		

	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
	/* szz(x     ,z     ) = c11(x     ,z     ) ezz(x     ,z     ) + */
	/*                      c13(x     ,z     ) exx(x     ,z     )   */
	/* sxx(x     ,z     ) = c13(x     ,z     ) ezz(x     ,z     ) + */
	/*                      c33(x     ,z     ) exx(x     ,z     )   */
	/* sxz(x-dx/2,z-dz/2) = c44(x-dx/2,z-dz/2) exz(x-dx/2,z-dz/2)   */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i1,i2,s11,s12,s22) shared(fdm,t11,t12,t22,c11,c13,c33,c44)
#endif
	for    (i2=0; i2<fdm->n2pad; i2++) {
	    for(i1=0; i1<fdm->n1pad; i1++) {
		
		s11 = c33[i2][i1] * t11[i2][i1] 
		    + c13[i2][i1] * t22[i2][i1];
		
		s22 = c13[i2][i1] * t11[i2][i1] 
		    + c11[i2][i1] * t22[i2][i1];
		
		s12 = c44[i2][i1] * t12[i2][i1];

		t11[i2][i1] = s11;
		t22[i2][i1] = s22;
		t12[i2][i1] = s12;
	    }
	}

	if(free) {
	    for    (i2=0; i2<fdm->n2pad; i2++) {
		for(i1=0; i1<fdm->nb; i1++) {
		    t11[i2][i1]=0;
		    t12[i2][i1]=0;
		    t22[i2][i1]=0;
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(ssou) {
/*	    lint2d_inject(t11,ww[it][0],cs);*/
/*	    lint2d_inject(t22,ww[it][0],cs);*/
	    lint2d_bell(t11,ww[it][0],cs);
	    lint2d_bell(t22,ww[it][0],cs);
	}

	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* az(x,z-dz/2) = Fx( sxz(x-dx/2,z-dz/2) ) + */
	/*                Bz( szz(x     ,z     ) )   */
	/* ax(x-dx/2,z) = Bx( sxx(x     ,z     ) ) + */
	/*                Fz( sxz(x-dx/2,z-dz/2) )   */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i1,i2) shared(fdm,t11,t12,t22,ua1,ua2,id1,id2)
#endif
	for    (i2=NOP; i2<fdm->n2pad-NOP; i2++) {
	    for(i1=NOP; i1<fdm->n1pad-NOP; i1++) {

		ua1[i2][i1] = F2( t12,i2,i1,id2 ) 
		    +         B1( t11,i2,i1,id1 );

		ua2[i2][i1] = B2( t22,i2,i1,id2 ) 
		    +         F1( t12,i2,i1,id1 );
	    }
	}

	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(!ssou) {
/*	    lint2d_inject(ua1,ww[it][0],cs);*/
/*	    lint2d_inject(ua2,ww[it][1],cs);*/
	    lint2d_bell(ua1,ww[it][0],cs);
	    lint2d_bell(ua2,ww[it][1],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i1,i2) shared(fdm,uo1,uo2,um1,um2,up1,up2,ua1,ua2,ro,dt2)
#endif
	for    (i2=0; i2<fdm->n2pad; i2++) {
	    for(i1=0; i1<fdm->n1pad; i1++) {
		up1[i2][i1] = 2*uo1[i2][i1] 
		    -           um1[i2][i1] 
		    +           ua1[i2][i1] / ro[i2][i1] * dt2; 

		up2[i2][i1] = 2*uo2[i2][i1] 
		    -           um2[i2][i1] 
		    +           ua2[i2][i1] / ro[i2][i1] * dt2; 
	    }
	}
	/* circulate wavefield arrays */
	ut1=um1; ut2=um2;
	um1=uo1; um2=uo2;
	uo1=up1; uo2=up2;
	up1=ut1; up2=ut2;

	/*------------------------------------------------------------*/
	/* one-way abc                                                */
	/*------------------------------------------------------------*/
	abcone2d_apply(uo1,um1,NOP,abcp,fdm);
	abcone2d_apply(uo2,um2,NOP,abcp,fdm);

	abcone2d_apply(uo1,um1,NOP,abcs,fdm);
	abcone2d_apply(uo2,um2,NOP,abcs,fdm);

	/*------------------------------------------------------------*/
	/* sponge abc                                                 */
	/*------------------------------------------------------------*/
	sponge2d_apply(um1,spo,fdm);
	sponge2d_apply(uo1,spo,fdm);
	sponge2d_apply(up1,spo,fdm);

	sponge2d_apply(um2,spo,fdm);
	sponge2d_apply(uo2,spo,fdm);
	sponge2d_apply(up2,spo,fdm);

	/*------------------------------------------------------------*/
	/* compute potentials                                         */
	/*------------------------------------------------------------*/
	/* qp(x     ,z     ) = Fx( ux(x-dx/2,z     ) ) + */
	/*                     Fz( uz(x     ,z-dz/2) )   */
	/* qs(x-dx/2,z-dz/2) = Bz( ux(x-dx/2,z     ) ) + */
	/*                     Bz( ux(x     ,z-dz/2) )   */

	if(opot) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i1,i2) shared(fdm,uo1,uo2,qp,qs,id1,id2)
#endif
	    for    (i2=NOP; i2<fdm->n2pad-NOP; i2++) {
		for(i1=NOP; i1<fdm->n1pad-NOP; i1++) {
		    
		    qp[i2][i1] = F1( uo1,i2,i1,id1 ) 
			+        F2( uo2,i2,i1,id2 );
		    
		    qs[i2][i1] = B1( uo2,i2,i1,id1 ) 
			-        B2( uo1,i2,i1,id2 );
		}
	    }

	    lint2d_extract(qp,dd[0],cr);
	    lint2d_extract(qs,dd[1],cr);
	    
	    if(snap && it%jsnap==0) {
		sf_floatwrite(qp[0],fdm->n1pad*fdm->n2pad,Fwfl);
		sf_floatwrite(qs[0],fdm->n1pad*fdm->n2pad,Fwfl);
	    }
	} else {

	    lint2d_extract(uo1,dd[0],cr);
	    lint2d_extract(uo2,dd[1],cr);
	    
	    if(snap && it%jsnap==0) {
		sf_floatwrite(uo1[0],fdm->n1pad*fdm->n2pad,Fwfl);
		sf_floatwrite(uo2[0],fdm->n1pad*fdm->n2pad,Fwfl);
	    }
	}
	sf_floatwrite(dd[0],nr*nc,Fdat);

    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
