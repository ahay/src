/* linearized acoustic time-domain FD modeling */
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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fdutil.h"

/* check: dt<= 0.2 * min(dx,dz)/vmin */

#define NOP 2 /* derivative operator half-size */

#define C0 -2.500000 /*    c0=-30./12.; */
#define CA +1.333333 /*    ca=+16./12.; */
#define CB -0.083333 /*    cb=- 1./12.; */

#define C1  0.66666666666666666666 /*  2/3  */	
#define C2 -0.08333333333333333333 /* -1/12 */

/* centered FD derivative stencils */
#define D1(a,i2,i1,s) (C2*(a[i2  ][i1+2] - a[i2  ][i1-2]) +  \
                       C1*(a[i2  ][i1+1] - a[i2  ][i1-1])  )*s
#define D2(a,i2,i1,s) (C2*(a[i2+2][i1  ] - a[i2-2][i1  ]) +  \
                       C1*(a[i2+1][i1  ] - a[i2-1][i1  ])  )*s

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl; 
    int  jsnap,ntsnap;
    int  jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */

    sf_file Fvel=NULL; /* velocity  */
    sf_file Fref=NULL; /* reflectivity */
    sf_file Fden=NULL; /* density   */

    sf_file Fdat=NULL; /* data (background)      */
    sf_file Fwfl=NULL; /* wavefield (background) */

    sf_file Flid=NULL; /* data (scattered)      */
    sf_file Fliw=NULL; /* wavefield (scattered) */

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */

    float **vpin=NULL;         /* velocity  */
    float **roin=NULL;         /* density   */
    float **rfin=NULL;         /* reflectivity */

    float **vp=NULL;           /* velocity     in expanded domain */
    float **ro=NULL;           /* density      in expanded domain */
    float **ro1=NULL;          /* normalized 1st derivative of density on axis 1 */
    float **ro2=NULL;          /* normalized 1st derivative of density on axis 2 */

    float **rf=NULL;           /* reflectivity in expanded domain */

    float  *bdd=NULL;          /* data (background) */
    float  *sdd=NULL;          /* data (scattered)  */

    float **vt=NULL;           /* temporary vp*vp * dt*dt */

    float **bum,**buo,**bup,**bua,**but; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **sum,**suo,**sup,**sua,**sut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* cube axes */
    sf_axis at,a1,a2,as,ar;
    int     nt,n1,n2,ns,nr,nb;
    int     it,i1,i2;
    float   dt,d1,d2,id1,id2,dt2;

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    fdm2d    fdm;
    abcone2d abc;     /* abc */
    sponge spo;

    /* FD operator size */
    float co,ca2,cb2,ca1,cb1;

    int ompchunk; 
#ifdef _OPENMP
    int ompnth,ompath;
#endif

    sf_axis   ac1=NULL,ac2=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  
    /* OpenMP data chunk size */
#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;  
    /* OpenMP available threads */

#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* "exploding reflector" */

    Fwav = sf_input ("in" ); /* wavelet   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */

    Fvel = sf_input ("vel"); /* velocity  */
    Fden = sf_input ("den"); /* density   */
    Fref = sf_input ("ref"); /* reflectivity */

    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    Fliw = sf_output("liw"); /* wavefield (scattered) */
    Flid = sf_output("lid"); /* data (scattered) */

    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    a1 = sf_iaxa(Fvel,1); sf_setlabel(a1,"z"); if(verb) sf_raxa(a1); /* depth */
    a2 = sf_iaxa(Fvel,2); sf_setlabel(a2,"x"); if(verb) sf_raxa(a2); /* space */

    nt = sf_n(at); dt = sf_d(at);
    ns = sf_n(as);
    nr = sf_n(ar);
    n1 = sf_n(a1); d1 = sf_d(a1);
    n2 = sf_n(a2); d2 = sf_d(a2);

    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }


    /*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Flid,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);
    sf_oaxa(Flid,at,2);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(a1);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(a2);
	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(a1);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(a2);
	dqz=sf_d(a1);
	dqx=sf_d(a2);

	ac1 = sf_maxa(nqz,oqz,dqz);
	ac2 = sf_maxa(nqx,oqx,dqx);

	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(ac1),sf_n(ac2));

	ntsnap=0;
        for(it=0; it<nt; it++) {
            if(it%jsnap==0) ntsnap++;
        }
        sf_setn(at,  ntsnap);
        sf_setd(at,dt*jsnap);
        if(verb) sf_raxa(at);

/*	sf_setn(at,nt/jsnap);
	sf_setd(at,dt*jsnap); */

	sf_oaxa(Fwfl,ac1,1);
	sf_oaxa(Fwfl,ac2,2);
	sf_oaxa(Fwfl,at, 3);

	sf_oaxa(Fliw,ac1,1);
	sf_oaxa(Fliw,ac2,2);
	sf_oaxa(Fliw,at, 3);
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,a1,a2,nb,ompchunk);

    sf_setn(a1,fdm->nzpad); sf_seto(a1,fdm->ozpad); if(verb) sf_raxa(a1);
    sf_setn(a2,fdm->nxpad); sf_seto(a2,fdm->oxpad); if(verb) sf_raxa(a2);

    /*------------------------------------------------------------*/
    if(expl) ww = sf_floatalloc( 1);
    else     ww = sf_floatalloc(ns);
    bdd =sf_floatalloc(nr);
    sdd =sf_floatalloc(nr);

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
    dt2 =    dt*dt;
    id1 = 1/d1;
    id2 = 1/d2;

    co = C0 * (id2*id2+id1*id1);
    ca2= CA *  id2*id2;
    cb2= CB *  id2*id2;
    ca1= CA *          id1*id1;
    cb1= CB *          id1*id1;

    /*------------------------------------------------------------*/ 
    /* input density */
    roin=sf_floatalloc2(n1,   n2   ); 
    ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    ro1 =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ro2 =sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    sf_floatread(roin[0],n1*n2,Fden); 
    expand(roin,ro,fdm);

    /* normalized density derivatives */
    for    (i2=NOP; i2<fdm->nxpad-NOP; i2++) {
	for(i1=NOP; i1<fdm->nzpad-NOP; i1++) {
	    ro1[i2][i1] = D1(ro,i2,i1,id1) / ro[i2][i1];
	    ro2[i2][i1] = D2(ro,i2,i1,id2) / ro[i2][i1];
	}
    }

    free(*roin); free(roin);

    /*------------------------------------------------------------*/
    /* input velocity */
    vpin=sf_floatalloc2(n1,   n2   ); 
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(vpin[0],n1*n2,Fvel);
    expand(vpin,vp,fdm);
    free(*vpin); free(vpin);

    /*------------------------------------------------------------*/
    /* input reflectivity */
    rfin=sf_floatalloc2(n1,   n2   ); 
    rf  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(rfin[0],n1*n2,Fref); 
    expand(rfin,rf,fdm);
    free(*rfin); free(rfin);

    for    (i2=0; i2<fdm->nxpad; i2++) {
	for(i1=0; i1<fdm->nzpad; i1++) {
	    vt[i2][i1] = vp[i2][i1] * vp[i2][i1] * dt2;
	}
    }

    /* free surface */
    if(fsrf) {
	for    (i2=0; i2<fdm->nxpad; i2++) {
	    for(i1=0; i1<fdm->nb; i1++) {
		vt[i2][i1]=0;
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    bum=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    buo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    bup=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    bua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    sum=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    suo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    sup=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    sua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (i2=0; i2<fdm->nxpad; i2++) {
	for(i1=0; i1<fdm->nzpad; i1++) {
	    bum[i2][i1]=0;
	    buo[i2][i1]=0;
	    bup[i2][i1]=0;
	    bua[i2][i1]=0;

	    sum[i2][i1]=0;
	    suo[i2][i1]=0;
	    sup[i2][i1]=0;
	    sua[i2][i1]=0;
	}
    }

    /*------------------------------------------------------------*/
    /* one-way abc setup */
    abc = abcone2d_make(NOP,dt,vp,fsrf,fdm);
    /* sponge abc setup */
    spo = sponge_make(fdm->nb);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i2,i1) shared(fdm,bua,buo,sua,suo,co,ca2,ca1,cb2,cb1,id2,id1)
#endif
	for    (i2=NOP; i2<fdm->nxpad-NOP; i2++) {
	    for(i1=NOP; i1<fdm->nzpad-NOP; i1++) {
		
		/* 4th order Laplacian operator */
		bua[i2][i1] = 
		    co * buo[i2  ][i1  ] + 
		    ca2*(buo[i2-1][i1  ] + buo[i2+1][i1  ]) +
		    cb2*(buo[i2-2][i1  ] + buo[i2+2][i1  ]) +
		    ca1*(buo[i2  ][i1-1] + buo[i2  ][i1+1]) +
		    cb1*(buo[i2  ][i1-2] + buo[i2  ][i1+2]);
		sua[i2][i1] = 
		    co * suo[i2  ][i1  ] + 
		    ca2*(suo[i2-1][i1  ] + suo[i2+1][i1  ]) +
		    cb2*(suo[i2-2][i1  ] + suo[i2+2][i1  ]) +
		    ca1*(suo[i2  ][i1-1] + suo[i2  ][i1+1]) +
		    cb1*(suo[i2  ][i1-2] + suo[i2  ][i1+2]);
		
		/* density term */
		bua[i2][i1] -= (
		    D1(buo,i2,i1,id1) * ro1[i2][i1] +
		    D2(buo,i2,i1,id2) * ro2[i2][i1] );
		sua[i2][i1] -= (
		    D1(suo,i2,i1,id1) * ro1[i2][i1] +
		    D2(suo,i2,i1,id2) * ro2[i2][i1] );
	    }
	}   
	
	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    lint2d_inject1(bua,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    lint2d_inject(bua,ww,cs);
	}

	/* single scattering */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(i2,i1) shared(fdm,buo,sua,rf)
#endif
	for     (i2=0; i2<fdm->nxpad; i2++) {
	    for (i1=0; i1<fdm->nzpad; i1++) {
		sua[i2][i1] -= bua[i2][i1] * 2*rf[i2][i1];
	    }
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(i2,i1) shared(fdm,bua,buo,bum,bup,sua,suo,sum,sup,vt,dt2)
#endif
	for    (i2=0; i2<fdm->nxpad; i2++) {
	    for(i1=0; i1<fdm->nzpad; i1++) {
		bup[i2][i1] = 2*buo[i2][i1] 
		    -           bum[i2][i1] 
		    +           bua[i2][i1] * vt[i2][i1];

		sup[i2][i1] = 2*suo[i2][i1] 
		    -           sum[i2][i1] 
		    +           sua[i2][i1] * vt[i2][i1];

	    }
	}
	/* circulate wavefield arrays */
	but=bum;
	bum=buo;
	buo=bup;
	bup=but;

	sut=sum;
	sum=suo;
	suo=sup;
	sup=sut;
	
	/* one-way abc apply*/
	abcone2d_apply(buo,bum,NOP,abc,fdm);
	sponge2d_apply(bum,        spo,fdm);
	sponge2d_apply(buo,        spo,fdm);

	abcone2d_apply(suo,sum,NOP,abc,fdm);
	sponge2d_apply(sum,        spo,fdm);
	sponge2d_apply(suo,        spo,fdm);

	/* extract data at receivers */
	lint2d_extract(buo,bdd,cr);
	lint2d_extract(suo,sdd,cr);
	if(        it%jdata==0) {
	    sf_floatwrite(bdd,nr,Fdat);
	    sf_floatwrite(sdd,nr,Flid);
	}

	/* extract wavefield in the "box" */
	if(snap && it%jsnap==0) {
	    cut2d(buo,uc,fdm,ac1,ac2);
	    sf_floatwrite(uc[0],sf_n(ac1)*sf_n(ac2),Fwfl);

	    cut2d(suo,uc,fdm,ac1,ac2);
	    sf_floatwrite(uc[0],sf_n(ac1)*sf_n(ac2),Fliw);
	}

    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
