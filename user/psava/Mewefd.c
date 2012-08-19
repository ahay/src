/* elastic time-domain FD modeling */
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

#define NOP 4 /* derivative operator half-size */

/* Muir's derivative operator */
#define C1 +0.598144 /*  1225/ 1024      /2 */
#define C2 -0.039876 /* -1225/(1024*  15)/2 */
#define C3 +0.004785 /*  1225/(1024* 125)/2 */
#define C4 -0.000348 /* -1225/(1024*1715)/2 */

/*  forward FD derivative stencils */
#define Fz(a,ix,iz,s) (C4*(a[ix  ][iz+4] - a[ix  ][iz-3]) +	\
		       C3*(a[ix  ][iz+3] - a[ix  ][iz-2]) +	\
		       C2*(a[ix  ][iz+2] - a[ix  ][iz-1]) +	\
                       C1*(a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Fx(a,ix,iz,s) (C4*(a[ix+4][iz  ] - a[ix-3][iz  ]) +	\
		       C3*(a[ix+3][iz  ] - a[ix-2][iz  ]) +	\
		       C2*(a[ix+2][iz  ] - a[ix-1][iz  ]) +	\
                       C1*(a[ix+1][iz  ] - a[ix  ][iz  ])  )*s

/* backward FD derivative stencils */
#define Bz(a,ix,iz,s) (C4*(a[ix  ][iz+3] - a[ix  ][iz-4]) +	\
		       C3*(a[ix  ][iz+2] - a[ix  ][iz-3]) +	\
		       C2*(a[ix  ][iz+1] - a[ix  ][iz-2]) +	\
                       C1*(a[ix  ][iz  ] - a[ix  ][iz-1])  )*s
#define Bx(a,ix,iz,s) (C4*(a[ix+3][iz  ] - a[ix-4][iz  ]) +	\
		       C3*(a[ix+2][iz  ] - a[ix-3][iz  ]) +	\
		       C2*(a[ix+1][iz  ] - a[ix-2][iz  ]) +	\
                       C1*(a[ix  ][iz  ] - a[ix-1][iz  ])  )*s

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,ssou,opot;
    int  jsnap;
    int  jdata;

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

    float **umz,**uoz,**upz,**ua1,**utz; /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **umx,**uox,**upx,**uax,**utx;

    float **tzz,**txz,**txx;       /* stress/strain tensor */ 
    float   szz,  sxz,  sxx;

    float **qp=NULL,**qs=NULL;     /* potential (P waves, S waves) */
    float **vp,     **vs;          /* velocity  (P waves, S waves) */

    /* cube axes */
    sf_axis at,az,ax,as,ar,ac;
    int     nt,n1,n2,ns,nr,nc,nb;
    int     it,iz,ix;
    float   dt,d1,d2,idz,idx,dt2;

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    fdm2d    fdm;
    abcone2d abcp,abcs;     /* abc */
    sponge2d spo;

    int nbell;

    int ompchunk;
#ifdef _OPENMP
    int ompnth,ompath;
#endif 

    sf_axis   acz=NULL,acx=NULL;
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
    if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
    if(! sf_getbool("opot",&opot)) opot=false; /* output potential */

    if(! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
    sf_warning("nbell=%d",nbell);

    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    Fden = sf_input ("den"); /* density   */

    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    /* 2D vector components */
    nc=2;
    ac=sf_maxa(nc,0,1);

    nt = sf_n(at); dt = sf_d(at);
    ns = sf_n(as);
    nr = sf_n(ar);
    n1 = sf_n(az); d1 = sf_d(az);
    n2 = sf_n(ax); d2 = sf_d(ax);

    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,ompchunk);
    fdbell_init(nbell);

    sf_setn(az,fdm->n1pad); sf_seto(az,fdm->o1pad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->n2pad); sf_seto(ax,fdm->o2pad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);
	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
	dqz=sf_d(az);
	dqx=sf_d(ax);

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);

	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

	sf_setn(at,nt/jsnap);
	sf_setd(at,dt*jsnap);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,ac, 3);
	sf_oaxa(Fwfl,at, 4);
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
    idz = 2/d1;
    idx = 2/d2;

    /*------------------------------------------------------------*/ 
    /* input density */
    roin=sf_floatalloc2(n1   ,n2   ); 
    ro  =sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    sf_floatread(roin[0],n1*n2,Fden); 
    expand(roin,ro,fdm);
    free(*roin); free(roin);

    /*------------------------------------------------------------*/
    /* input stiffness */
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

    free(*c11in); free(c11in);
    free(*c13in); free(c13in);
    free(*c33in); free(c33in);
    free(*c44in); free(c44in);

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    umz=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    uoz=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    upz=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    ua1=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    umx=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    uox=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    upx=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    uax=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    if(opot) {
	qp=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
	qs=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    }

    tzz=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    txz=sf_floatalloc2(fdm->n1pad,fdm->n2pad);
    txx=sf_floatalloc2(fdm->n1pad,fdm->n2pad);

    for    (ix=0; ix<fdm->n2pad; ix++) {
	for(iz=0; iz<fdm->n1pad; iz++) {
	    umz[ix][iz]=0; umx[ix][iz]=0;
	    uoz[ix][iz]=0; uox[ix][iz]=0;
	    upz[ix][iz]=0; upx[ix][iz]=0;
	    ua1[ix][iz]=0; uax[ix][iz]=0;
	}
    }

    /* one-way abc setup   */
    vp = sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    vs = sf_floatalloc2(fdm->n1pad,fdm->n2pad); 
    for    (ix=0; ix<fdm->n2pad; ix++) {
	for(iz=0; iz<fdm->n1pad; iz++) {
	    vp[ix][iz] = sqrt( c11[ix][iz]/ro[ix][iz] );
	    vs[ix][iz] = sqrt( c13[ix][iz]/ro[ix][iz] );
	}
    }
    abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
    abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);

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
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(iz,ix) shared(fdm,tzz,txz,txx,uoz,uox,idz,idx)
#endif
	for    (ix=NOP; ix<fdm->n2pad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->n1pad-NOP; iz++) {
		
		tzz[ix][iz] = Fz(uoz,ix,iz,idz);
		
		txx[ix][iz] = Fx(uox,ix,iz,idx);

	        txz[ix][iz] = Bx(uoz,ix,iz,idx) 
		    +         Bz(uox,ix,iz,idz);
	    }
	}		

	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
	/* szz(x     ,z     ) = c33(x     ,z     ) ezz(x     ,z     ) + */
	/*                      c13(x     ,z     ) exx(x     ,z     )   */
	/* sxx(x     ,z     ) = c13(x     ,z     ) ezz(x     ,z     ) + */
	/*                      c11(x     ,z     ) exx(x     ,z     )   */
	/* sxz(x-dx/2,z-dz/2) = c44(x-dx/2,z-dz/2) exz(x-dx/2,z-dz/2)   */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(iz,ix,szz,sxz,sxx) shared(fdm,tzz,txz,txx,c11,c13,c33,c44)
#endif
	for    (ix=0; ix<fdm->n2pad; ix++) {
	    for(iz=0; iz<fdm->n1pad; iz++) {
		
		szz = c33[ix][iz] * tzz[ix][iz] 
		    + c13[ix][iz] * txx[ix][iz];
		
		sxx = c13[ix][iz] * tzz[ix][iz] 
		    + c11[ix][iz] * txx[ix][iz];
		
		sxz = c44[ix][iz] * txz[ix][iz];

		tzz[ix][iz] = szz;
		txx[ix][iz] = sxx;
		txz[ix][iz] = sxz;
	    }
	}

	if(fsrf) {
	    for    (ix=0; ix<fdm->n2pad; ix++) {
		for(iz=0; iz<fdm->nb; iz++) {
		    tzz[ix][iz]=0;
		    txz[ix][iz]=0;
		    txx[ix][iz]=0;
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(ssou) {
/*	    lint2d_inject(tzz,ww[it][0],cs);*/
/*	    lint2d_inject(txx,ww[it][0],cs);*/
	    lint2d_bell(tzz,ww[it][0],cs);
	    lint2d_bell(txx,ww[it][0],cs);
	}

	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* az(x,z-dz/2) = Fx( sxz(x-dx/2,z-dz/2) ) + */
	/*                Bz( szz(x     ,z     ) )   */
	/* ax(x-dx/2,z) = Bx( sxx(x     ,z     ) ) + */
	/*                Fz( sxz(x-dx/2,z-dz/2) )   */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(iz,ix) shared(fdm,tzz,txz,txx,ua1,uax,idz,idx)
#endif
	for    (ix=NOP; ix<fdm->n2pad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->n1pad-NOP; iz++) {

		ua1[ix][iz] = Fx( txz,ix,iz,idx ) 
		    +         Bz( tzz,ix,iz,idz );

		uax[ix][iz] = Bx( txx,ix,iz,idx ) 
		    +         Fz( txz,ix,iz,idz );
	    }
	}

	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(!ssou) {
/*	    lint2d_inject(ua1,ww[it][0],cs);*/
/*	    lint2d_inject(uax,ww[it][1],cs);*/
	    lint2d_bell(ua1,ww[it][0],cs);
	    lint2d_bell(uax,ww[it][1],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(iz,ix) shared(fdm,uoz,uox,umz,umx,upz,upx,ua1,uax,ro,dt2)
#endif
	for    (ix=0; ix<fdm->n2pad; ix++) {
	    for(iz=0; iz<fdm->n1pad; iz++) {
		upz[ix][iz] = 2*uoz[ix][iz] 
		    -           umz[ix][iz] 
		    +           ua1[ix][iz] / ro[ix][iz] * dt2; 

		upx[ix][iz] = 2*uox[ix][iz] 
		    -           umx[ix][iz] 
		    +           uax[ix][iz] / ro[ix][iz] * dt2; 
	    }
	}
	/* circulate wavefield arrays */
	utz=umz; utx=umx;
	umz=uoz; umx=uox;
	uoz=upz; uox=upx;
	upz=utz; upx=utx;

	/*------------------------------------------------------------*/
	/* one-way abc                                                */
	/*------------------------------------------------------------*/
	abcone2d_apply(uoz,umz,NOP,abcp,fdm);
	abcone2d_apply(uox,umx,NOP,abcp,fdm);

	abcone2d_apply(uoz,umz,NOP,abcs,fdm);
	abcone2d_apply(uox,umx,NOP,abcs,fdm);

	/*------------------------------------------------------------*/
	/* sponge abc                                                 */
	/*------------------------------------------------------------*/
	sponge2d_apply(umz,spo,fdm);
	sponge2d_apply(uoz,spo,fdm);

	sponge2d_apply(umx,spo,fdm);
	sponge2d_apply(uox,spo,fdm);

	/*------------------------------------------------------------*/
	/* compute potentials                                         */
	/*------------------------------------------------------------*/
	/* qp(x     ,z     ) = Fx( ux(x-dx/2,z     ) ) + */
	/*                     Fz( uz(x     ,z-dz/2) )   */
	/* qs(x-dx/2,z-dz/2) = Bz( ux(x-dx/2,z     ) ) + */
	/*                     Bx( uz(x     ,z-dz/2) )   */

	if(opot) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(iz,ix) shared(fdm,uoz,uox,qp,qs,idz,idx)
#endif
	    for    (ix=NOP; ix<fdm->n2pad-NOP; ix++) {
		for(iz=NOP; iz<fdm->n1pad-NOP; iz++) {
		    
		    qp[ix][iz] = Fz( uoz,ix,iz,idz ) 
			+        Fx( uox,ix,iz,idx );
		    
		    qs[ix][iz] = Bz( uox,ix,iz,idz ) 
			-        Bx( uoz,ix,iz,idx );
		}
	    }

	    lint2d_extract(qp,dd[0],cr);
	    lint2d_extract(qs,dd[1],cr);
	    
	    if(snap && it%jsnap==0) {
		cut2d(qp,uc,fdm,acz,acx);
		sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
		
		cut2d(qs,uc,fdm,acz,acx);
		sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	    }
	} else {

	    lint2d_extract(uoz,dd[0],cr);
	    lint2d_extract(uox,dd[1],cr);
	    
	    if(snap && it%jsnap==0) {
		cut2d(uoz,uc,fdm,acz,acx);
		sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

		cut2d(uox,uc,fdm,acz,acx);
		sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	    }
	}
	if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
