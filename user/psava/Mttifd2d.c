/* 2D TTI time-domain FD modeling */
/*
  Copyright (C) 2011 Colorado School of Mines
  
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

#include <math.h>
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

#define Dxx(a,ix,iz,co,ca,cb)	\
    (co* a[ix  ][iz  ] +	\
     ca*(a[ix-1][iz  ] +	\
	 a[ix+1][iz  ])+	\
     cb*(a[ix-2][iz  ] +	\
	 a[ix+2][iz  ]))

#define Dzz(a,ix,iz,co,ca,cb)	\
    (co* a[ix  ][iz  ] +	\
     ca*(a[ix  ][iz-1] +	\
	 a[ix  ][iz+1])+	\
     cb*(a[ix  ][iz-2] +	\
	 a[ix  ][iz+2]))

#define D1x(a,ix,iz,c1x,c2x) \
    (c2x*(a[ix+2][iz  ] -    \
	  a[ix-2][iz  ])+    \
     c1x*(a[ix+1][iz  ] -    \
	  a[ix-1][iz  ]) )

#define D1z(a,ix,iz,c1z,c2z) \
    (c2z*(a[ix  ][iz+2] -    \
	  a[ix  ][iz-2])+    \
     c1z*(a[ix  ][iz+1] -    \
	  a[ix  ][iz-1]) )

#define Dzx(a,ix,iz,c1z,c2z,c1x,c2x) \
    (c2x *(D1z(a,ix+2,iz,c1z,c2z) -  \
	   D1z(a,ix-2,iz,c1z,c2z))+  \
     c1x *(D1z(a,ix+1,iz,c1z,c2z) -  \
	   D1z(a,ix-1,iz,c1z,c2z)))

#define H1(a,ix,iz,				\
	   st,ct,				\
	   cox,cax,cbx,c1x,c2x,			\
	   coz,caz,cbz,c1z,c2z)			\
    (st*st   * Dxx(a,ix,iz,cox,cax,cbx) +	\
     ct*ct   * Dzz(a,ix,iz,coz,caz,cbz) +	\
     2*st*ct * Dzx(a,ix,iz,c1z,c2z,c1x,c2x))

#define H2(a,ix,iz,				\
	   st,ct,				\
	   cox,cax,cbx,c1x,c2x,			\
	   coz,caz,cbz,c1z,c2z)			\
    (ct*ct   * Dxx(a,ix,iz,cox,cax,cbx) +	\
     st*st   * Dzz(a,ix,iz,coz,caz,cbz) -	\
     2*st*ct * Dzx(a,ix,iz,c1z,c2z,c1x,c2x))

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc,sout,uses;
    int  jsnap,ntsnap,jdata;
    char *atype;
#ifdef _OPENMP
    int ompnth=1;
#endif

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fvel=NULL; /* velocity  */
    sf_file Fang=NULL; /* angles    */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,az,ax;
    sf_axis as,ar;

    int     nt,nz,nx,ns,nr,nb;
    int     it,iz,ix;
    float   dt,dz,dx,dt2;

    /* FDM structure */
    fdm2d    fdm=NULL;
    abcone2d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
    float  *dd=NULL;           /* data      */

    float **tt=NULL;
    float **vp=NULL;           /* velocity */

    float **vpn=NULL;
    float **vpz=NULL;
    float **vpx=NULL;
    float **vsz=NULL;

    float **tht=NULL,**sit=NULL,**cot=NULL;
    float st,ct;

    float **pm=NULL,**po=NULL,**pp=NULL,**pa=NULL,**pt=NULL; /*      main wavefield */
    float **qm=NULL,**qo=NULL,**qp=NULL,**qa=NULL,**qt=NULL; /* auxiliary wavefield */
    float **sf=NULL; /* "stress" field */

    /* linear inteppolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float cox,cax,cbx,c1x,c2x;
    float coz,caz,cbz,c1z,c2z;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **pc=NULL;

    float H1p,H2p,H1q,H2q;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* select anisotropy model */
    if (NULL == (atype = sf_getstring("atype"))) atype = "i";
    switch(atype[0]) {
	case 't':
	    sf_warning("TTI model");
	    break;

	case 'v':
	    sf_warning("VTI model");
	    break;

	case 'i':
	default:
	    sf_warning("ISO model");
	    break;
    }

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
	if(!ompnth)
		abort();
#endif
	/*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface */
    if(! sf_getbool("expl",&expl)) expl=false; /* "exploding reflector" */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    if(! sf_getbool("sout",&sout)) sout=false; /* stress output */
    if(! sf_getbool("uses",&uses)) uses=false; /* use vsz */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fang = sf_input ("ang"); /* angles    */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

    ns = sf_n(as);
    nr = sf_n(ar);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;        
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);

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

	pc=sf_floatalloc2(sf_n(acz),sf_n(acx));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,at, 3);
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/

    if(expl) {
	ww = sf_floatalloc( 1);
    } else {
	ww = sf_floatalloc(ns);
    }
    dd = sf_floatalloc(nr);

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
    cox = C0 / (dx*dx);
    cax = CA / (dx*dx);
    cbx = CB / (dx*dx);
    c1x = C1 / dx;
    c2x = C2 / dx;

    coz = C0 / (dz*dz);
    caz = CA / (dz*dz);
    cbz = CB / (dz*dz);
    c1z = C1 / dz;
    c2z = C2 / dz;

    /* precompute dt^2*/
    dt2 = dt*dt;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 
    /*------------------------------------------------------------*/

    /* input velocity */
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

    vpz =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(tt[0],nz*nx,Fvel ); 
    expand(tt,vpz,fdm); /* VPz */

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {	    
	    vp [ix][iz] = vpz[ix][iz];
	    vpz[ix][iz] = vpz[ix][iz] * vpz[ix][iz];
	}
    }
    if(fsrf) { /* free surface */
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nb; iz++) {
		vpz[ix][iz]=0;
	    }
	}
    }

    if(atype[0] != 'i') {
	vpn =sf_floatalloc2(fdm->nzpad,fdm->nxpad);     
	sf_floatread(tt[0],nz*nx,Fvel );    
	expand(tt,vpn,fdm); /* VPn */

	vpx =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	sf_floatread(tt[0],nz*nx,Fvel );    
	expand(tt,vpx,fdm); /* VPx */

	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {	    
		vpn[ix][iz] = vpn[ix][iz] * vpn[ix][iz];
		vpx[ix][iz] = vpx[ix][iz] * vpx[ix][iz];
	    }
	}

	if(fsrf) { /* free surface */
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb; iz++) {
		    vpn[ix][iz]=0;
		    vpx[ix][iz]=0;
		}
	    }
	}

	if(uses) {
	    vsz =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	    sf_floatread(tt[0],nz*nx,Fvel );    
	    expand(tt,vsz,fdm); /* VSz */
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    vsz[ix][iz] = vsz[ix][iz] * vsz[ix][iz];
		}
	    }
	}
    }

    /*------------------------------------------------------------*/
    if( atype[0]=='t') {
	/* input tilt angle */
	tht =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

	sit =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	cot =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

	sf_floatread(tt[0],nz*nx,Fang); 
	expand(tt,tht,fdm);
	
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {	    
		tht[ix][iz] *= SF_PI/180.;
		sit[ix][iz] =   sinf(tht[ix][iz]);
		cot[ix][iz] =   cosf(tht[ix][iz]);
	    }
	}

	free(*tht); free(tht);
    }

    /*------------------------------------------------------------*/
    free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    pm=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    po=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    pp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    pa=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    pm[ix][iz]=0;
	    po[ix][iz]=0;
	    pp[ix][iz]=0;
	    pa[ix][iz]=0;
	}
    }
    
    if(atype[0] != 'i') {
	qm=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	qo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	qp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	qa=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		qm[ix][iz]=0;
		qo[ix][iz]=0;
		qp[ix][iz]=0;
		qa[ix][iz]=0;
	    }
	}

	if(sout) sf=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    }

    /*------------------------------------------------------------*/
    if(dabc) {
	abc = abcone2d_make(NOP,dt,vp,fsrf,fdm); /* one-way */
	spo = sponge_make(fdm->nb);              /* sponge  */
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

	/* compute acceleration */
	switch(atype[0]) {
	    case 't':

		if(uses) {
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iz,H1p,H2p,H1q,H2q,st,ct)				\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpz,vpx,vsz,						\
	   sit,cot)
#endif
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    st=sit[ix][iz];
			    ct=cot[ix][iz];
			    
			    H1p = H1(po,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    H2p = H2(po,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    H1q = H1(qo,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    H2q = H2(qo,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    /* p - main field */
			    pa[ix][iz] = 
				H1p * vsz[ix][iz] +
				H2p * vpx[ix][iz] + 
				H1q * vpz[ix][iz] -
				H1q * vsz[ix][iz];
			    
			    /* q - auxiliary field */
			    qa[ix][iz] = 
				H2p * vpn[ix][iz] -
				H2p * vsz[ix][iz] +
				H1q * vpz[ix][iz] +
				H2q * vsz[ix][iz];
			    
			}
		    }
		} else {
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iz,H2p,H1q,st,ct)					\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpz,vpx,							\
	   sit,cot)
#endif
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    st=sit[ix][iz];
			    ct=cot[ix][iz];
			    
			    H2p = H2(po,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    H1q = H1(qo,ix,iz,				\
				     st,ct,				\
				     cox,cax,cbx,c1x,c2x,		\
				     coz,caz,cbz,c1z,c2z);
			    
			    /* p - main field */
			    pa[ix][iz] = 
				H2p * vpx[ix][iz] + 
				H1q * vpz[ix][iz];
			    
			    /* q - auxiliary field */
			    qa[ix][iz] = 
				H2p * vpn[ix][iz] +
				H1q * vpz[ix][iz];
			}
		    }

		}
		break;
		    
	    case 'v':

		if(uses) {
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iz,H1p,H2p,H1q,H2q)					\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,							\
	   coz,caz,cbz,							\
	   vpn,vpz,vpx,vsz)
#endif
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    H1p = Dzz(po,ix,iz,coz,caz,cbz);
			    H1q = Dzz(qo,ix,iz,coz,caz,cbz);
			    
			    H2p = Dxx(po,ix,iz,cox,cax,cbx);
			    H2q = Dxx(qo,ix,iz,cox,cax,cbx);
			    
			    /* p - main field */
			    pa[ix][iz] = 
				H1p * vsz[ix][iz] +
				H2p * vpx[ix][iz] + 
				H1q * vpz[ix][iz] -
				H1q * vsz[ix][iz];
			    
			    /* q - auxiliary field */
			    qa[ix][iz] = 
				H2p * vpn[ix][iz] -
				H2p * vsz[ix][iz] +
				H1q * vpz[ix][iz] +
				H2q * vsz[ix][iz];
			}
		    } 
		} else {
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iz,H2p,H1q)					\
    shared(fdm,pa,po,qa,qo,					\
	   cox,cax,cbx,						\
	   coz,caz,cbz,						\
	   vpn,vpx,vpz)
#endif
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    H1q = Dzz(qo,ix,iz,coz,caz,cbz);			    
			    H2p = Dxx(po,ix,iz,cox,cax,cbx);
			    
			    /* p - main field */
			    pa[ix][iz] = 
				H2p * vpx[ix][iz] + 
				H1q * vpz[ix][iz];
			    
			    /* q - auxiliary field */
			    qa[ix][iz] = 
				H2p * vpn[ix][iz] +
				H1q * vpz[ix][iz];
			}
		    } 
		}
		break;
		
	    case 'i':
	    default:
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iz)						\
    shared(fdm,pa,po,						\
	   cox,cax,cbx,						\
	   coz,caz,cbz,						\
	   vpz)
#endif
		for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			
			pa[ix][iz] = ( Dxx(po,ix,iz,cox,cax,cbx) + 
				       Dzz(po,ix,iz,coz,caz,cbz) ) * vpz[ix][iz];
			
		    }
		}   
		break;
	}

	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    ;                   lint2d_inject1(pa,ww[0],cs);
	    if(atype[0] != 'i') lint2d_inject1(qa,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    ;                   lint2d_inject(pa,ww,cs);
	    if(atype[0] != 'i') lint2d_inject(qa,ww,cs);
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iz)		    \
    shared(fdm,pa,po,pm,pp,dt2)
#endif
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		pp[ix][iz] = 2*po[ix][iz] 
		    -          pm[ix][iz] 
		    +          pa[ix][iz] * dt2;
	    }
	}
	/* circulate wavefield arrays */
	pt=pm;
	pm=po;
	po=pp;
	pp=pt;
	
	if(atype[0] != 'i') {
	    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iz)				\
    shared(fdm,qa,qo,qm,qp,dt2)
#endif
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    qp[ix][iz] = 2*qo[ix][iz] 
			-          qm[ix][iz] 
			+          qa[ix][iz] * dt2;
		}
	    }
	    /* circulate wavefield arrays */
	    qt=qm;
	    qm=qo;
	    qo=qp;
	    qp=qt;
	}

	/* one-way abc apply */
	if(dabc) {
	    abcone2d_apply(po,pm,NOP,abc,fdm);
	    sponge2d_apply(pm,spo,fdm);
	    sponge2d_apply(po,spo,fdm);
	    
	    if(atype[0] != 'i') {
		abcone2d_apply(qo,qm,NOP,abc,fdm);
		sponge2d_apply(qm,spo,fdm);
		sponge2d_apply(qo,spo,fdm);
	    }
	}
	
	/* compute stress */
	if(sout && (atype[0] != 'i')) {
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iz)				\
    shared(fdm,po,qo,sf)
#endif
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    sf[ix][iz] = po[ix][iz] + qo[ix][iz];
		}
	    }
	}

	/* extract data at receivers */
	if(sout && (atype[0] != 'i')) {lint2d_extract(sf,dd,cr);
	} else {                       lint2d_extract(po,dd,cr);}
	if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

	/* extract wavefield in the "box" */
	if(snap && it%jsnap==0) {
	    if(sout && (atype[0] != 'i')) {cut2d(sf,pc,fdm,acz,acx);
	    } else {                       cut2d(po,pc,fdm,acz,acx);}
	    sf_floatwrite(pc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}

    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */

    free(*pm); free(pm);
    free(*pp); free(pp);
    free(*po); free(po);
    free(*pa); free(pa);
    free(*pc); free(pc);

    free(*vp);  free(vp);
    free(*vpz); free(vpz);

    if(atype[0] != 'i') {
	free(*qm); free(qm);
	free(*qp); free(qp);
	free(*qo); free(qo);
	free(*qa); free(qa);

	free(*vpn); free(vpn);
	free(*vpx); free(vpx);

	if(uses){ free(*vsz); free(vsz); }
	if(sout){ free(*sf);  free(sf);  }
    }

    if(atype[0] == 't') {
	free(*sit); free(sit);
	free(*cot); free(cot);
    }

    free(ww);
    free(ss);
    free(rr);
    free(dd);
    /*------------------------------------------------------------*/

    exit (0);
}

