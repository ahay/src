/* 3D TTI time-domain FD modeling */
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

#define Dxx(a,ix,iy,iz,co,ca,cb)	\
    (co* a[iy  ][ix  ][iz  ] +		\
     ca*(a[iy  ][ix-1][iz  ] +	       	\
	 a[iy  ][ix+1][iz  ])+	       	\
     cb*(a[iy  ][ix-2][iz  ] + 		\
	 a[iy  ][ix+2][iz  ]))

#define Dyy(a,ix,iy,iz,co,ca,cb)	\
    (co* a[iy  ][ix  ][iz  ] +		\
     ca*(a[iy-1][ix  ][iz  ] +		\
	 a[iy+1][ix  ][iz  ])+		\
     cb*(a[iy-2][ix  ][iz  ] +		\
	 a[iy+2][ix  ][iz  ]))

#define Dzz(a,ix,iy,iz,co,ca,cb)	\
    (co* a[iy  ][ix  ][iz  ] +		\
     ca*(a[iy  ][ix  ][iz-1] +		\
	 a[iy  ][ix  ][iz+1])+		\
     cb*(a[iy  ][ix  ][iz-2] +		\
	 a[iy  ][ix  ][iz+2]))

#define D1x(a,ix,iy,iz,c1x,c2x)			\
    (c2x*(a[iy  ][ix+2][iz  ] -			\
	  a[iy  ][ix-2][iz  ])+			\
     c1x*(a[iy  ][ix+1][iz  ] -			\
	  a[iy  ][ix-1][iz  ]) )

#define D1y(a,ix,iy,iz,c1y,c2y)			\
    (c2y*(a[iy+2][ix  ][iz  ] -			\
	  a[iy-2][ix  ][iz  ])+			\
     c1y*(a[iy+1][ix  ][iz  ] -			\
	  a[iy-1][ix  ][iz  ]) )

#define D1z(a,ix,iy,iz,c1z,c2z)			\
    (c2z*(a[iy  ][ix  ][iz+2] -			\
	  a[iy  ][ix  ][iz-2])+			\
     c1z*(a[iy  ][ix  ][iz+1] -			\
	  a[iy  ][ix  ][iz-1]) )

#define Dxy(a,ix,iy,iz,c1x,c2x,c1y,c2y)		\
    (c2y *(D1x(a,ix,iy+2,iz,c1x,c2x) -		\
	   D1x(a,ix,iy-2,iz,c1x,c2x))+		\
     c1y *(D1x(a,ix,iy+1,iz,c1x,c2x) -		\
	   D1x(a,ix,iy-1,iz,c1x,c2x)))

#define Dyz(a,ix,iy,iz,c1y,c2y,c1z,c2z)		\
    (c2z *(D1y(a,ix,iy,iz+2,c1y,c2y) -		\
	   D1y(a,ix,iy,iz-2,c1y,c2y))+		\
     c1z *(D1y(a,ix,iy,iz+1,c1y,c2y) -		\
	   D1y(a,ix,iy,iz-1,c1y,c2y)))

#define Dzx(a,ix,iy,iz,c1z,c2z,c1x,c2x)		\
    (c2x *(D1z(a,ix+2,iy,iz,c1z,c2z) -		\
	   D1z(a,ix-2,iy,iz,c1z,c2z))+		\
     c1x *(D1z(a,ix+1,iy,iz,c1z,c2z) -		\
	   D1z(a,ix-1,iy,iz,c1z,c2z)))

#define H1(a,ix,iy,iz,						\
	   st,ct,sp,cp,						\
	   cox,cax,cbx,c1x,c2x,					\
	   coy,cay,cby,c1y,c2y,					\
	   coz,caz,cbz,c1z,c2z)					\
    (st*st   *  cp*cp * Dxx(a,ix,iy,iz,cox,cax,cbx) +		\
     st*st   *  sp*sp * Dyy(a,ix,iy,iz,coy,cay,cby) +		\
     ct*ct            * Dzz(a,ix,iy,iz,coz,caz,cbz) +		\
     st*st   *2*sp*cp * Dxy(a,ix,iy,iz,c1x,c2x,c1y,c2y) +	\
     2*st*ct *  sp    * Dyz(a,ix,iy,iz,c1y,c2y,c1z,c2z) +	\
     2*st*ct *  cp    * Dzx(a,ix,iy,iz,c1z,c2z,c1x,c2x))

#define H2(a,ix,iy,iz,						\
	   st,ct,sp,cp,						\
	   cox,cax,cbx,c1x,c2x,					\
	   coy,cay,cby,c1y,c2y,					\
	   coz,caz,cbz,c1z,c2z)					\
    ((1-st*st*  cp*cp)* Dxx(a,ix,iy,iz,cox,cax,cbx) +		\
     (1-st*st*  sp*sp)* Dyy(a,ix,iy,iz,coy,cay,cby) +		\
     st*st            * Dzz(a,ix,iy,iz,coz,caz,cbz) -		\
     st*st   *2*sp*cp * Dxy(a,ix,iy,iz,c1x,c2x,c1y,c2y) -	\
     2*st*ct *  sp    * Dyz(a,ix,iy,iz,c1y,c2y,c1z,c2z) -	\
     2*st*ct *  cp    * Dzx(a,ix,iy,iz,c1z,c2z,c1x,c2x))

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
    sf_axis at,az,ax,ay;
    sf_axis as,ar;

    int     nt,nz,nx,ny,ns,nr,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy,dt2;

    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    float  *dd=NULL;           /* data      */

    float ***tt=NULL;
    float ***vp=NULL;           /* velocity */

    float ***vpn=NULL;
    float ***vpz=NULL;
    float ***vpx=NULL;
    float ***vsz=NULL;

    float ***tht=NULL,***sit=NULL,***cot=NULL;
    float ***phi=NULL,***sip=NULL,***cop=NULL;
    float st,ct,sp,cp;

    float ***pm=NULL,***po=NULL,***pp=NULL,***pa=NULL,***pt=NULL; /*      main wavefield */
    float ***qm=NULL,***qo=NULL,***qp=NULL,***qa=NULL,***qt=NULL; /* auxiliary wavefield */
    float ***sf=NULL; /* "stress" field */

    /* linear inteppolation weights/indices */
    lint3d cs,cr;

    /* FD operator size */
    float cox,cax,cbx,c1x,c2x;
    float coy,cay,cby,c1y,c2y;
    float coz,caz,cbz,c1z,c2z;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***pc=NULL;

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
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* t */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* y */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

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
	if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
	if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);

	dqz=sf_d(az);
	dqx=sf_d(ax);
	dqy=sf_d(ay);

	acz = sf_maxa(nqz,oqz,dqz);
	acx = sf_maxa(nqx,oqx,dqx);
	acy = sf_maxa(nqy,oqy,dqy);
	/* check if the imaging window fits in the wavefield domain */

	pc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,acy,3);
	sf_oaxa(Fwfl,at, 4);
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad);
    /*------------------------------------------------------------*/

    if(expl) {
	ww = sf_floatalloc( 1);
    } else {
	ww = sf_floatalloc(ns);
    }
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    cox = C0 / (dx*dx);
    cax = CA / (dx*dx);
    cbx = CB / (dx*dx);
    c1x = C1 / dx;
    c2x = C2 / dx;

    coy = C0 / (dy*dy);
    cay = CA / (dy*dy);
    cby = CB / (dy*dy);
    c1y = C1 / dy;
    c2y = C2 / dy;

    coz = C0 / (dz*dz);
    caz = CA / (dz*dz);
    cbz = CB / (dz*dz);
    c1z = C1 / dz;
    c2z = C2 / dz;

    /* precompute dt^2*/
    dt2 = dt*dt;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 
    /*------------------------------------------------------------*/

    /* input velocity */
    vp  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

    vpz =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    sf_floatread(tt[0][0],nz*nx*ny,Fvel ); 
    expand3d(tt,vpz,fdm); /* VPz */

    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {	    
		vp [iy][ix][iz] = vpz[iy][ix][iz];
		vpz[iy][ix][iz] = vpz[iy][ix][iz] * vpz[iy][ix][iz];
	    }
	}
    }
    if(fsrf) { /* free surface */
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb;    iz++) {
		    vpz[iy][ix][iz]=0;
		}
	    }
	}
    }

    if(atype[0] != 'i') {
	vpn =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);     
	sf_floatread(tt[0][0],nz*nx*ny,Fvel );    
	expand3d(tt,vpn,fdm); /* VPn */

	vpx =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	sf_floatread(tt[0][0],nz*nx*ny,Fvel );    
	expand3d(tt,vpx,fdm); /* VPx */

	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {	    
		    vpn[iy][ix][iz] = vpn[iy][ix][iz] * vpn[iy][ix][iz];
		    vpx[iy][ix][iz] = vpx[iy][ix][iz] * vpx[iy][ix][iz];
		}
	    }
	}

	if(fsrf) { /* free surface */
	    for        (iy=0; iy<fdm->nypad; iy++) {
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nb;    iz++) {
			vpn[iy][ix][iz]=0;
			vpx[iy][ix][iz]=0;
		    }
		}
	    }
	}

	if(uses) {
	    vsz =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	    sf_floatread(tt[0][0],nz*nx*ny,Fvel );    
	    expand3d(tt,vsz,fdm); /* VSz */
	    for        (iy=0; iy<fdm->nypad; iy++) {
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nzpad; iz++) {
			vsz[iy][ix][iz] = vsz[iy][ix][iz] * vsz[iy][ix][iz];
		    }
		}
	    }
	}
    }

    /*------------------------------------------------------------*/
    if( atype[0]=='t') {
	/* input tilt angles */

	tht =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	sf_floatread(tt[0][0],nz*nx*ny,Fang); 
	expand3d(tt,tht,fdm);

	phi =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	sf_floatread(tt[0][0],nz*nx*ny,Fang); 
	expand3d(tt,phi,fdm);

	sit =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	cot =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

	sip =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	cop =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {	    
		    tht[iy][ix][iz] *= SF_PI/180.;
		    sit[iy][ix][iz] =   sinf(tht[iy][ix][iz]);
		    cot[iy][ix][iz] =   cosf(tht[iy][ix][iz]);

		    phi[iy][ix][iz] *= SF_PI/180.;
		    sip[iy][ix][iz] =   sinf(phi[iy][ix][iz]);
		    cop[iy][ix][iz] =   cosf(phi[iy][ix][iz]);
		}
	    }
	}

	free(**tht); free(*tht); free(tht);
	free(**phi); free(*phi); free(phi);
    }

    /*------------------------------------------------------------*/
    free(**tt); free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    pm=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    po=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    pp=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    pa=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		pm[iy][ix][iz]=0;
		po[iy][ix][iz]=0;
		pp[iy][ix][iz]=0;
		pa[iy][ix][iz]=0;
	    }
	}
    }
    
    if(atype[0] != 'i') {
	qm=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qp=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qa=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    qm[iy][ix][iz]=0;
		    qo[iy][ix][iz]=0;
		    qp[iy][ix][iz]=0;
		    qa[iy][ix][iz]=0;
		}
	    }
	}

	if(sout) sf=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    }

    /*------------------------------------------------------------*/
    if(dabc) {
	abc = abcone3d_make(NOP,dt,vp,fsrf,fdm); /* one-way */
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
    private(ix,iy,iz,H1p,H2p,H1q,H2q,st,ct,sp,cp)			\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coy,cay,cby,c1y,c2y,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpz,vpx,vsz,						\
	   sit,cot,sip,cop)
#endif
		    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {		    
			for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
				
				st=sit[iy][ix][iz];
				ct=cot[iy][ix][iz];
				sp=sip[iy][ix][iz];
				cp=cop[iy][ix][iz];
				
				H1p = H1(po,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);
				
				H2p = H2(po,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);

				H1q = H1(qo,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);
				
				H2q = H2(qo,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);
				
				/* p - main field */
				pa[iy][ix][iz] = 
				    H1p * vsz[iy][ix][iz] +
				    H2p * vpx[iy][ix][iz] + 
				    H1q * vpz[iy][ix][iz] -
				    H1q * vsz[iy][ix][iz];
				
				/* q - auxiliary field */
				qa[iy][ix][iz] = 
				    H2p * vpn[iy][ix][iz] -
				    H2p * vsz[iy][ix][iz] +
				    H1q * vpz[iy][ix][iz] +
				    H2q * vsz[iy][ix][iz];
				
			    }
			}
		    }

		} else {

#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,H2p,H1q,st,ct,sp,cp)				\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coy,cay,cby,c1y,c2y,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpz,vpx,							\
	   sit,cot,sip,cop)
#endif
		    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {	
			for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
				
				st=sit[iy][ix][iz];
				ct=cot[iy][ix][iz];
				sp=sip[iy][ix][iz];
				cp=cop[iy][ix][iz];
				
				H2p = H2(po,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);

				H1q = H1(qo,ix,iy,iz,			\
					 st,ct,sp,cp,			\
					 cox,cax,cbx,c1x,c2x,		\
					 coy,cay,cby,c1y,c2y,		\
					 coz,caz,cbz,c1z,c2z);
				
				/* p - main field */
				pa[iy][ix][iz] = 
				    H2p * vpx[iy][ix][iz] + 
				    H1q * vpz[iy][ix][iz];
				
				/* q - auxiliary field */
				qa[iy][ix][iz] = 
				    H2p * vpn[iy][ix][iz] +
				    H1q * vpz[iy][ix][iz];
			    }
			}
		    }

		}
		break;
		    
	    case 'v':

		if(uses) {

#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,H1p,H2p,H1q,H2q)					\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coy,cay,cby,c1y,c2y,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpz,vpx,vsz)
#endif
		    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {	
			for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
				H1p = Dzz(po,ix,iy,iz,coz,caz,cbz);
				H1q = Dzz(qo,ix,iy,iz,coz,caz,cbz);
				
				H2p = Dxx(po,ix,iy,iz,cox,cax,cbx) \
				    + Dyy(po,ix,iy,iz,coy,cay,cby);
				H2q = Dxx(qo,ix,iy,iz,cox,cax,cbx) \
				    + Dyy(qo,ix,iy,iz,coy,cay,cby);
				
				/* p - main field */
				pa[iy][ix][iz] = 
				    H1p * vsz[iy][ix][iz] +
				    H2p * vpx[iy][ix][iz] + 
				    H1q * vpz[iy][ix][iz] -
				    H1q * vsz[iy][ix][iz];
				
				/* q - auxiliary field */
				qa[iy][ix][iz] = 
				    H2p * vpn[iy][ix][iz] -
				    H2p * vsz[iy][ix][iz] +
				    H1q * vpz[iy][ix][iz] +
				    H2q * vsz[iy][ix][iz];
			    }
			} 
		    }

		} else {

#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,H2p,H1q)						\
    shared(fdm,pa,po,qa,qo,						\
	   cox,cax,cbx,c1x,c2x,						\
	   coy,cay,cby,c1y,c2y,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpn,vpx,vpz)
#endif
		    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {	
			for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
				H1q = Dzz(qo,ix,iy,iz,coz,caz,cbz);
			    
				H2p = Dxx(po,ix,iy,iz,cox,cax,cbx) \
				    + Dyy(po,ix,iy,iz,coy,cay,cby);
				
				/* p - main field */
				pa[iy][ix][iz] = 
				    H2p * vpx[iy][ix][iz] + 
				    H1q * vpz[iy][ix][iz];
				
				/* q - auxiliary field */
				qa[iy][ix][iz] = 
				    H2p * vpn[iy][ix][iz] +
				    H1q * vpz[iy][ix][iz];
			    }
			} 
		    }
		}
		break;
		
	    case 'i':
	    default:
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,pa,po,							\
	   cox,cax,cbx,c1x,c2x,						\
	   coy,cay,cby,c1y,c2y,						\
	   coz,caz,cbz,c1z,c2z,						\
	   vpz)
#endif
		for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    pa[iy][ix][iz] = ( Dxx(po,ix,iy,iz,cox,cax,cbx) +
					       Dyy(po,ix,iy,iz,coy,cay,cby) +
					       Dzz(po,ix,iy,iz,coz,caz,cbz) ) * vpz[iy][ix][iz];
			    
			}
		    }   
		}
		break;
	}

	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    ;                   lint3d_inject1(pa,ww[0],cs);
	    if(atype[0] != 'i') lint3d_inject1(qa,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    ;                   lint3d_inject(pa,ww,cs);
	    if(atype[0] != 'i') lint3d_inject(qa,ww,cs);
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iy,iz)		    \
    shared(fdm,pa,po,pm,pp,dt2)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    pp[iy][ix][iz] = 2*po[iy][ix][iz] 
			-              pm[iy][ix][iz] 
			+              pa[iy][ix][iz] * dt2;
		}
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
    private(ix,iy,iz)				\
    shared(fdm,qa,qo,qm,qp,dt2)
#endif
	    for        (iy=0; iy<fdm->nypad; iy++) {
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nzpad; iz++) {
			qp[iy][ix][iz] = 2*qo[iy][ix][iz] 
			    -              qm[iy][ix][iz] 
			    +              qa[iy][ix][iz] * dt2;
		    }
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
	    abcone3d_apply(po,pm,NOP,abc,fdm);
	    sponge3d_apply(pm,spo,fdm);
	    sponge3d_apply(po,spo,fdm);
	    
	    if(atype[0] != 'i') {
		abcone3d_apply(qo,qm,NOP,abc,fdm);
		sponge3d_apply(qm,spo,fdm);
		sponge3d_apply(qo,spo,fdm);
	    }
	}
	
	/* compute stress */
	if(sout && (atype[0] != 'i')) {
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iy,iz)				\
    shared(fdm,po,qo,sf)
#endif
	    for        (iy=0; iy<fdm->nypad; iy++) {	    
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nzpad; iz++) {
			sf[iy][ix][iz] = 2*po[iy][ix][iz] + qo[iy][ix][iz];
		    }
		}
	    }
	}

	/* extract data at receivers */
	if(sout && (atype[0] != 'i')) {lint3d_extract(sf,dd,cr);
	} else {                       lint3d_extract(po,dd,cr);}
	if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

	/* extract wavefield in the "box" */
	if(snap && it%jsnap==0) {
	    if(sout && (atype[0] != 'i')) {cut3d(sf,pc,fdm,acz,acx,acy);
	    } else {                       cut3d(po,pc,fdm,acz,acx,acy);}
	    sf_floatwrite(pc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
	}

    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */

    free(**pm); free(*pm); free(pm);
    free(**pp); free(*pp); free(pp);
    free(**po); free(*po); free(po);
    free(**pa); free(*pa); free(pa);
    free(**pc); free(*pc); free(pc);

    free(**vp);  free(*vp);  free(vp);
    free(**vpz); free(*vpz); free(vpz);

    if(atype[0] != 'i') {
	free(**qm); free(*qm); free(qm);
	free(**qp); free(*qp); free(qp);
	free(**qo); free(*qo); free(qo);
	free(**qa); free(*qa); free(qa);

	free(**vpn); free(*vpn); free(vpn);
	free(**vpx); free(*vpx); free(vpx);

	if(uses){ free(**vsz); free(*vsz); free(vsz); }
	if(sout){ free(**sf);  free(*sf);  free(sf);  }
    }

    if(atype[0] == 't') {
	free(**sit); free(*sit); free(sit);
	free(**cot); free(*cot); free(cot);
	free(**sip); free(*sip); free(sip);
	free(**cop); free(*cop); free(cop);
    }

    free(ww);
    free(ss);
    free(rr);
    free(dd);
    /*------------------------------------------------------------*/

    exit (0);
}

