/* 2D acoustic time-domain FD modeling  for source perturbation -first order*/
/*
  Copyright (C) 2010 KAUST
  
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
#include "omputil.h"
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
#define DX(a,ix,iz,s) (C2*(a[ix+2][iz  ] - a[ix-2][iz  ]) +  \
                       C1*(a[ix+1][iz  ] - a[ix-1][iz  ])  )*s
#define DZ(a,ix,iz,s) (C2*(a[ix  ][iz+2] - a[ix  ][iz-2]) +  \
                       C1*(a[ix  ][iz+1] - a[ix  ][iz-1])  )*s

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc; 
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fvel=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,az,ax;
    sf_axis as,ar;

    int     nt,nz,nx,ns,nr,nb;
    int     it,iz,ix;
    float   dt,dz,dx,idz,idx;
    float   dl; /*perturbation length*/

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
    float **ro=NULL;           /* density */
    float **roz=NULL;          /* normalized 1st derivative of density on axis 1 */
    float **rox=NULL;          /* normalized 1st derivative of density on axis 2 */
    float **vp=NULL;           /* velocity */
    float **vt=NULL;           /* temporary vp*vp * dt*dt */
    float **dv=NULL;           /* dwdx over w */

    float **um,**uo,**up,**ua,**ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **pm,**po,**pp,**pa,**pt; /* perturbation wavefield: pm = P @ t-1; po = P @ t; pp = P @ t+1 */

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float co,cax,cbx,caz,cbz;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL,**pc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* "exploding reflector" */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    Fden = sf_input ("den"); /* density   */

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
    /* perterbation distance*/
    if(! sf_getfloat("dl",&dl)) dl=0.0; /* dl=0.0  perturbation distance */
    sf_warning("dl=%f",dl);
    /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;        
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
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

	uc=sf_floatalloc2(sf_n(acz),sf_n(acx));
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
    idz = 1/dz;
    idx = 1/dx;

    co = C0 * (idx*idx+idz*idz);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 

    ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    roz =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    rox =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    dv  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

    /* input density */
    sf_floatread(tt[0],nz*nx,Fden);     expand(tt,ro ,fdm);
    /* normalized density derivatives */
    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
	    roz[ix][iz] = DZ(ro,ix,iz,idz) / ro[ix][iz];
	    rox[ix][iz] = DX(ro,ix,iz,idx) / ro[ix][iz];
	}
    }   
    free(*ro); free(ro);

    /* input velocity */
    sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vp,fdm);
    /* precompute vp^2 * dt^2 */
    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    vt[ix][iz] = vp[ix][iz] * vp[ix][iz] * dt*dt;
	    dv[ix][iz] = 0.0;
	}
    }
    /* precompute dwdx over w */
    /*    for(iz=0; iz<fdm->nzpad; iz++){
      dv[0][iz] = 2*(vp[1][iz]-vp[0][iz])/(dx*vp[0][iz]);
      dv[fdm->nxpad-1][iz] = 2*(vp[fdm->nxpad-1][iz]-vp[fdm->nxpad-2][iz])/(dx*vp[fdm->nxpad-1][iz]);
      }*/
    for(ix=NOP; ix<fdm->nxpad-NOP; ix++)
	for(iz=0; iz<fdm->nzpad; iz++) {
	  dv[ix][iz] = 2*DX(vp,ix,iz,idx)/vp[ix][iz];
	};

    if(fsrf) { /* free surface */
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nb; iz++) {
		vt[ix][iz]=0;
		dv[ix][iz] =0;
	    }
	}
    }

    free(*tt); free(tt);    
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    um=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    up=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    pm=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    po=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    pp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    pa=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ua[ix][iz]=0;
	    pm[ix][iz]=0;
	    po[ix][iz]=0;
	    pp[ix][iz]=0;
	    pa[ix][iz]=0;
	}
    }

    /*------------------------------------------------------------*/
    if(dabc) {
	/* one-way abc setup */
	abc = abcone2d_make(NOP,dt,vp,fsrf,fdm);
	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

#ifdef _OPENMP
#pragma omp parallel for				\
    schedule(dynamic,fdm->ompchunk)			\
    private(ix,iz)					\
    shared(fdm,ua,uo,co,cax,caz,cbx,cbz,idx,idz)
#endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
		/* 4th order Laplacian operator */
		ua[ix][iz] = 
		    co * uo[ix  ][iz  ] + 
		    cax*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
		    cbx*(uo[ix-2][iz  ] + uo[ix+2][iz  ]) +
		    caz*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
		    cbz*(uo[ix  ][iz-2] + uo[ix  ][iz+2]);
		
		/* density term */
	/*	ua[ix][iz] -= (
		    DZ(uo,ix,iz,idz) * roz[ix][iz] +
		    DX(uo,ix,iz,idx) * rox[ix][iz] );*/
	    }
	}

	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    lint2d_inject1(ua,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    lint2d_inject(ua,ww,cs);
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iz)		    \
    shared(fdm,ua,uo,um,up,vt)
#endif
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		up[ix][iz] = 2*uo[ix][iz] 
		    -          um[ix][iz] 
		    +          ua[ix][iz] * vt[ix][iz];
	    }
	}
	/* circulate wavefield arrays */
	ut=um;
	um=uo;
	uo=up;
	up=ut;

	if(dabc) {
	    /* one-way abc apply */
	    abcone2d_apply(uo,um,NOP,abc,fdm);
	    sponge2d_apply(um,spo,fdm);
	    sponge2d_apply(uo,spo,fdm);
	    sponge2d_apply(up,spo,fdm);
	}

	/*Perturbation part ------------------------------------------------------------------------------------------------------------*/
	#ifdef _OPENMP
#pragma omp parallel for				\
    schedule(dynamic,fdm->ompchunk)			\
    private(ix,iz)					\
    shared(fdm,pa,po,co,cax,caz,cbx,cbz,idx,idz)
#endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
		/* 4th order Laplacian operator */
		pa[ix][iz] = 
		    co * po[ix  ][iz  ] + 
		    cax*(po[ix-1][iz  ] + po[ix+1][iz  ]) +
		    cbx*(po[ix-2][iz  ] + po[ix+2][iz  ]) +
		    caz*(po[ix  ][iz-1] + po[ix  ][iz+1]) +
		    cbz*(po[ix  ][iz-2] + po[ix  ][iz+2]);
		
		/* density term */
		/*		pa[ix][iz] -= (
		    DZ(po,ix,iz,idz) * roz[ix][iz] +
		    DX(po,ix,iz,idx) * rox[ix][iz] );*/

		/*adding the perturbation source function*/
		pa[ix][iz] -= dv[ix][iz]*ua[ix][iz];
	    }
	}


	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iz)		    \
    shared(fdm,pa,po,pm,pp,vt)
#endif
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		pp[ix][iz] = 2*po[ix][iz] 
		    -          pm[ix][iz] 
		    +          pa[ix][iz] * vt[ix][iz];
		ua[ix][iz] = uo[ix][iz]-pp[ix][iz]*dl;
	    }
	}
	/* circulate wavefield arrays */
	pt=pm;
	pm=po;
	po=pp;
	pp=pt;

	if(dabc) {
	    /* one-way abc apply */
	    abcone2d_apply(po,pm,NOP,abc,fdm);
	    sponge2d_apply(pm,spo,fdm);
	    sponge2d_apply(po,spo,fdm);
	    sponge2d_apply(pp,spo,fdm);
	}
	/*Perturbation part end ------------------------------------------------------------------------------------------------------------*/

	/* extract data */
	lint2d_extract(ua,dd,cr);

	if(snap && it%jsnap==0) {
	    cut2d(ua,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
	if(        it%jdata==0) 
	sf_floatwrite(dd,nr,Fdat);

	/*lint2d_extract(po,dd,cr);

	if(snap && it%jsnap==0) {
	    cut2d(po,pc,fdm,acz,acx);
	    sf_floatwrite(pc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
	if(        it%jdata==0) 
	sf_floatwrite(dd,nr,Fdat);*/

    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
   
    /* deallocate arrays */
    free(*um); free(um);
    free(*up); free(up);
    free(*uo); free(uo);
    free(*ua); free(ua);
    free(*uc); free(uc);

    free(*pm); free(pm);
    free(*pp); free(pp);
    free(*po); free(po);
    free(*pa); free(pa);
    free(*pc); free(pc);

    free(*rox); free(rox);
    free(*roz); free(roz);
    free(*vp);  free(vp);
    free(*vt);  free(vt);
    free(*dv);  free(dv);

    free(ww);
    free(ss);
    free(rr);
    free(dd);
    /*------------------------------------------------------------*/

    if (Fwav!=NULL) sf_fileclose(Fwav);
    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Frec!=NULL) sf_fileclose(Frec);
    if (Fvel!=NULL) sf_fileclose(Fvel);
    if (Fden!=NULL) sf_fileclose(Fden);
    if (Fdat!=NULL) sf_fileclose(Fdat);
    if (Fwfl!=NULL) sf_fileclose(Fwfl);

    /*------------------------------------------------------------*/
    exit (0);
}

