/* 3D acoustic time-domain FD modeling */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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
#define DX(a,ix,iy,iz,s) (C2*(a[iy  ][ix+2][iz  ] - a[iy  ][ix-2][iz  ]) +  \
                          C1*(a[iy  ][ix+1][iz  ] - a[iy  ][ix-1][iz  ])  )*s

#define DY(a,ix,iy,iz,s) (C2*(a[iy+2][ix  ][iz  ] - a[iy-2][ix  ][iz  ]) +  \
                          C1*(a[iy+1][ix  ][iz  ] - a[iy-1][ix  ][iz  ])  )*s

#define DZ(a,ix,iy,iz,s) (C2*(a[iy  ][ix  ][iz+2] - a[iy  ][ix  ][iz-2]) +  \
                          C1*(a[iy  ][ix  ][iz+1] - a[iy  ][ix  ][iz-1])  )*s

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc; 
    int  jsnap,ntsnap,jdata;

    /* OMP parameters */
#ifdef _OPENMP
    int ompnth;
#endif 

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fvel=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,az,ax,ay;
    sf_axis as,ar;

    int     nt,nz,nx,ny,ns,nr,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy,idz,idx,idy;

    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    float  *dd=NULL;           /* data      */

    float***tt=NULL;
    float***ro=NULL;           /* density */
    float***roz=NULL;          /* normalized 1st derivative of density on axis 1 */
    float***rox=NULL;          /* normalized 1st derivative of density on axis 2 */
    float***roy=NULL;          /* normalized 1st derivative of density on axis 3 */
    float***vp=NULL;           /* velocity */
    float***vt=NULL;           /* temporary vp*vp * dt*dt */

    float***um,***uo,***up,***ua,***ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint3d cs,cr;

    /* FD operator size */
    float co,cax,cbx,cay,cby,caz,cbz;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***uc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
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
    ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space */

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
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);
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

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); sf_raxa(acy);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

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
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;

    co = C0 * (idx*idx+idy*idy+idz*idz);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;
    cay= CA *  idy*idy;
    cby= CB *  idy*idy;
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 

    ro  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    roz =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    rox =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    roy =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vp  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    vt  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

    /* input density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden);     expand3d(tt,ro ,fdm);

    /* normalized density derivatives */
    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		roz[iy][ix][iz] = DZ(ro,ix,iy,iz,idz) / ro[iy][ix][iz];
		rox[iy][ix][iz] = DX(ro,ix,iy,iz,idx) / ro[iy][ix][iz];
		roy[iy][ix][iz] = DY(ro,ix,iy,iz,idy) / ro[iy][ix][iz];
	    }
	}   
    }
    free(**ro);  free(*ro); free(ro);  

    /* input velocity */
    sf_floatread(tt[0][0],nz*nx*ny,Fvel );    expand3d(tt,vp,fdm);
    /* precompute vp^2 * dt^2 */
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		vt[iy][ix][iz] = vp[iy][ix][iz] * vp[iy][ix][iz] * dt*dt;
	    }
	}
    }
    if(fsrf) { /* free surface */
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb; iz++) {
		    vt[iy][ix][iz]=0;
		}
	    }
	}
    }

    free(**tt);  free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    um=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    up=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ua=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		um[iy][ix][iz]=0;
		uo[iy][ix][iz]=0;
		up[iy][ix][iz]=0;
		ua[iy][ix][iz]=0;
	    }
	}
    }

    /*------------------------------------------------------------*/
    if(dabc) {
	/* one-way abc setup */
	abc = abcone3d_make(NOP,dt,vp,fsrf,fdm);
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
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,ua,uo,co,cax,cay,caz,cbx,cby,cbz,idx,idy,idz)
#endif
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		    
		    /* 4th order Laplacian operator */
		    ua[iy][ix][iz] = 
			co * uo[iy  ][ix  ][iz  ] + 
			cax*(uo[iy  ][ix-1][iz  ] + uo[iy  ][ix+1][iz  ]) +
			cbx*(uo[iy  ][ix-2][iz  ] + uo[iy  ][ix+2][iz  ]) +
			cay*(uo[iy-1][ix  ][iz  ] + uo[iy+1][ix  ][iz  ]) +
			cby*(uo[iy-2][ix  ][iz  ] + uo[iy+2][ix  ][iz  ]) +
			caz*(uo[iy  ][ix  ][iz-1] + uo[iy  ][ix  ][iz+1]) +
			cbz*(uo[iy  ][ix  ][iz-2] + uo[iy  ][ix  ][iz+2]);
		    
		    /* density term */
		    ua[iy][ix][iz] -= (
			DZ(uo,ix,iy,iz,idz) * roz[iy][ix][iz] +
			DX(uo,ix,iy,iz,idx) * rox[iy][ix][iz] +
			DY(uo,ix,iy,iz,idy) * roy[iy][ix][iz] );
		}
	    }   
	}

	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    lint3d_inject1(ua,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    lint3d_inject(ua,ww,cs);
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iy,iz)		    \
    shared(fdm,ua,uo,um,up,vt)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    up[iy][ix][iz] = 2*uo[iy][ix][iz] 
			-              um[iy][ix][iz] 
			+              ua[iy][ix][iz] * vt[iy][ix][iz];
		}
	    }
	}
	/* circulate wavefield arrays */
	ut=um;
	um=uo;
	uo=up;
	up=ut;
	
	if(dabc) {
	    /* one-way abc apply */
	    abcone3d_apply(uo,um,NOP,abc,fdm);
	    sponge3d_apply(um,spo,fdm);
	    sponge3d_apply(uo,spo,fdm);
	    sponge3d_apply(up,spo,fdm);
	}

	/* extract data */
	lint3d_extract(uo,dd,cr);

	if(snap && it%jsnap==0) {
	    cut3d(uo,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
	}
	if(        it%jdata==0) 
	    sf_floatwrite(dd,nr,Fdat);
    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(**um); free(*um); free(um);
    free(**up); free(*up); free(up);
    free(**uo); free(*uo); free(uo);
    free(**ua); free(*ua); free(ua);
    free(**uc); free(*uc); free(uc);

    free(**rox); free(*rox); free(rox);
    free(**roy); free(*roy); free(roy);
    free(**roz); free(*roz); free(roz);

    free(**vp); free(*vp); free(vp);
    free(**vt); free(*vt); free(vt);

    free(ss);
    free(rr);
    free(dd);
    free(ww);

    exit (0);
}
