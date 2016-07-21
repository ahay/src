/* 3D elastic time-domain pseudo-spectral (k-space) modeling using KISS-FFT
   sou wavelet  (nx,ny,nc,nt)
   rec data     (nx,ny,nc,nt)
   sou geometry (nc,nx,ny)
   rec geometry (nc,nx,ny)
*/
/*
  Copyright (C) 2016 The University of Texas at Austin
  
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

  The main framework is following Paul Sava's sfewefd3d program.
  This program uses pseudo-spectral method to calculate spatial derivatives.
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ksutil.h"

int main(int argc, char* argv[])
{
    /*------------------------------------------------------------*/
    /* Execution control, I/O files and geometry                  */
    /*------------------------------------------------------------*/
    bool verb,fsrf,snap,dabc,opot,kspace,back,secd; /* execution flags */
    int  jsnap,ntsnap,jdata; /* jump along axes */
    int  ssou; /* source type */
    int  shft; /* time shift for wavefield matching in RTM */

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,ax,ay,az; /* time, x, y, z */ 
    sf_axis asx,asy,arx,ary,ac;    /* sou, rec-x, rec-y, component */ 

    /* dimension, index and interval */
    int     nt,nz,nx,ny,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy,idz,idx,idy;

    /* FDM and KSP structure */ //!!!JS
    fdm3d    fdm=NULL;
    dft3d    dft=NULL;
    ksp3d    ksp=NULL;
    abcone3d /* abcp=NULL, */ abcs=NULL;
    sponge   spo=NULL;

    /* I/O arrays for sou & rec */
    float***ww=NULL;           /* wavelet   */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    float **dd=NULL;           /* data      */

    /*------------------------------------------------------------*/
    /* model parameter arays                                      */
    /*------------------------------------------------------------*/
    float ***tt=NULL;          /* temporary array */
    float ***ro=NULL;          /* density */

    /* orthorombic footprint - 9 coefficients */
    /* c11 c12 c13 
       .   c22 c23 
       .   .   c33 
                  c44
                     c55
                        c66 */
    float ***c11=NULL;
    float ***c22=NULL;
    float ***c33=NULL;
    float ***c44=NULL;
    float ***c55=NULL;
    float ***c66=NULL;
    float ***c12=NULL;
    float ***c13=NULL;
    float ***c23=NULL;
    float ***vp,***vs; /* velocity */
    float ***qp=NULL,***qsx=NULL,***qsy=NULL,***qsz=NULL; /* potentials */ //!!!JS

    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1       */
    /*------------------------------------------------------------*/
    float ***umz,***uoz,***upz,***uaz,***utz;  /* dis*3, acc, tmp (for rotation) */
    float ***umx,***uox,***upx,***uax,***utx;
    float ***umy,***uoy,***upy,***uay,***uty;

    /* stress/strain tensor */ 
    float ***tzz,***txx,***tyy,***txy,***tyz,***tzx; /* strain then stress (in-place) */
    float    szz,   sxx,   syy,   sxy,   syz,   szx; /* tmp var for storing stress */

    /*------------------------------------------------------------*/
    /* spatial derivatives from pseudo-spectral method            */
    /*------------------------------------------------------------*/
    float ***xdx, ***xdy, ***xdz;
    float ***ydx, ***ydy, ***ydz;
    float ***zdx, ***zdy, ***zdz;

    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices                       */
    /*------------------------------------------------------------*/
    lint3d cs,cr; /* for injecting source and extracting data */

    /* Gaussian bell */
    int nbell;
    
    /*------------------------------------------------------------*/
    /* wavefield cut params                                       */
    /*------------------------------------------------------------*/
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***uc=NULL; /* tmp array for output wavefield snaps */

    /*------------------------------------------------------------*/
    /* init RSF                                                   */
    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters                                             */
    /*------------------------------------------------------------*/
#ifdef _OPENMP
    omp_init();
#endif

    /*------------------------------------------------------------*/
    /* read execution flags                                       */
    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getint( "ssou",&ssou)) ssou=0;     /* 0 -> acceleration source; 1 -> stress source; 2 -> displacement source */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    if(! sf_getbool("opot",&opot)) opot=false; /* output potentials -> 1*scalar, 3*vector potentials */
    if(! sf_getbool("back",&back)) back=false; /* backward extrapolation flag (for rtm) */
    if(! sf_getbool("kspace",&kspace)) kspace=false; /* k-space method (ps) flag */
    if(! sf_getbool("secd",&secd)) secd=false; /* second order weq flag */
    if (secd) {
        if (!kspace || ssou==1) sf_error("This option is not yet supported!");
    }

    /*------------------------------------------------------------*/
    /* I/O files                                                  */
    /*------------------------------------------------------------*/
    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    /*------------------------------------------------------------*/
    /* axes                                                       */
    /*------------------------------------------------------------*/
    at = sf_iaxa(Fwav,4); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space y */

    asx = sf_iaxa(Fsou,2); sf_setlabel(asx,"sx"); if(verb) sf_raxa(asx); /* sources x */
    asy = sf_iaxa(Fsou,3); sf_setlabel(asy,"sy"); if(verb) sf_raxa(asy); /* sources y */
    arx = sf_iaxa(Frec,2); sf_setlabel(arx,"rx"); if(verb) sf_raxa(arx); /* receivers x */
    ary = sf_iaxa(Frec,3); sf_setlabel(ary,"ry"); if(verb) sf_raxa(ary); /* receivers y */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    ns = sf_n(asx)*sf_n(asy);
    nr = sf_n(arx)*sf_n(ary);

    /*------------------------------------------------------------*/
    /* other execution parameters                                 */
    /*------------------------------------------------------------*/
    if(! sf_getint("nbell",&nbell)) nbell=5;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }
    if(back) {
        shft = (nt-1)%jsnap;
        sf_warning("For backward extrapolation, make sure nbell(%d)=0 and ssou(%d)=2",nbell,ssou);
    } else shft = 0;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC                     */
    /*------------------------------------------------------------*/
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP; //!!!JS

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    if(nbell) fdbell3d_init(nbell);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);

    /*------------------------------------------------------------*/
    /* 3D vector components                                       */
    /*------------------------------------------------------------*/
    nc=3;
    if(opot) { /* output 1 scalar potential and 3 vector potentials */
	ac=sf_maxa(nc+1,0,1);
    } else {   /* output 3 cartesian components */
	ac=sf_maxa(nc  ,0,1);
    }

    /*------------------------------------------------------------*/
    /* setup output data header                                   */
    /*------------------------------------------------------------*/
    sf_oaxa(Fdat,arx,1);
    sf_oaxa(Fdat,ary,2);
    sf_oaxa(Fdat,ac,3);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,4);

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
	/* TODO: check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0; /* ntsnap = it/jsnap+1; */
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,acy,3);
	sf_oaxa(Fwfl,ac, 4);
	sf_oaxa(Fwfl,at, 5);
    }

    /*------------------------------------------------------------*/
    /* source and data array                                      */
    /*------------------------------------------------------------*/
    ww=sf_floatalloc3(ns,nc,nt); /* Fast axis: n_sou > n_comp > n_time */
    sf_floatread(ww[0][0],nt*nc*ns,Fwav);

    if(opot) {
	dd=sf_floatalloc2(nr,nc+1);
    } else {
	dd=sf_floatalloc2(nr,nc  );
    }

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates                          */
    /*------------------------------------------------------------*/
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    /* calculate 3d linear interpolation coef for sou & rec */
    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup derivative coefficients                              */
    /*------------------------------------------------------------*/
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;

    /*------------------------------------------------------------*/ 
    /* allocation I/O arrays                                      */
    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 
    
    ro =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    c11=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c22=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c33=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c44=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c55=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c66=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c12=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c13=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    c23=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);     

    /* input density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden);     expand3d(tt,ro ,fdm);

    /* input stiffness */
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,fdm);    
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,fdm);
    sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,fdm);

    free(**tt); free(*tt); free(tt);

    /*------------------------------------------------------------*/
    /* setup absorbing boundary condition                         */
    /*------------------------------------------------------------*/
    if (dabc) {
	/* one-way abc setup   */
	vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	vs = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    vp[iy][ix][iz] = sqrt( c11[iy][ix][iz]/ro[iy][ix][iz] );
		    vs[iy][ix][iz] = sqrt( c55[iy][ix][iz]/ro[iy][ix][iz] );
		    /* vs[iy][ix][iz] = sqrt( c13[iy][ix][iz]/ro[iy][ix][iz] ); */
		}
	    }
	}
	/* abcp = abcone3d_make(NOP,dt,vp,fsrf,fdm); */
	abcs = abcone3d_make(NOP,dt,vs,fsrf,fdm);
	free(**vp); free(*vp); free(vp);
	free(**vs); free(*vs); free(vs);

	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 for updating displacement in time   */
    /*------------------------------------------------------------*/
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		ro[iy][ix][iz] = dt*dt/ro[iy][ix][iz];
	    }
	}
     }

    /*------------------------------------------------------------*/
    /* allocate and initialize wavefield arrays                   */
    /*------------------------------------------------------------*/
    /* z-component */
    umz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uaz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* x-component */
    umx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uax=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* y-component */
    umy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uay=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* stress/strain tensor */
    tzz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tzx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* derivative array for k-space method */
    if (kspace) {
        dft = dft3d_init(1,true,false,fdm);
        ksp = ksp3d_make(dft);
        xdx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        xdy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        xdz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        ydx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        ydy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        ydz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        zdx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        zdy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
        zdz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    }

    /* initialize to zero */
    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		umz[iy][ix][iz]=0; umx[iy][ix][iz]=0; umy[iy][ix][iz]=0;
		uoz[iy][ix][iz]=0; uox[iy][ix][iz]=0; uoy[iy][ix][iz]=0;
		upz[iy][ix][iz]=0; upx[iy][ix][iz]=0; upy[iy][ix][iz]=0;
		uaz[iy][ix][iz]=0; uax[iy][ix][iz]=0; uay[iy][ix][iz]=0;
	    }
	}
    }

    /* allocate arrays for calculating 4 potentials */
    if(opot) {
	qp =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qsx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qsy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	qsz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    }

    /*------------------------------------------------------------*/ 
    /*------------------------ MAIN LOOP -------------------------*/ 
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
        if(verb) sf_warning("it=%d/%d;",it,nt); /*fprintf(stderr,"\b\b\b\b\b%d",it);*/

        if (secd) {
            ksp3d_apply2(txx[0][0],tyy[0][0],tzz[0][0],txy[0][0],tyz[0][0],tzx[0][0],uox[0][0],dft,ksp);
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,c11,c22,c33,c44,c55,c66,c12,c13,c23)
#endif
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {
                        uax[iy][ix][iz] = c11[iy][ix][iz]*txx[iy][ix][iz] + c66[iy][ix][iz]*tyy[iy][ix][iz] + c55[iy][ix][iz]*tzz[iy][ix][iz];
                        uay[iy][ix][iz] = (c12[iy][ix][iz] + c66[iy][ix][iz])*txy[iy][ix][iz];
                        uaz[iy][ix][iz] = (c13[iy][ix][iz] + c55[iy][ix][iz])*tzx[iy][ix][iz];
                    }
                }
            }

            ksp3d_apply2(txx[0][0],tyy[0][0],tzz[0][0],txy[0][0],tyz[0][0],tzx[0][0],uoy[0][0],dft,ksp);
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,c11,c22,c33,c44,c55,c66,c12,c13,c23)
#endif
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {
                        uax[iy][ix][iz] += (c12[iy][ix][iz] + c66[iy][ix][iz])*txy[iy][ix][iz];
                        uay[iy][ix][iz] += c66[iy][ix][iz]*txx[iy][ix][iz] + c22[iy][ix][iz]*tyy[iy][ix][iz] + c44[iy][ix][iz]*tzz[iy][ix][iz];
                        uaz[iy][ix][iz] += (c23[iy][ix][iz] + c44[iy][ix][iz])*tyz[iy][ix][iz];
                    }
                }
            }

            ksp3d_apply2(txx[0][0],tyy[0][0],tzz[0][0],txy[0][0],tyz[0][0],tzx[0][0],uoz[0][0],dft,ksp);
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,c11,c22,c33,c44,c55,c66,c12,c13,c23)
#endif
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {
                        uax[iy][ix][iz] += (c13[iy][ix][iz] + c55[iy][ix][iz])*tzx[iy][ix][iz];
                        uay[iy][ix][iz] += (c23[iy][ix][iz] + c44[iy][ix][iz])*tyz[iy][ix][iz];
                        uaz[iy][ix][iz] += c55[iy][ix][iz]*txx[iy][ix][iz] + c44[iy][ix][iz]*tyy[iy][ix][iz] + c33[iy][ix][iz]*tzz[iy][ix][iz];
                    }
                }
            }

        } else {

	/*------------------------------------------------------------*/
	/* from displacement to strain                                */
        /*------------------------------------------------------------*/
        /* 
	 * exx = Fx(ux)
	 * eyy = Fy(uy)
	 * ezz = Fz(uz)
	 * exy = By(ux) + Bx(uy)
	 * eyz = Bz(uy) + By(uz)
	 * ezx = Bx(uz) + Bz(ux)
	 */
        if (!kspace) { /* finite-difference */
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,idx,idy,idz)
#endif
            for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
                for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
                    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {

                        txx[iy][ix][iz] = Dx(uox,ix,iy,iz,idx);
                        tyy[iy][ix][iz] = Dy(uoy,ix,iy,iz,idy);
                        tzz[iy][ix][iz] = Dz(uoz,ix,iy,iz,idz);

                        txy[iy][ix][iz] = Dy(uox,ix,iy,iz,idy) + Dx(uoy,ix,iy,iz,idx);
                        tyz[iy][ix][iz] = Dz(uoy,ix,iy,iz,idz) + Dy(uoz,ix,iy,iz,idy);
                        tzx[iy][ix][iz] = Dx(uoz,ix,iy,iz,idx) + Dz(uox,ix,iy,iz,idz);
                    }
                }
            }
        } else { /* pseudo-spectral */
            ksp3d_apply(xdx[0][0],xdy[0][0],xdz[0][0],uox[0][0],dft,ksp);
            ksp3d_apply(ydx[0][0],ydy[0][0],ydz[0][0],uoy[0][0],dft,ksp);
            ksp3d_apply(zdx[0][0],zdy[0][0],zdz[0][0],uoz[0][0],dft,ksp);
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,xdx,xdy,xdz,ydx,ydy,ydz,zdx,zdy,zdz)
#endif
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {

                        txx[iy][ix][iz] = xdx[iy][ix][iz];
                        tyy[iy][ix][iz] = ydy[iy][ix][iz];
                        tzz[iy][ix][iz] = zdz[iy][ix][iz];

                        txy[iy][ix][iz] = xdy[iy][ix][iz] + ydx[iy][ix][iz];
                        tyz[iy][ix][iz] = ydz[iy][ix][iz] + zdy[iy][ix][iz];
                        tzx[iy][ix][iz] = zdx[iy][ix][iz] + xdz[iy][ix][iz];
                    }
                }
            }
        }

	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,c11,c22,c33,c44,c55,c66,c12,c13,c23)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    
		    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
			+ c12[iy][ix][iz] * tyy[iy][ix][iz]
			+ c13[iy][ix][iz] * tzz[iy][ix][iz];
		    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
			+ c22[iy][ix][iz] * tyy[iy][ix][iz]
			+ c23[iy][ix][iz] * tzz[iy][ix][iz];
		    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
			+ c23[iy][ix][iz] * tyy[iy][ix][iz]
			+ c33[iy][ix][iz] * tzz[iy][ix][iz];
		    
		    sxy = c66[iy][ix][iz] * txy[iy][ix][iz];
		    syz = c44[iy][ix][iz] * tyz[iy][ix][iz];
		    szx = c55[iy][ix][iz] * tzx[iy][ix][iz];
		    
		    txx[iy][ix][iz] = sxx;
		    tyy[iy][ix][iz] = syy;
		    tzz[iy][ix][iz] = szz;

		    txy[iy][ix][iz] = sxy;
		    tyz[iy][ix][iz] = syz;
		    tzx[iy][ix][iz] = szx;
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	if(fsrf) {
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx)
#endif
	    for        (iy=0; iy<fdm->nypad; iy++) {
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nb;    iz++) {
			txx[iy][ix][iz]=0;
			tyy[iy][ix][iz]=0;
			tzz[iy][ix][iz]=0;

			txy[iy][ix][iz]=0;
			tyz[iy][ix][iz]=0;
			tzx[iy][ix][iz]=0;
		    }
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(ssou==1) {
            if(nbell) {
                lint3d_bell(tzz,ww[it][0],cs);
                lint3d_bell(txx,ww[it][1],cs);
                lint3d_bell(tyy,ww[it][2],cs);
            } else {
                lint3d_inject(tzz,ww[it][0],cs);
                lint3d_inject(txx,ww[it][1],cs);
                lint3d_inject(tyy,ww[it][2],cs);
            }
	}
	
	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* 
	 * ax = Bx(txx) + Fy(txy) + Fz(txz)
	 * ay = Fx(txy) + By(tyy) + Fz(tyz)
	 * az = Fx(txz) + Fy(tyz) + Bz(tzz)
	 */	
        if (!kspace) { /* finite-difference */
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,idx,idy,idz)
#endif
            for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
                for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
                    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		    
                        uax[iy][ix][iz] = Dx( txx,ix,iy,iz,idx ) + Dy( txy,ix,iy,iz,idy ) + Dz( tzx,ix,iy,iz,idz ) ;
                        uay[iy][ix][iz] = Dx( txy,ix,iy,iz,idx ) + Dy( tyy,ix,iy,iz,idy ) + Dz( tyz,ix,iy,iz,idz ) ;
                        uaz[iy][ix][iz] = Dx( tzx,ix,iy,iz,idx ) + Dy( tyz,ix,iy,iz,idy ) + Dz( tzz,ix,iy,iz,idz ) ;		    
                    }
                }
            }
        } else { /* pseudo-spectral */
            ksp3d_apply1(xdx[0][0],xdy[0][0],xdz[0][0],txx[0][0],txy[0][0],tzx[0][0],dft,ksp);
            ksp3d_apply1(ydx[0][0],ydy[0][0],ydz[0][0],txy[0][0],tyy[0][0],tyz[0][0],dft,ksp);
            ksp3d_apply1(zdx[0][0],zdy[0][0],zdz[0][0],tzx[0][0],tyz[0][0],tzz[0][0],dft,ksp);
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,uax,uay,uaz,xdx,xdy,xdz,ydx,ydy,ydz,zdx,zdy,zdz)
#endif
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {		    
                        uax[iy][ix][iz] = xdx[iy][ix][iz] + xdy[iy][ix][iz] + xdz[iy][ix][iz]; 
                        uay[iy][ix][iz] = ydx[iy][ix][iz] + ydy[iy][ix][iz] + ydz[iy][ix][iz];
                        uaz[iy][ix][iz] = zdx[iy][ix][iz] + zdy[iy][ix][iz] + zdz[iy][ix][iz];
                    }
                }
            }
        }

        } /* if secd */

	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(ssou==0) {
            if(nbell) {
                lint3d_bell(uaz,ww[it][0],cs);
                lint3d_bell(uax,ww[it][1],cs);
                lint3d_bell(uay,ww[it][2],cs);
            } else {
                lint3d_inject(uaz,ww[it][0],cs);
                lint3d_inject(uax,ww[it][1],cs);
                lint3d_inject(uay,ww[it][2],cs);
            }
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,uox,uoy,uoz,umx,umy,umz,upx,upy,upz,uax,uay,uaz,ro)
#endif
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    upx[iy][ix][iz] = 2*uox[iy][ix][iz] 
			-               umx[iy][ix][iz] 
			+               uax[iy][ix][iz] * ro[iy][ix][iz]; 

		    upy[iy][ix][iz] = 2*uoy[iy][ix][iz] 
			-               umy[iy][ix][iz] 
			+               uay[iy][ix][iz] * ro[iy][ix][iz]; 

		    upz[iy][ix][iz] = 2*uoz[iy][ix][iz] 
			-               umz[iy][ix][iz] 
			+               uaz[iy][ix][iz] * ro[iy][ix][iz]; 
		    
		}
	    }
	}
	/* circulate wavefield arrays */
	utz=umz; uty=umy; utx=umx;
	umz=uoz; umy=uoy; umx=uox;
	uoz=upz; uoy=upy; uox=upx;
	upz=utz; upy=uty; upx=utx;

	/*------------------------------------------------------------*/
	/* absorbing boundary condition                               */
	/*------------------------------------------------------------*/
	if(dabc) {
	    /* one-way ABC */
            /*
	    abcone3d_apply(uoz,umz,NOP,abcp,fdm);
	    abcone3d_apply(uox,umx,NOP,abcp,fdm);
	    abcone3d_apply(uoy,umy,NOP,abcp,fdm);
            */
	    
	    abcone3d_apply(uoz,umz,NOP,abcs,fdm);
	    abcone3d_apply(uox,umx,NOP,abcs,fdm);
	    abcone3d_apply(uoy,umy,NOP,abcs,fdm);

	    /* sponge ABC */
	    sponge3d_apply(umz,spo,fdm);
	    sponge3d_apply(uoz,spo,fdm);
	    
	    sponge3d_apply(umx,spo,fdm);
	    sponge3d_apply(uox,spo,fdm);

	    sponge3d_apply(umy,spo,fdm);
	    sponge3d_apply(uoy,spo,fdm);
	}	    

        /*------------------------------------------------------------*/
	/* inject displacement source                                 */
	/*------------------------------------------------------------*/
	if(ssou==2) {
            if(nbell) {
                lint3d_bell(uoz,ww[it][0],cs);
                lint3d_bell(uox,ww[it][1],cs);
                lint3d_bell(uoy,ww[it][2],cs);
            } else {
                lint3d_inject(uoz,ww[it][0],cs);
                lint3d_inject(uox,ww[it][1],cs);
                lint3d_inject(uoy,ww[it][2],cs);
            }
	}

	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
	if(opot) {
		
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iy,iz)				\
    shared(fdm,uox,uoy,uoz,idx,idy,idz)
#endif
	    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
		for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
			
			qp [iy][ix][iz] = Dx( uox,ix,iy,iz,idx )
			    +             Dy( uoy,ix,iy,iz,idy )
			    +             Dz( uoz,ix,iy,iz,idz );
			
			qsx[iy][ix][iz] = Dy( uoz,ix,iy,iz,idy ) - Dz( uoy,ix,iy,iz,idz );
			qsy[iy][ix][iz] = Dz( uox,ix,iy,iz,idz ) - Dx( uoz,ix,iy,iz,idx );
			qsz[iy][ix][iz] = Dx( uoy,ix,iy,iz,idx ) - Dy( uox,ix,iy,iz,idy );
		    }
		}
	    }

	    if(snap && (it-shft)%jsnap==0) {
		cut3d(qp ,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

		cut3d(qsz,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

		cut3d(qsx,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

		cut3d(qsy,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
	    }
	    
	    lint3d_extract(qp , dd[0],cr);
	    lint3d_extract(qsx, dd[1],cr);
	    lint3d_extract(qsy, dd[2],cr);
	    lint3d_extract(qsz, dd[3],cr);
	    if(it%jdata==0) sf_floatwrite(dd[0],nr*(nc+1),Fdat);

	} else {

	    if(snap && (it-shft)%jsnap==0) {
		cut3d(uoz,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
		
		cut3d(uox,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
		
		cut3d(uoy,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
	    }
	    
	    lint3d_extract(uoz,dd[0],cr);
	    lint3d_extract(uox,dd[1],cr);
	    lint3d_extract(uoy,dd[2],cr);
	    if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
	}

    }
    if(verb) sf_warning(".");
    if(verb) fprintf(stderr,"\n");    
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    if (kspace) {
        dft3d_finalize();
        free(dft);
        ksp3d_finalize();
        free(ksp);
        free(**xdx); free(*xdx); free(xdx);
        free(**xdy); free(*xdy); free(xdy);
        free(**xdz); free(*xdz); free(xdz);
        free(**ydx); free(*ydx); free(ydx);
        free(**ydy); free(*ydy); free(ydy);
        free(**ydz); free(*ydz); free(ydz);
        free(**zdx); free(*zdx); free(zdx);
        free(**zdy); free(*zdy); free(zdy);
        free(**zdz); free(*zdz); free(zdz);
    }

    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);

    free(**ro);  free(*ro);  free(ro);
    free(**c11); free(*c11); free(c11);
    free(**c22); free(*c22); free(c22);
    free(**c33); free(*c33); free(c33);
    free(**c44); free(*c44); free(c44);
    free(**c55); free(*c55); free(c55);
    free(**c66); free(*c66); free(c66);
    free(**c12); free(*c12); free(c12);
    free(**c13); free(*c13); free(c13);
    free(**c23); free(*c23); free(c23);

    free(**umz); free(*umz); free(umz);
    free(**uoz); free(*uoz); free(uoz);
    free(**upz); free(*upz); free(upz);
    free(**uaz); free(*uaz); free(uaz);

    free(**umx); free(*umx); free(umx);
    free(**uox); free(*uox); free(uox);
    free(**upx); free(*upx); free(upx);
    free(**uax); free(*uax); free(uax);

    free(**umy); free(*umy); free(umy);
    free(**uoy); free(*uoy); free(uoy);
    free(**upy); free(*upy); free(upy);
    free(**uay); free(*uay); free(uay);

    free(**tzz); free(*tzz); free(tzz);
    free(**txx); free(*txx); free(txx);
    free(**tyy); free(*tyy); free(tyy);
    free(**txy); free(*txy); free(txy);
    free(**tyz); free(*tyz); free(tyz);
    free(**tzx); free(*tzx); free(tzx);

    if (snap) {
       free(**uc);  free(*uc);  free(uc);    
    }

    if(opot) {
	free(**qp);  free(*qp);  free(qp);    
	free(**qsx); free(*qsx); free(qsx);
	free(**qsy); free(*qsy); free(qsy);
	free(**qsz); free(*qsz); free(qsz);
    }
    /*------------------------------------------------------------*/


    exit (0);
}

