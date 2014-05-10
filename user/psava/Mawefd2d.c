/* 2D acoustic time-domain FD modeling
4th order in space, 2nd order in time. Absorbing boundary conditions.
*/
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
#include <time.h>
/* check: dt<= 0.2 * min(dx,dz)/vmin */

#define NOP 2 /* derivative operator half-size */

#define C0 -2.500000 /*    c0=-30./12.; */
#define CA +1.333333 /*    ca=+16./12.; */
#define CB -0.083333 /*    cb=- 1./12.; */

/* LS-optimized coefficients */
#define F1  +1.16303535		
#define F2  -0.05686715	

/* FD derivative stencils */
#define FX(a,ix,iz,s) (F1*(a[ix  ][iz  ] - a[ix-1][iz  ]) + \
                       F2*(a[ix+1][iz  ] - a[ix-2][iz  ]))*s
#define FZ(a,ix,iz,s) (F1*(a[ix  ][iz  ] - a[ix  ][iz-1]) + \
                       F2*(a[ix  ][iz+1] - a[ix  ][iz-2]))*s

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc,cden,adj; 
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
    float **iro=NULL;			/* inverse density */
    float **uat=NULL;          /* auxiliary wavefield */
    float **vp=NULL;           /* velocity */
    float **vt=NULL;           /* temporary vp*vp * dt*dt */

    float **um,**uo,**up,**ua,**ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float co,cax,cbx,caz,cbz,f1x,f2x,f1z,f2z;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL;

    /* for benchmarking */
    clock_t start_t, end_t;
    float total_t;

    float ox;
    float mx;
    float oz;
    float mz;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* Wavefield snapshots flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* Multiple sources, one wvlt */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
    if(! sf_getbool("cden",&cden)) cden=false; /* Constant density */
    if(! sf_getbool("adj", &adj))   adj=false; /* adjoint flag */

    if(! sf_getbool("free",&fsrf)) fsrf=false; /* Free surface flag */
    if(! sf_getbool("fsrf",&fsrf)) fsrf=false; /* Free surface flag */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fdat = sf_output("out"); /* data      */
    if( snap) Fwfl = sf_output("wfl"); /* wavefield */
    if(!cden) Fden = sf_input ("den"); /* density   */

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
    /* # of t steps at which to save receiver data */
    if(snap) {
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;        
        /* # of t steps at which to save wavefield */ 
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
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az); /* Saved wfld window nz */
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax); /* Saved wfld window nx */

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az); /* Saved wfld window oz */
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax); /* Saved wfld window ox */

	if(!sf_getfloat("dqz",&dqz)) dqz=sf_d(az); /* Saved wfld window dz */
	if(!sf_getfloat("dqx",&dqx)) dqx=sf_d(ax); /* Saved wfld window dx */

	acz = sf_maxa(nqz,oqz,dqz); if(verb) sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); if(verb) sf_raxa(acx);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

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
    /*( nb=2 boundary padding in grid points )*/

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/

    if(expl) ww = sf_floatalloc( 1);
    else     ww = sf_floatalloc(ns);
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    /*------------------------------------------------------------*/
    /* coordinate check: if point is outside the grid, an error is thrown */
    ox = sf_o(ax); mx = ox + (sf_n(ax)-1)*sf_d(ax);
    oz = sf_o(az); mz = oz + (sf_n(az)-1)*sf_d(az);

    for (ix=0; ix<ns; ++ix){
      if (ss[ix].x < ox || ss[ix].x > mx){
        sf_error("fatal error: source coordinate x is outside the grid for point %d",ix);      
      }
      if (ss[ix].z < oz || ss[ix].z > mz){
        sf_error("fatal error: source coordinate z is outside the grid for point %d",ix);      
      }
    }

    for (ix=0; ix<nr; ++ix){
      if (rr[ix].x < ox || rr[ix].x > mx){
        sf_error("fatal error: coordinate x is outside the grid for point %d",ix);      
      }
      if (rr[ix].z < oz || rr[ix].z > mz){
        sf_error("fatal error: receiver coordinate z is outside the grid for point %d",ix);      
      }
    }
    /* End of coordinate consistency check */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);
    fdbell_init(5);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;

    co = C0 * (idx*idx+idz*idz);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    f1x = F1*idx;
    f2x = F2*idx;
    
    f1z = F1*idz;
    f2z = F2*idz;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 

    vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

    if (!cden) {

        /* input density */
        ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        iro =sf_floatalloc2(fdm->nzpad,fdm->nxpad);        
        sf_floatread(tt[0],nz*nx,Fden); expand(tt,ro ,fdm);
	
        /* auxiliary vector */
        uat =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	
	/*inverse density, to avoid division in the extrapolation */
	iro[0][0] = 1/ro[0][0];
	for    (ix=1; ix<fdm->nxpad; ix++){
	    for(iz=1; iz<fdm->nzpad; iz++){
		iro[ix][iz] = 4./( 2*ro[ix  ][iz  ] + 
				     ro[ix-1][iz  ] + 
				     ro[ix  ][iz-1]);
	    }
	}
	
    }

    /* input velocity */
    sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vp,fdm);
    /* precompute vp^2 * dt^2 */
    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    vt[ix][iz] = vp[ix][iz] * vp[ix][iz] * dt*dt;
	}
    }
    if(fsrf) { /* free surface */
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nb; iz++) {
		vt[ix][iz]=0;
	    }
	}
    }

    free(*tt); free(tt);    

    /*------------------------------------------------------------*/

    if(dabc) {
	/* one-way abc setup */
	abc = abcone2d_make(NOP,dt,vp,fsrf,fdm);
	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    free(*vp);  free(vp);

    /*------------------------------------------------------------*/

    /* allocate wavefield arrays */
    um=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    up=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ua[ix][iz]=0;
	}
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if (cden){
	/* begin constant density */
	
	if(verb) fprintf(stderr,"\n");
	start_t = clock();
	for (it=0; it<nt; it++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	    
#ifdef _OPENMP
#pragma omp parallel						\
    private(ix,iz)						\
    shared(fdm,ua,uo,co,cax,caz,cbx,cbz,idx,idz)
#endif
	    {	
#ifdef _OPENMP
#pragma omp for	schedule(dynamic,fdm->ompchunk)	
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
		    }
		}
		
		/* step forward in time */
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		for    (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nzpad; iz++) {
			up[ix][iz] = 2*uo[ix][iz] 
			    -          um[ix][iz] 
			    +          ua[ix][iz] * vt[ix][iz];
		    }
		}
		
	    } /* end parallel section */
	    
	    /* inject displacement source */
	    if(expl) {
		if(adj) sf_seek(Fwav,(off_t)(nt-1-it)*1 *sizeof(float),SEEK_SET);
		sf_floatread(ww, 1,Fwav);
		lint2d_bell1(up,ww[0],cs);
	    } else {
		if(adj) sf_seek(Fwav,(off_t)(nt-1-it)*ns*sizeof(float),SEEK_SET);
		sf_floatread(ww,ns,Fwav);
		lint2d_bell(up,ww,cs);
	    }
	    
	    /* extract data at receivers */
	    lint2d_extract(up,dd,cr);
	    if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);
	    
	    /* extract wavefield in the "box" */
	    if(snap && it%jsnap==0) {
		cut2d(up,uc,fdm,acz,acx);
		sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
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
	    }
	    
	}	/* end time loop */
	end_t = clock();
	if(verb) fprintf(stderr,"\n");
	
	if (verb){	
	    total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
	    fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
	    fprintf(stderr,"Exiting of the program...\n");
	}

	/* end constant density section */

    } else {

	uat=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	if(verb) fprintf(stderr,"\n");
	start_t = clock();
	    for (it=0; it<nt; it++) {
		if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
		
#ifdef _OPENMP
#pragma omp parallel private(ix,iz)	
#endif
		{		
		    /* Z spatial derivatives */
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		
			    
			    /* gather */
			    uat[ix][iz]  = iro[ix][iz]*(
				f1z*(uo[ix  ][iz  ] - uo[ix  ][iz-1]) +
				f2z*(uo[ix  ][iz+1] - uo[ix  ][iz-2]));
			    
			}
		    }
		    
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		    for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    /* scatter */
			    ua[ix][iz  ]  =    
				f1z*(uat[ix][iz  ] - 
				     uat[ix][iz+1]) + 
				f2z*(uat[ix][iz-1] -
				     uat[ix][iz+2]);
			}
		    }
	
		    /* X spatial derivatives */	    
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		    for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    /* gather */
			    uat[ix][iz]  = iro[ix][iz]*(
				f1x*(uo[ix  ][iz  ] - uo[ix-1][iz  ]) +
				f2x*(uo[ix+1][iz  ] - uo[ix-2][iz  ]));
			    
			}
		    }
		    
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		    for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    /* scatter */
			    ua[ix  ][iz]  +=    
				f1x*(uat[ix  ][iz] -
				     uat[ix+1][iz]) +
				f2x*(uat[ix-1][iz] -
				     uat[ix+2][iz]);
			}
		    }
		    
		    /* step forward in time */
#ifdef _OPENMP
#pragma omp for schedule(dynamic,fdm->ompchunk)
#endif
		    for    (ix=0; ix<fdm->nxpad; ix++) {
			for(iz=0; iz<fdm->nzpad; iz++) {
			    up[ix][iz] = 2*uo[ix][iz] 
				-  um[ix][iz] 
				-  ro[ix][iz]*vt[ix][iz]*ua[ix][iz];
			}
		    }	
		}	/* end parallel section */

		/* inject displacement source */
		if(expl) {
		    if(adj) sf_seek(Fwav,(off_t)(nt-1-it)*1 *sizeof(float),SEEK_SET);
		    sf_floatread(ww, 1,Fwav);
		    lint2d_bell1(up,ww[0],cs);
		} else {
		    if(adj) sf_seek(Fwav,(off_t)(nt-1-it)*ns*sizeof(float),SEEK_SET);
		    sf_floatread(ww,ns,Fwav);	
		    lint2d_bell(up,ww,cs);
		}
		
		/* extract data at receivers */
		lint2d_extract(up,dd,cr);
		if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);
		
		/* extract wavefield in the "box" */
		if(snap && it%jsnap==0) {
		    cut2d(up,uc,fdm,acz,acx);
		    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
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
		}
		
	    }	/* end time loop*/
	    end_t = clock();
	    if(verb) fprintf(stderr,"\n");
	    
	    if (verb){	
		total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
		fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
		fprintf(stderr,"Exiting of the program...\n");
	    }	

	    /* end variable density section*/
    }
        
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*um); free(um);
    free(*up); free(up);
    free(*uo); free(uo);
    free(*ua); free(ua);
    if(snap) {
        free(*uc); free(uc);
    }
    
    if(!cden) {
        free(*uat); free(uat);
        free(*ro); free(ro);
        free(*iro); free(iro);        
    }
    
    free(*vt);  free(vt);
    
    free(ww);
    free(ss);
    free(rr);
    free(dd);
    
    /*------------------------------------------------------------*/ 
    /* CLOSE FILES AND EXIT */
    if (Fwav!=NULL) sf_fileclose(Fwav); 
    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Frec!=NULL) sf_fileclose(Frec);
    if (Fvel!=NULL) sf_fileclose(Fvel);
    if (Fden!=NULL) sf_fileclose(Fden);
    if (Fdat!=NULL) sf_fileclose(Fdat);
    if (Fwfl!=NULL) sf_fileclose(Fwfl);

    exit (0);
}
