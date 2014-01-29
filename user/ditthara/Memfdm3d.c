/* 3D Electromagnetic time-domain FD modeling

4th order in space, 2nd order in time. Absorbing boundary conditions. */
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
#endif

#include "fdutil.h"
/* check: dt<= 0.2 * min(dx,dz)/vmin */

#define NOP 2 /* derivative operator half-size */

#define C0 -2.500000 /*    c0=-30./12.; */
#define CA +1.333333 /*    ca=+16./12.; */
#define CB -0.083333 /*    cb=- 1./12.; */

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc; 
    int  jsnap,ntsnap,jdata;

    /* I/O files */

    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */

    sf_file Fmag=NULL; /* magnetic permitivity */
    sf_file Fele=NULL; /* electric susceptibility */
    sf_file Fcdt=NULL; /* conductivity */

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

    float ***tt=NULL;

    float ***vel=NULL;         /* velocity */
    float ***mag=NULL;
    float ***ele=NULL;
    float ***cdt=NULL;

    float ***cdtele=NULL;      /* temporary cdt*dt/2*ele */
    float ***magele=NULL;	   /* temporary dt*dt/mag*ele */

    float ***vt=NULL;	       /* temporary dt*dt/vel*vel */

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
    omp_init();
#endif
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* Wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* Free surface flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* Multiple sources, one wvlt */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */

    if(snap) Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    /*------------------------------------------------------------*/
    Fmag = sf_input ("mag"); /* magnetic permitivity */
    Fele = sf_input ("ele"); /* electric susceptibility */
    Fcdt = sf_input ("cdt"); /* conductivity */
    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fmag,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fmag,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
    ay = sf_iaxa(Fmag,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space */

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
	if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay); /* Saved wfld window ny */

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az); /* Saved wfld window oz */
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax); /* Saved wfld window ox */
	if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay); /* Saved wfld window oy */

	if(!sf_getfloat("dqz",&dqz)) dqz=sf_d(az); /* Saved wfld window dz */
	if(!sf_getfloat("dqx",&dqx)) dqx=sf_d(ax); /* Saved wfld window dx */
	if(!sf_getfloat("dqy",&dqy)) dqy=sf_d(ay); /* Saved wfld window dy */

	acz = sf_maxa(nqz,oqz,dqz); if(verb) sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); if(verb) sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); if(verb) sf_raxa(acy);
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

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;
    /*( nb=2 boundary padding in grid points )*/

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

    if(expl) ww = sf_floatalloc( 1);
    else     ww = sf_floatalloc(ns);
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */


    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);
    fdbell3d_init(5);
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
    /*------------------------------------------------------------*/ 
    vel = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    mag = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ele = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    cdt = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

    cdtele =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    magele =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    vt =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /*------------------------------------------------------------*/ 

    /* input magnetic susceptibility*/
    sf_floatread(tt[0][0],nz*nx*ny,Fmag );    expand3d(tt,mag,fdm);
    /* input electric susceptibility*/
    sf_floatread(tt[0][0],nz*nx*ny,Fele );    expand3d(tt,ele,fdm);    
    /* input conductivity*/
    sf_floatread(tt[0][0],nz*nx*ny,Fcdt );    expand3d(tt,cdt,fdm);

    /* cdtele = sigma*dt/2*epsilon */
    /* magele = dt*dt/mu*epsilon */	
		for (iy=0; iy<fdm->nypad; iy++) {
		  for (ix=0; ix<fdm->nxpad; ix++) {
			for (iz=0; iz<fdm->nzpad; iz++) {
                cdtele[iy][ix][iz] = cdt[iy][ix][iz]*dt/(2*(ele[iy][ix][iz]));
	            magele[iy][ix][iz] = dt*dt/(mag[iy][ix][iz]*ele[iy][ix][iz]);
	               vel[iy][ix][iz] = 1./(sqrt(mag[iy][ix][iz]*ele[iy][ix][iz]));
	               vt[iy][ix][iz]  = dt*dt*vel[iy][ix][iz]*vel[iy][ix][iz];
			  }
		    }
		  }		

    if(fsrf) { /* free surface */
	for (iy=0; iy<fdm->nypad; iy++) {
	  for (ix=0; ix<fdm->nxpad; ix++) {
		for (iz=0; iz<fdm->nb; iz++) {
		    cdtele[iy][ix][iz]=0;
	      	magele[iy][ix][iz]=0;
	      	vt[iy][ix][iz]    =0;

		}
	    }
	}
    }

    free(**tt);  free(*tt); free(tt);   
    /*------------------------------------------------------------*/
    if(dabc) {
	/* one-way abc setup */
	abc = abcone3d_make(NOP,dt,vel,fsrf,fdm);
	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }
    /*------------------------------------------------------------*/

    /* allocate wavefield arrays */
    um=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    up=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ua=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    for(iy=0; iy<fdm->nypad; iy++) {
	  for(ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		um[iy][ix][iz]=0;
		uo[iy][ix][iz]=0;
		up[iy][ix][iz]=0;
		ua[iy][ix][iz]=0;
	    }
	   }
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
#pragma omp parallel							\
    private(ix,iy,iz)							\
    shared(fdm,ua,uo,co,cax,cay,caz,cbx,cby,cbz,idx,idy,idz,cdtele,magele)
#endif
	    {
#ifdef _OPENMP
#pragma omp for	schedule(dynamic,fdm->ompchunk)	
#endif

		for (iy=NOP; iy<fdm->nypad-NOP; iy++) {
		  for (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
			for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			    
			    /* 4th order Laplacian operator */
			    ua[iy][ix][iz] = 
				co * uo[iy  ][ix  ][iz  ] + 
				cax*(uo[iy  ][ix-1][iz  ] + uo[iy  ][ix+1][iz  ]) +
				cbx*(uo[iy  ][ix-2][iz  ] + uo[iy  ][ix+2][iz  ]) +
				cay*(uo[iy-1][ix  ][iz  ] + uo[iy+1][ix  ][iz  ]) +
				cby*(uo[iy-2][ix  ][iz  ] + uo[iy+2][ix  ][iz  ]) +
				caz*(uo[iy  ][ix  ][iz-1] + uo[iy  ][ix  ][iz+1]) +
				cbz*(uo[iy  ][ix  ][iz-2] + uo[iy  ][ix  ][iz+2]);
			    
			}
		  }
		}

#ifdef _OPENMP
#pragma omp for	schedule(dynamic,fdm->ompchunk)	
#endif

	/* step forward in time */
		for        (iy=0; iy<fdm->nypad; iy++) {
		  for    (ix=0; ix<fdm->nxpad; ix++) {
			for(iz=0; iz<fdm->nzpad; iz++) {
			    up[iy][ix][iz] =(2*uo[iy][ix][iz] 
				-               (1-cdtele[iy][ix][iz])*um[iy][ix][iz] 
				+                ua[iy][ix][iz]*(magele[iy][ix][iz]))/(1+cdtele[iy][ix][iz]);
			  } 
		     }
		    }
           } /* end parallel section */

    	/* inject displacement source */
	    if(expl) {
		sf_floatread(ww, 1,Fwav);
		lint3d_bell1(up,ww[0],cs);
	    } else {
		sf_floatread(ww,ns,Fwav);	
		lint3d_bell(up,ww,cs);
	    }
	    
	    /* extract data at receivers */
	    lint3d_extract(up,dd,cr);
	    if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);
	    
	    /* extract wavefield in the "box" */
	    if(snap && it%jsnap==0) {
		cut3d(up,uc,fdm,acz,acx,acy);
		sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
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
	    }

	} /* end time loop */
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(**um); free(*um); free(um);
    free(**up); free(*up); free(up);
    free(**uo); free(*uo); free(uo);
    free(**ua); free(*ua); free(ua);
    if (snap) {
       free(**uc); free(*uc); free(uc);
    }
        
    free(ww);
    free(ss);
    free(rr);
    free(dd);

    free(**mag);    free(*mag);    free(mag);
    free(**ele);    free(*ele);    free(ele);
    free(**cdt);    free(*cdt);    free(cdt);
    free(**vel);    free(*vel);    free(vel);  
    free(**cdtele); free(*cdtele); free(cdtele);
    free(**magele); free(*magele); free(magele);
    free(**vt);     free(*vt);     free(vt);
    
    /*------------------------------------------------------------*/
    /* CLOSE FILES AND EXIT */
    if (Fwav!=NULL) sf_fileclose(Fwav);  
    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Frec!=NULL) sf_fileclose(Frec);
    
    if (Fmag!=NULL) sf_fileclose(Fmag);
    if (Fele!=NULL) sf_fileclose(Fele);
    if (Fcdt!=NULL) sf_fileclose(Fcdt);

    if (Fdat!=NULL) sf_fileclose(Fdat); 
    if (Fwfl!=NULL) sf_fileclose(Fwfl);

    exit (0);
}
