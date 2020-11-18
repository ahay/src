/* Heat diffusion equation FD modeling */
/*
  Copyright (C) 2010 Colorado School of Mines
  
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

#define NOP 2 /* derivative operator half-size */

#define C0 -2.500000 /*    c0=-30./12.; */
#define CA +1.333333 /*    ca=+16./12.; */
#define CB -0.083333 /*    cb=- 1./12.; */

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap; 
    int  jsnap,ntsnap,jdata;
#ifdef _OPENMP
    int ompnth=1;
#endif

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fcon=NULL; /* conductivity  */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,az,ax,as,ar;
    int     nt,nz,nx,ns,nr,nb;
    int     it,iz,ix;
    float   dt,dz,dx,idz,idx;

    /* FDM structures */
    fdm2d    fdm;

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    float  *dd=NULL;           /* data      */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */

    float **kkin=NULL;         /* conductivity  */
    float **kk=NULL;           /* conductivity in expanded domain */

    float **uo,**up,**ua,**ut; /* wavefield: uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float co,cax,cbx,caz,cbz;

    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
    if(!ompnth)
        abort();
#endif
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag*/

    Fwav = sf_input ("in" ); /* wavelet   */
    Fcon = sf_input ("con"); /* conductivity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    az = sf_iaxa(Fcon,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fcon,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    nt = sf_n(at); dt = sf_d(at);
    ns = sf_n(as);
    nr = sf_n(ar);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

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
    /* expand domain for FD operators */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
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

	acz = sf_maxa(nqz,oqz,dqz);
	acx = sf_maxa(nqx,oqx,dqx);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,at,3);
    }

    ww  =sf_floatalloc(ns);
    dd  =sf_floatalloc(nr);

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
    caz= CA *          idz*idz;
    cbz= CB *          idz*idz;

    /*------------------------------------------------------------*/ 
    /* setup model */
    kkin=sf_floatalloc2(nz,   nx   ); 
    kk  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(kkin[0],nz*nx,Fcon);
    expand(kkin,kk,fdm);
    free(*kkin); free(kkin);

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    up=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
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
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b%d",it);

	/* inject heat source */
	sf_floatread(ww,ns,Fwav);	
	lint2d_inject(uo,ww,cs);
	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(ix,iz) shared(fdm,ua,uo,co,cax,caz,cbx,cbz,idx,idz)
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
#pragma omp parallel for schedule(dynamic,fdm->ompchunk) private(ix,iz) shared(fdm,ua,uo,up,kk,dt)
#endif
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		up[ix][iz] = uo[ix][iz] 
		    +        ua[ix][iz] * kk[ix][iz] * dt;
	    }
	}
	/* circulate wavefield arrays */
	ut=uo;
	uo=up;
	up=ut;

	/* re-fill boundaries */
	bfill(uo,fdm);

	/* extract data */
	lint2d_extract(uo,dd,cr);
	if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

	/* extract wavefield in the "box" */
	if(snap && it%jsnap==0) {
	    cut2d(uo,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
    }
    if(verb) fprintf(stderr,"\n");    

/* deallocate arrays */
    free(*up); free(up);
    free(*uo); free(uo);
    free(*ua); free(ua);
    if(snap) {
        free(*uc); free(uc);
    }

    free(ww);
    free(ss);
    free(rr);
    free(dd);

    exit (0);
}
