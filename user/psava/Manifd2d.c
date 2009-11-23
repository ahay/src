/* 2D anisotropic time-domain FD modeling */
/*
  Copyright (C) 2009 Colorado School of Mines
  
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

#define Dxx(a,ix,iz,co,ca,cb) (co* a[ix  ][iz  ] +			\
			       ca*(a[ix-1][iz  ] + a[ix+1][iz  ]) +	\
			       cb*(a[ix-2][iz  ] + a[ix+2][iz  ]) )

#define Dzz(a,ix,iz,co,ca,cb) (co* a[ix  ][iz  ] +			\
			       ca*(a[ix  ][iz-1] + a[ix  ][iz+1]) +	\
			       cb*(a[ix  ][iz-2] + a[ix  ][iz+2]) )

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

#define Dxz(a,ix,iz,c1x,c2x,c1z,c2z) \
    (c2x * D1z(a,ix+2,iz,c1z,c2z) -  \
     c2x * D1z(a,ix-2,iz,c1z,c2z) +  \
     c1x * D1z(a,ix+1,iz,c1z,c2z) -  \
     c1x * D1z(a,ix-1,iz,c1z,c2z))

#define Dzx(a,ix,iz,c1x,c2x,c1z,c2z) \
    (c2z * D1x(a,ix,iz+2,c1x,c2x) -  \
     c2z * D1x(a,ix,iz-2,c1x,c2x) +  \
     c1z * D1x(a,ix,iz+1,c1x,c2x) -  \
     c1z * D1x(a,ix,iz-1,c1x,c2x))

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc; 
    int  jsnap,ntsnap,jdata;
    char *atype;

    /* OMP parameters */
#ifdef _OPENMP
    int ompnth;
#endif 

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fvel=NULL; /* velocity  */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,az,ax;
    sf_axis as,ar;

    int     nt,nz,nx,ns,nr,nb;
    int     it,iz,ix;
    float   dt,dz,dx,idz,idx,dt2;

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

    float **vn2=NULL;
    float **vv2=NULL;
    float **vh2=NULL;

    float **ang=NULL;
    float **sia=NULL;
    float **coa=NULL;

    float **rm=NULL,**ro=NULL,**rp=NULL,**ra=NULL,**rt; /*      main wavefield */
    float **qm,**qo,**qp,**qa,**qt; /* auxiliary wavefield */

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float cox,coz,cax,cbx,caz,cbz;
    float c1x,c2x,c1z,c2z;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **qc=NULL;

    float H2q,H1r;
    float s2,c2,sc;

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
	    sf_warning("ISO model");
	    break;

	default:
	    sf_warning("ISO modeling - by default!");
	    break;
    }

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif
    /*------------------------------------------------------------*/

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

	qc=sf_floatalloc2(sf_n(acz),sf_n(acx));

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

    cox= C0 * (idx*idx);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;

    coz= C0 * (idz*idz);
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    c1x = C1 * idx;
    c2x = C2 * idx;

    c1z = C1 * idz;
    c2z = C2 * idz;

    /* precompute dt^2*/
    dt2 = dt*dt;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 
    /*------------------------------------------------------------*/

    /* input velocity */
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

    vv2 =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vv2,fdm); /* vertical v */

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {	    
	    vp [ix][iz] = vv2[ix][iz];
	    vv2[ix][iz] = vv2[ix][iz] * vv2[ix][iz];
	}
    }
    if(fsrf) { /* free surface */
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nb; iz++) {
		vv2[ix][iz]=0;
	    }
	}
    }

    if(atype[0] != 'i') {
	vn2 =sf_floatalloc2(fdm->nzpad,fdm->nxpad);     
	vh2 =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	
	sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vn2,fdm); /* NMO v */
	sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vh2,fdm); /* horizontal v */
	
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {	    
		vn2[ix][iz] = vn2[ix][iz] * vn2[ix][iz];
		vh2[ix][iz] = vh2[ix][iz] * vh2[ix][iz];
	    }
	}

	if(fsrf) { /* free surface */
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb; iz++) {
		    vn2[ix][iz]=0;
		    vh2[ix][iz]=0;
		}
	    }
	}
    }

    /*------------------------------------------------------------*/

    if( atype[0]=='t') {
	/* input tilt angle */
	ang =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	sia =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	coa =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	
	sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,ang,fdm);
	
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {	    
		ang[ix][iz] *= SF_PI/180.;
		sia[ix][iz] = sinf(ang[ix][iz]);
		coa[ix][iz] = cosf(ang[ix][iz]);
	    }
	}
    }

    /*------------------------------------------------------------*/
    free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
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

    if(atype[0] != 'i') {
	rm=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	ro=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	rp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	ra=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		rm[ix][iz]=0;
		ro[ix][iz]=0;
		rp[ix][iz]=0;
		ra[ix][iz]=0;
	    }
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

	/* compute acceleration */
	switch(atype[0]) {
	    case 't':

#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iz,H2q,H1r,s2,c2,sc)					\
    shared(fdm,ra,ro,qa,qo,cox,coz,cax,caz,cbx,cbz,c1x,c2x,c1z,c2z,idx,idz,vn2,vv2,vh2)
#endif
		for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			
			s2 =  sia[ix][iz]*sia[ix][iz];
			c2 =  coa[ix][iz]*coa[ix][iz];
			sc =2*sia[ix][iz]*coa[ix][iz];
			
			H2q = c2*Dxx(qo,ix,iz,cox,cax,cbx) 
			    + s2*Dzz(qo,ix,iz,coz,caz,cbz)
			    - sc*Dxz(qo,ix,iz,c1x,c2x,c1z,c2z);
			
			H1r = s2*Dxx(ro,ix,iz,cox,cax,cbx)
			    + c2*Dzz(ro,ix,iz,coz,caz,cbz)
			    + sc*Dxz(ro,ix,iz,c1x,c2x,c1z,c2z);
			
			/* main field - q */
			qa[ix][iz] = H2q * vh2[ix][iz] + H1r * vv2[ix][iz] ;			
			/* auxiliary field - r */
			ra[ix][iz] = H2q * vn2[ix][iz] + H1r * vv2[ix][iz] ;
			
		    }
		}   
		break;
		
	    case 'v':

#ifdef _OPENMP
#pragma omp parallel for				\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iz,H2q,H1r)					\
    shared(fdm,ra,ro,qa,qo,cox,coz,cax,caz,cbx,cbz,idx,idz,vn2)
#endif
		for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			
			H2q = Dxx(qo,ix,iz,cox,cax,cbx);
			
			H1r = Dzz(ro,ix,iz,coz,caz,cbz);
			
			/* main field - q */
			qa[ix][iz] = H2q * vh2[ix][iz] + H1r * vv2[ix][iz] ;			
			/* auxiliary field - r */
			ra[ix][iz] = H2q * vn2[ix][iz] + H1r * vv2[ix][iz] ;
			
		    }
		}   
		break;

	    case 'i':
	    default:
	
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iz)						\
    shared(fdm,qa,qo,cox,coz,cax,caz,cbx,cbz,idx,idz,vn2)
#endif
		for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			
			qa[ix][iz] = ( Dxx(qo,ix,iz,cox,cax,cbx) + 
				       Dzz(qo,ix,iz,coz,caz,cbz) ) * vv2[ix][iz];
			
		    }
		}   
		break;
	}

	/* inject acceleration source */
	if(expl) {
	    sf_floatread(ww, 1,Fwav);
	    lint2d_inject1(qa,ww[0],cs);
	    if(atype[0] != 'i') lint2d_inject1(ra,ww[0],cs);
	} else {
	    sf_floatread(ww,ns,Fwav);	
	    lint2d_inject(qa,ww,cs);
	    if(atype[0] != 'i') lint2d_inject(ra,ww,cs);
	}

	/* step forward in time */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk) \
    private(ix,iz)		    \
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
	
	if(atype[0] != 'i') {
	    
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(ix,iz)				\
    shared(fdm,ra,ro,rm,rp,dt2)
#endif
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    rp[ix][iz] = 2*ro[ix][iz] 
			-          rm[ix][iz] 
			+          ra[ix][iz] * dt2;
		}
	    }
	    /* circulate wavefield arrays */
	    rt=rm;
	    rm=ro;
	    ro=rp;
	    rp=rt;
	}

	if(dabc) {
	    /* one-way abc apply */
	    abcone2d_apply(qo,qm,NOP,abc,fdm);
	    sponge2d_apply(qm,spo,fdm);
	    sponge2d_apply(qo,spo,fdm);
	    sponge2d_apply(qp,spo,fdm);

	    if(atype[0] != 'i') {
		/* one-way abc apply */
		abcone2d_apply(ro,rm,NOP,abc,fdm);
		sponge2d_apply(rm,spo,fdm);
		sponge2d_apply(ro,spo,fdm);
		sponge2d_apply(rp,spo,fdm);
	    }
	}

	/* extract data */
	lint2d_extract(qo,dd,cr);

	if(snap && it%jsnap==0) {
	    cut2d(qo,qc,fdm,acz,acx);
	    sf_floatwrite(qc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
	if(        it%jdata==0) 
	    sf_floatwrite(dd,nr,Fdat);
    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */

    if(atype[0] != 'i') {
	free(*rm); free(rm);
	free(*rp); free(rp);
	free(*ro); free(ro);
	free(*ra); free(ra);
    }

    free(*qm); free(qm);
    free(*qp); free(qp);
    free(*qo); free(qo);
    free(*qa); free(qa);
    free(*qc); free(qc);

    free(*vp);  free(vp);
    free(*vv2); free(vv2);
    if(atype[0] != 'i') {
	free(*vn2); free(vn2);
	free(*vh2); free(vh2);
    }

    if(atype[0] == 't') {
	free(*ang); free(ang);
	free(*sia); free(sia);
	free(*coa); free(coa);
    }

    free(ww);
    free(ss);
    free(rr);
    free(dd);
    /*------------------------------------------------------------*/


    exit (0);
}

