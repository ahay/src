/* 3D AWE modeling
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

#define WFLLOOP(a)				\
    for        (iy=0; iy<fdm->ny; iy++) {	\
	for    (ix=0; ix<fdm->nx; ix++) {	\
	    for(iz=0; iz<fdm->nz; iz++) {	\
		{a}				\
	    }					\
	}					\
    }

#define PADLOOP(a)				\
    for        (iy=0; iy<fdm->nypad; iy++) {	\
	for    (ix=0; ix<fdm->nxpad; ix++) {	\
	    for(iz=0; iz<fdm->nzpad; iz++) {	\
		{a}				\
	    }					\
	}					\
    }

#define FDMLOOP(a)						   \
    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {		   \
        for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {		   \
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		   \
		{a}						   \
	    }							   \
	}							   \
    }

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool verb,fsrf,dabc; 
    bool adj;            /* adjoint operator flag */
    bool spoint, rpoint; /* indicates point sources/receivers */

    /* I/O files */
    sf_file Fvel=NULL; /* velocity */
    sf_file Fm=NULL;   /* model    */
    sf_file Fd=NULL;   /* data     */
    sf_file Fs=NULL;   /* sources   */
    sf_file Fr=NULL;   /* receivers */

    /* cube axes */
    sf_axis at,az,ay,ax;
    sf_axis as,ar,aa;
    int     nt,nz,ny,nx,nb;
    int     it,iz,iy,ix,is,ir;
    float   dt,dz,dy,dx,idz,idy,idx;
    
    pt3d   *ss=NULL,*rr=NULL; /* source/receiver coordinates */
    lint3d  cs=NULL, cr=NULL; /* weights/indices */

    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float ***tt=NULL;
    float ***vp=NULL;           /* velocity */
    float ***ww=NULL;           /* padded inject/extract wavefield   */
    float ***um,***uo,***up,***ua,***ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float *ws=NULL,*wr=NULL; /* inject/extract data */
    

    /* FD operator */
    float co,cax,cbx,cay,cby,caz,cbz;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity  */
    if(! sf_getbool("fsrf",&fsrf)) fsrf=false; /* free surface  */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
    if(! sf_getbool("adj" ,&adj))   adj=false; /* adjoint flag */
    
    aa = sf_maxa(1,0,1);
    
    /*------------------------------------------------------------*/
    /* check if input/output is at points */
    
    spoint=false;
    if(NULL != sf_getstring("sou")){ spoint=true; Fs = sf_input("sou"); }
    if(spoint) {
        as = sf_iaxa(Fs,2); sf_setlabel(as,"s"); sf_setunit(as,"");
        ss = (pt3d*) sf_alloc(sf_n(as),sizeof(*ss));
        pt3dread1(Fs,ss,sf_n(as),3);
    }
    if(verb) sf_warning("spoint=%d",spoint);
    
    rpoint=false;
    if(NULL != sf_getstring("rec")){ rpoint=true; Fr = sf_input("rec"); }
    if(rpoint) {
        ar = sf_iaxa(Fr,2); sf_setlabel(ar,"r"); sf_setunit(ar,"");
        rr = (pt3d*) sf_alloc(sf_n(ar),sizeof(*rr));
        pt3dread1(Fr,rr,sf_n(ar),3);
    }
    if(verb) sf_warning("rpoint=%d",rpoint);

    /*------------------------------------------------------------*/
    /* I/O files */
    Fvel = sf_input ("vel"); /* velocity  */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* y */

    if(adj) {
	Fd = sf_input ("in" ); /*  data */
	Fm = sf_output("out"); /* model */
    } else {
	Fm = sf_input ("in" ); /* model */
	Fd = sf_output("out"); /*  data */
    }
    
    /*------------------------------------------------------------*/
    
    if(adj) {
        if(rpoint) at = sf_iaxa(Fd,2);
        else       at = sf_iaxa(Fd,4);
    } else {
        if(spoint) at = sf_iaxa(Fm,2);
        else       at = sf_iaxa(Fm,4);
    }
    sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* t */
    
    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    /*------------------------------------------------------------*/
    
    if(adj) {
        if(spoint) {
            sf_oaxa(Fm,as,1);
            sf_oaxa(Fm,at,2);
            sf_oaxa(Fm,aa,3);
	    sf_oaxa(Fm,aa,4);
        } else {
            sf_oaxa(Fm,az,1);
            sf_oaxa(Fm,ax,2);
            sf_oaxa(Fm,ay,3);
            sf_oaxa(Fm,at,4);
        }
    } else {
        if(rpoint) {
            sf_oaxa(Fd,ar,1);
            sf_oaxa(Fd,at,2);
            sf_oaxa(Fd,aa,3);
	    sf_oaxa(Fd,aa,4);
        } else {
            sf_oaxa(Fd,az,1);
            sf_oaxa(Fd,ax,2);
            sf_oaxa(Fd,ay,3);
            sf_oaxa(Fd,at,4);
        }
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;
    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/
    /* interpolation coefficients */
    if(spoint) cs = lint3d_make(sf_n(as),ss,fdm);
    if(rpoint) cr = lint3d_make(sf_n(ar),rr,fdm);

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
    vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

    /* input velocity */
    sf_floatread(tt[0][0],nz*nx*ny,Fvel );    expand3d(tt,vp,fdm);

    if(dabc) {
        abc = abcone3d_make(NOP,dt,vp,fsrf,fdm); /* one-way abc setup */
        spo = sponge_make(fdm->nb);              /* sponge abc setup */
    }

    /* precompute vp^2 * dt^2 */
    PADLOOP( vp[iy][ix][iz] = vp[iy][ix][iz]*vp[iy][ix][iz] * dt*dt; );
    if(fsrf) { /* free surface */
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb;    iz++) {
		    vp[iy][ix][iz]=0;
		}
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* allocate arrays */
    if(spoint) {ws=sf_floatalloc(sf_n(as)); for(is=0; is<sf_n(as); is++) ws[is]=0;}
    if(rpoint) {wr=sf_floatalloc(sf_n(ar)); for(ir=0; ir<sf_n(ar); ir++) wr[ir]=0;}
    ww=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    
    um=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    up=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ua=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    PADLOOP( ww[iy][ix][iz]=um[iy][ix][iz]=uo[iy][ix][iz]=up[iy][ix][iz]=ua[iy][ix][iz]=0; );

    WFLLOOP( tt[iy][ix][iz]=0; );
    if(adj) {   /* reserve output disk space */
        if(spoint) for(it=0;it<nt;it++) sf_floatwrite(ws,sf_n(as),Fm);
        else for(it=0;it<nt;it++) sf_floatwrite(tt[0][0],nz*nx*ny,Fm);
        sf_seek(Fm,0,SEEK_SET);
    } else {
        if(rpoint)  for(it=0;it<nt;it++) sf_floatwrite(wr,sf_n(ar),Fd);
        else        for(it=0;it<nt;it++) sf_floatwrite(tt[0][0],nz*nx*ny,Fd);
        sf_seek(Fd,0,SEEK_SET);
    }

    /* 
     * MAIN LOOP
     */
    for (it=0;it<nt;it++){ if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
        
    /* read data/wavefield */
        if(adj) {
            if(rpoint) {
                sf_seek(Fd,(off_t)(nt-1-it)*sf_n(ar)*sizeof(float),SEEK_SET);
                sf_floatread(wr,sf_n(ar),Fd);
                PADLOOP( ww[iy][ix][iz]=0;);
                lint3d_inject(ww,wr,cr);
            } else {
                sf_seek(Fd,(off_t)(nt-1-it)*nz*nx*ny*sizeof(float),SEEK_SET);
                sf_floatread(tt[0][0],nz*nx*ny,Fd);
                wpad3d(ww,tt,fdm);
            }
        } else {
            if(spoint) {
                sf_floatread(ws,sf_n(as),Fm);
                PADLOOP( ww[iy][ix][iz]=0;);
                lint3d_inject(ww,ws,cs);
            }
            else {
                sf_floatread(tt[0][0],nz*nx*ny,Fm);
                wpad3d(ww,tt,fdm);
            }
        }
        
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic)						\
    private(ix,iy,iz)						\
    shared(fdm,ua,uo,um,up,aa,vp,co,cax,cbx,cay,cby,caz,cbz)
#endif
        FDMLOOP( /*  Laplacian */
	    ua[iy][ix][iz] =       co * uo[iy  ][ix  ][iz  ] + 
	    cax*(uo[iy  ][ix-1][iz  ] + uo[iy  ][ix+1][iz  ]) +
	    cbx*(uo[iy  ][ix-2][iz  ] + uo[iy  ][ix+2][iz  ]) +
	    cay*(uo[iy-1][ix  ][iz  ] + uo[iy+1][ix  ][iz  ]) +
	    cby*(uo[iy-2][ix  ][iz  ] + uo[iy+2][ix  ][iz  ]) +
	    caz*(uo[iy  ][ix  ][iz-1] + uo[iy  ][ix  ][iz+1]) +
	    cbz*(uo[iy  ][ix  ][iz-2] + uo[iy  ][ix  ][iz+2]);
	    );
	
	/* sponge abc */
        if(dabc) sponge3d_apply(ua,spo,fdm);

#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic)							\
    private(ix,iz)							\
    shared(fdm,ua,uo,up,um,aa,vp)
#endif
	/* time step */
        PADLOOP( up[iy][ix][iz] = 2*uo[iy][ix][iz] - um[iy][ix][iz] + (ww[iy][ix][iz]+ua[iy][ix][iz]) * vp[iy][ix][iz]; );

    /* write data/wavefield */
        if(adj) {
            if(spoint) {
                lint3d_extract(up,ws,cs);
                sf_seek(Fm,(off_t)(nt-1-it)*sf_n(as)*sizeof(float),SEEK_SET);
                sf_floatwrite(ws,sf_n(as),Fm);
            } else {
                wwin3d(tt,up,fdm);
                sf_seek(Fm,(off_t)(nt-1-it)*nz*nx*ny*sizeof(float),SEEK_SET);
                sf_floatwrite(tt[0][0],nz*nx*ny,Fm);
            }
        } else {
            if(rpoint) {
                lint3d_extract(up,wr,cr);
                sf_floatwrite(wr,sf_n(ar),Fd);
            }
            else {
                wwin3d(tt,up,fdm);
                sf_floatwrite(tt[0][0],nz*nx*ny,Fd);
            }
        }

        /* circulate wavefield arrays */
        ut=um; um=uo; uo=up; up=ut;
        /* one-way abc */
        /* if(dabc) abcone3d_apply(uo,um,NOP,abc,fdm); */
    } /* it */
    if(verb) fprintf(stderr,"\n");

    /*------------------------------------------------------------*/
    free(**ww); free(*ww); free(ww); /* inject/extract data */
    free(ws);
    free(wr);
    
    free(**um); free(*um); free(um); /* wavefield */
    free(**up); free(*up); free(up);
    free(**uo); free(*uo); free(uo);
    free(**ua); free(*ua); free(ua); /* acceleration */
    
    free(**tt); free(*tt); free(tt); /* temporary */
    free(**vp); free(*vp); free(vp); /* velocity */
    
    exit (0);
}
