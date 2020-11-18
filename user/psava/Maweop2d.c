/* 2D AWE modeling */
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
    for    (ix=0; ix<fdm->nx; ix++) {		\
	for(iz=0; iz<fdm->nz; iz++) {		\
	    {a}					\
	}					\
    }

#define PADLOOP(a)				\
    for    (ix=0; ix<fdm->nxpad; ix++) {	\
	for(iz=0; iz<fdm->nzpad; iz++) {	\
	    {a}					\
	}					\
    }

#define FDMLOOP(a)				\
    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {	\
	for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {	\
	    {a}					\
	}					\
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
    sf_axis at,az,ax,as,ar,aa;
    int     nt,nz,nx,nb;
    int     it,iz,ix,is,ir;
    float   dt,dz,dx,idz,idx;
    
    pt2d   *ss=NULL,*rr=NULL; /* source/receiver coordinates */
    lint2d  cs=NULL, cr=NULL; /* weights/indices */
    
    /* FDM structure */
    fdm2d    fdm=NULL;
    // abcone2d abc=NULL;
    sponge   spo=NULL;
    
    /* I/O arrays */
    float **tt=NULL;
    float **vp=NULL;                    /* velocity */
    float **ww=NULL, *ws=NULL,*wr=NULL; /* inject/extract data */
    float **um,**uo,**up,**ua,**ut;     /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    
    /* FD operator */
    float co,cax,cbx,caz,cbz;
    
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
        ss = (pt2d*) sf_alloc(sf_n(as),sizeof(*ss));
        pt2dread1(Fs,ss,sf_n(as),2);
    }
    if(verb) sf_warning("spoint=%d",spoint);

    rpoint=false;
    if(NULL != sf_getstring("rec")){ rpoint=true; Fr = sf_input("rec"); }
    if(rpoint) {
        ar = sf_iaxa(Fr,2); sf_setlabel(ar,"r"); sf_setunit(ar,"");
        rr = (pt2d*) sf_alloc(sf_n(ar),sizeof(*rr));
        pt2dread1(Fr,rr,sf_n(ar),2);
    }    
    if(verb) sf_warning("rpoint=%d",rpoint);

    /*------------------------------------------------------------*/
    /* I/O files */
    Fvel = sf_input ("vel"); /* velocity  */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    
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
        else       at = sf_iaxa(Fd,3);
    } else {
        if(spoint) at = sf_iaxa(Fm,2);
        else       at = sf_iaxa(Fm,3);
    }
    sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* t */
    
    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

    /*------------------------------------------------------------*/
    
    if(adj) {
        if(spoint) {
            sf_oaxa(Fm,as,1);
            sf_oaxa(Fm,at,2);
            sf_oaxa(Fm,aa,3);
        } else {
            sf_oaxa(Fm,az,1);
            sf_oaxa(Fm,ax,2);
            sf_oaxa(Fm,at,3);
        }
    } else {
        if(rpoint) {
            sf_oaxa(Fd,ar,1);
            sf_oaxa(Fd,at,2);
            sf_oaxa(Fd,aa,3);
        } else {
            sf_oaxa(Fd,az,1);
            sf_oaxa(Fd,ax,2);
            sf_oaxa(Fd,at,3);
        }
    }
    
    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;
    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);
    
    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/
    
    /* interpolation coefficients */
    if(spoint) cs = lint2d_make(sf_n(as),ss,fdm);
    if(rpoint) cr = lint2d_make(sf_n(ar),rr,fdm);
    
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
    vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    
    /* input velocity */
    sf_floatread(tt[0],nz*nx,Fvel); expand(tt,vp,fdm);
    
    if(dabc) {
        // abc = abcone2d_make(NOP,dt,vp,fsrf,fdm); /* one-way abc setup */
        spo = sponge_make(fdm->nb);              /* sponge abc setup */
    }
    
    /* precompute vp^2 * dt^2 */
    PADLOOP( vp[ix][iz] = vp[ix][iz]*vp[ix][iz] * dt*dt; );
    if(fsrf) { /* free surface */
        for    (ix=0; ix<fdm->nxpad; ix++) {
            for(iz=0; iz<fdm->nb; iz++) {
                vp[ix][iz]=0;
            }
        }
    }
    
    /*------------------------------------------------------------*/
    /* allocate arrays */
    if(spoint) {ws=sf_floatalloc(sf_n(as)); for(is=0; is<sf_n(as); is++) ws[is]=0;}
    if(rpoint) {wr=sf_floatalloc(sf_n(ar)); for(ir=0; ir<sf_n(ar); ir++) wr[ir]=0;}    
    ww=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    um=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    up=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    PADLOOP( ww[ix][iz]=um[ix][iz]=uo[ix][iz]=up[ix][iz]=ua[ix][iz]=0; );
    
    WFLLOOP( tt[ix][iz]=0; );
    if(adj) { /* reserve output disk space */
        if(spoint) for(it=0;it<nt;it++) sf_floatwrite(ws,sf_n(as),Fm);
        else       for(it=0;it<nt;it++) sf_floatwrite(tt[0],nz*nx,Fm);
        sf_seek(Fm,0,SEEK_SET);
    } else {
        if(rpoint) for(it=0;it<nt;it++) sf_floatwrite(wr,sf_n(ar),Fd);
        else       for(it=0;it<nt;it++) sf_floatwrite(tt[0],nz*nx,Fd);
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
                PADLOOP( ww[ix][iz]=0;);
                lint2d_inject(ww,wr,cr);
            } else {
                sf_seek(Fd,(off_t)(nt-1-it)*nz*nx*sizeof(float),SEEK_SET);
                sf_floatread(tt[0],nz*nx,Fd);
		wpad2d(ww,tt,fdm);
            }
        } else {
            if(spoint) {
                sf_floatread(ws,sf_n(as),Fm);
                PADLOOP( ww[ix][iz]=0;);
		lint2d_inject(ww,ws,cs);
            }
            else {
		sf_floatread(tt[0],nz*nx,Fm);
		wpad2d(ww,tt,fdm);
	    }
        }
        
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ix,iz)				\
    shared(fdm,ua,uo,co,cax,cbx,caz,cbz)
#endif
        
        FDMLOOP( /* Laplacian */
	    ua[ix][iz] =     co * uo[ix  ][iz  ] +
	    cax*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
	    cbx*(uo[ix-2][iz  ] + uo[ix+2][iz  ]) +
	    caz*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
	    cbz*(uo[ix  ][iz-2] + uo[ix  ][iz+2]);
	    );
        
        /* sponge abc */
        if(dabc) sponge2d_apply(ua,spo,fdm);
        
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ix,iz)				\
    shared(fdm,ua,uo,up,um,ww,vp)
#endif
        /* time step */
        PADLOOP( up[ix][iz] = 2*uo[ix][iz] - um[ix][iz] + (ww[ix][iz]+ua[ix][iz]) * vp[ix][iz]; );
        
	/* write data/wavefield */
        if(adj) {
            if(spoint) {
                lint2d_extract(up,ws,cs);
                sf_seek(Fm,(off_t)(nt-1-it)*sf_n(as)*sizeof(float),SEEK_SET);
                sf_floatwrite(ws,sf_n(as),Fm);
            } else {
		wwin2d(tt,up,fdm);
                sf_seek(Fm,(off_t)(nt-1-it)*nz*nx*sizeof(float),SEEK_SET);
                sf_floatwrite(tt[0],nz*nx,Fm);
            }
        } else {
            if(rpoint) {
                lint2d_extract(up,wr,cr);
                sf_floatwrite(wr,sf_n(ar),Fd);
            }
            else {
		wwin2d(tt,up,fdm);
		sf_floatwrite(tt[0],nz*nx,Fd);
	    }
        }
        
        /* circulate wavefield arrays */
        ut=um; um=uo; uo=up; up=ut;
        /* one-way abc */
        /* if(dabc) abcone2d_apply(uo,um,NOP,abc,fdm); */
    }
    if(verb) fprintf(stderr,"\n");
    
    /*------------------------------------------------------------*/
    free(*ww); free(ww); /* inject/extract data */
    free(ws);
    free(wr);

    free(*um); free(um); /* wavefield */
    free(*up); free(up);
    free(*uo); free(uo);
    free(*ua); free(ua); /* acceleration */

    free(*tt); free(tt); /* temporary */
    free(*vp); free(vp); /* velocity */
    
    exit (0);
}

