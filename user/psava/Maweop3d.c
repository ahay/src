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

    /* I/O files */
    sf_file Fvel=NULL; /* velocity */
    sf_file Fm=NULL;   /* model    */
    sf_file Fd=NULL;   /* data     */

    /* cube axes */
    sf_axis at,az,ay,ax;
    int     nt,nz,ny,nx,nb;
    int     it,iz,iy,ix;
    float   dt,dz,dy,dx,idz,idy,idx;

    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float ***tt=NULL;
    float ***vp=NULL;           /* velocity */
    float ***aa=NULL;           /* source   */
    float ***um,***uo,***up,***ua,***ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* FD operator */
    float co,cax,cbx,cay,cby,caz,cbz;

    int nslice;  /* wavefield slice size */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity  */
    if(! sf_getbool("fsrf",&fsrf)) fsrf=false; /* free surface  */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
    if(! sf_getbool("adj" ,&adj))   adj=false; /* adjoint flag */

    /*------------------------------------------------------------*/
    /* I/O files */
    Fvel = sf_input ("vel"); /* velocity  */

    if(adj) {
	Fd = sf_input ("in" ); /*  data */
	Fm = sf_output("out"); /* model */
    } else {
	Fm = sf_input ("in" ); /* model */
	Fd = sf_output("out"); /*  data */
    }

    /*------------------------------------------------------------*/
    /* axes */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* y */

    if(adj) at = sf_iaxa(Fd,4); 
    else    at = sf_iaxa(Fm,4);
    sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* t */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    nslice = nz*nx*ny*sizeof(float); /* wavefield slice */
    /*------------------------------------------------------------*/

    if(adj) {
	sf_oaxa(Fm,az,1);
	sf_oaxa(Fm,ax,2);
	sf_oaxa(Fm,ay,3);
	sf_oaxa(Fm,at,4);
    } else {
	sf_oaxa(Fd,az,1);
	sf_oaxa(Fd,ax,2);
	sf_oaxa(Fd,ay,3);
	sf_oaxa(Fd,at,4);
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
    aa=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    um=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    up=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ua=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    PADLOOP( aa[iy][ix][iz]=um[iy][ix][iz]=uo[iy][ix][iz]=up[iy][ix][iz]=ua[iy][ix][iz]=0; );

    WFLLOOP( tt[iy][ix][iz]=0; );
    if(adj) {
	for(it=0;it<nt;it++) sf_floatwrite(tt[0],nz*nx,Fm); /* reserve wfl */ 
	sf_seek(Fm,0,SEEK_SET);                             /* seek back */
    } else {
	for(it=0;it<nt;it++) sf_floatwrite(tt[0],nz*nx,Fd); /* reserve wfl */ 
	sf_seek(Fd,0,SEEK_SET);                             /* seek back */
    }

    /* 
     * MAIN LOOP
     */
    for (it=0;it<nt;it++){ if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	if(adj) {
	    sf_seek(Fd,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatread(tt[0][0],nz*nx*ny,Fd);
	} else {
	    sf_floatread(tt[0][0],nz*nx*ny,Fm);
	}
	wpad3d(aa,tt,fdm);

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
	PADLOOP( up[iy][ix][iz] = 2*uo[iy][ix][iz] - um[iy][ix][iz] + (aa[iy][ix][iz]+ua[iy][ix][iz]) * vp[iy][ix][iz]; );  

	wwin3d(tt,up,fdm);
	if(adj) {
	    sf_seek(Fm,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatwrite(tt[0][0],nz*nx*ny,Fm); 
	} else {
	    sf_floatwrite(tt[0][0],nz*nx*ny,Fd); 
	}

	/* circulate wavefield arrays */
	ut=um; um=uo; uo=up; up=ut; 
	/* one-way abc */
	if(dabc) abcone3d_apply(uo,um,NOP,abc,fdm);
    } /* it */
    if(verb) fprintf(stderr,"\n");

    /*------------------------------------------------------------*/
    free(**aa); free(*aa); free(aa); /* source */
    free(**um); free(*um); free(um); /* wavefield */
    free(**up); free(*up); free(up);
    free(**uo); free(*uo); free(uo);
    free(**ua); free(*ua); free(ua); /* acceleration */
    free(**tt); free(*tt); free(tt); /* temporary */
    free(**vp); free(*vp); free(vp); /* velocity */
    
    exit (0);
}
