/* 2-D finite-difference zero-offset exploding reflector modeling/migration */

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*******************************************************/
/* wave propagation struct */
typedef struct Par {
    /*geometry parameters*/
    int nx,nz,gpz;
    float dx, dz;
    /*time parameters*/
    int nt,snap,wfnt;
    float dt;
    /*ABC*/
    int nb,nxb,nzb;
    float c;
    /*verbosity*/
    bool verb;
} * par; /*geometry parameters*/

int fdexp(float **img, float **dat, bool adj, par pars, float **vv, float ***wvfld)
/*< zero-offset exploding reflector modeling/migration >*/
{
    int nx, nz, gpz;
    float dx, dz;
    int nt, snap, wfnt;
    float dt;
    int nb, nxb, nzb;
    float c;
    bool verb;

    int it, iz, ix, wfit;                      /* looping index */
    float idz, idx;                            /* inverse */
    float **next,**curr,**lapl,**dcur,**dprv;  /* tmp arrays */
    float *wb;                              /* absorbing boundary parameters */

    float c0=-30./12.,c1=+16./12.,c2=-1./12.;  /* FD coefficients */

    nx  = pars->nx;
    nz  = pars->nz;
    gpz = pars->gpz;
    dx  = pars->dx;
    dz  = pars->dz;
    nt  = pars->nt;
    snap= pars->snap;
    wfnt= pars->wfnt;
    dt  = pars->dt;
    nb  = pars->nb;
    nxb = pars->nxb;
    nzb = pars->nzb;
    c   = pars->c;
    verb= pars->verb;

    idz = 1/(dz*dz);
    idx = 1/(dx*dx);

    /* allocate temporary arrays */
    next=sf_floatalloc2(nzb,nxb);
    curr=sf_floatalloc2(nzb,nxb);
    dprv=sf_floatalloc2(nzb,nxb);
    dcur=sf_floatalloc2(nzb,nxb);
    lapl=sf_floatalloc2(nzb,nxb);
    for (iz=0; iz<nzb; iz++) {
	for (ix=0; ix<nxb; ix++) {
	    next[ix][iz]=0.f;
	    curr[ix][iz]=0.f;
	    dprv[ix][iz]=0.f;
	    dcur[ix][iz]=0.f;
	    lapl[ix][iz]=0.f;
	}
    }

    /* setting up ABC */
    wb =  nb? sf_floatalloc(nb): NULL;
    if (nb) { /* find absorbing coefficients */
	for(iz=0; iz<nb; iz++){
	    wb[iz]=exp(-c*c*(nb-1-iz)*(nb-1-iz));
	}
    }

    /* modeling/migration */
    if (adj) { /* migration <- read data */

    if (snap>0) wfit = (int)(nt-1)/snap; // wfnt-1
    /* time stepping */
    for (it=nt-1; it > -1; it--) {
      if (verb) sf_warning("it=%d/%d;",it,nt-1);

      #ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(curr,lapl)
#endif
	/* 4th order laplacian */
	for (iz=2; iz<nzb-2; iz++) {
	    for (ix=2; ix<nxb-2; ix++) {
		lapl[ix][iz] = 
		    c0* curr[ix  ][iz  ] * (idx+idz) + 
		    c1*(curr[ix-1][iz  ] + curr[ix+1][iz  ])*idx +
		    c2*(curr[ix-2][iz  ] + curr[ix+2][iz  ])*idx +
		    c1*(curr[ix  ][iz-1] + curr[ix  ][iz+1])*idz +
		    c2*(curr[ix  ][iz-2] + curr[ix  ][iz+2])*idz;
	    }
	}

	/* time stepping */
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(dcur,dprv,lapl,next,curr)
#endif
	for (iz=0; iz<nzb; iz++) {
	    for (ix=0; ix<nxb; ix++) {
		dcur[ix][iz] = dprv[ix][iz] + lapl[ix][iz]*vv[ix][iz]*vv[ix][iz]*dt;
		next[ix][iz] = curr[ix][iz] + dcur[ix][iz]*dt;
		dprv[ix][iz] = dcur[ix][iz]; 
		curr[ix][iz] = next[ix][iz];
 	    }
	}
	
	if (nb) {
	    /* absorbing boundary */
	    for (ix=0; ix < nb; ix++) {  
		for (iz=2; iz < nzb-2; iz++) {
		    curr[ix+2][iz] *= wb[ix];
		    curr[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		    dprv[ix+2][iz] *= wb[ix];
		    dprv[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		}
	    }
	    for (ix=2; ix < nxb-2; ix++) {  
		for (iz=0; iz < nb; iz++) {
		    curr[ix][iz+2] *= wb[iz];
		    curr[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		    dprv[ix][iz+2] *= wb[iz];
		    dprv[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		}
	    }
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	/* data injection */
	for (ix=0; ix < nx; ix++) {
	    curr[ix+nb+2][gpz+nb+2] += dat[ix][it];
	}

	if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	    for ( ix = 0; ix < nx; ix++) {
	        for ( iz = 0; iz < nz; iz++ ) { 
		    wvfld[wfit][ix][iz] = curr[ix+nb+2][iz+nb+2];
		}
	    }
	    if (snap>0) wfit--;
	}
    } /*time stepping*/

    /*generate image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
	    img[ix][iz] = curr[ix+nb+2][iz+nb+2];
	}
    }

    } else { /* modeling -> write data */

    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    curr[ix+nb+2][iz+nb+2]=img[ix][iz];
	}
    }
    if (snap>0) wfit = 0;

    /* Main loop: propagation in time */
    for (it=0; it < nt; it++) {
        if (verb) sf_warning("it=%d/%d;",it,nt-1);

	/* record data */
	for (ix=0; ix<nx; ix++){
	    dat[ix][it] = curr[ix+nb+2][gpz+nb+2];
	}

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(curr,lapl)
#endif
	/* 4th order laplacian */
	for (iz=2; iz<nzb-2; iz++) {
	    for (ix=2; ix<nxb-2; ix++) {
		lapl[ix][iz] = 
		    c0* curr[ix  ][iz  ] * (idx+idz) + 
		    c1*(curr[ix-1][iz  ] + curr[ix+1][iz  ])*idx +
		    c2*(curr[ix-2][iz  ] + curr[ix+2][iz  ])*idx +
		    c1*(curr[ix  ][iz-1] + curr[ix  ][iz+1])*idz +
		    c2*(curr[ix  ][iz-2] + curr[ix  ][iz+2])*idz;
	    }
	}

	/* time step */
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(dcur,dprv,lapl,next,curr)
#endif
	for (iz=0; iz<nzb; iz++) {
	    for (ix=0; ix<nxb; ix++) {
		dcur[ix][iz] = dprv[ix][iz] + lapl[ix][iz]*vv[ix][iz]*vv[ix][iz]*dt;
		next[ix][iz] = curr[ix][iz] + dcur[ix][iz]*dt;
		dprv[ix][iz] = dcur[ix][iz]; 
		curr[ix][iz] = next[ix][iz];
 	    }
	}
	
	if (nb) {
	    /* absorbing boundary */
	    for (ix=0; ix < nb; ix++) {  
		for (iz=2; iz < nzb-2; iz++) {
		    curr[ix+2][iz] *= wb[ix];
		    curr[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		    dprv[ix+2][iz] *= wb[ix];
		    dprv[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		}
	    }
	    for (ix=2; ix < nxb-2; ix++) {  
		for (iz=0; iz < nb; iz++) {
		    curr[ix][iz+2] *= wb[iz];
		    curr[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		    dprv[ix][iz+2] *= wb[iz];
		    dprv[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		}
	    }
	}
	
	if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	    for ( ix = 0; ix < nx; ix++) {
	        for ( iz = 0; iz < nz; iz++ ) { 
		    wvfld[wfit][ix][iz] = curr[ix+nb+2][iz+nb+2];
		}
	    }
	    wfit++;
	}

    }
    if (verb) sf_warning(".");

    } /* mig/mod */

    return 0;
}

int main(int argc, char* argv[])
{
    int nx, nz, gpz;
    float dx, dz;
    int nt, snap, wfnt;
    float dt;
    int nb, nxb, nzb;
    float c;
    bool verb;

    bool adj;
    int ntx,nzx,ix,iz;
    float ox, oz;
    float **vel, **img, **dat, ***wvfld;
    sf_file velocity, data, image, snaps;
    par pars;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&adj)) adj=false; /* if n, modeling; if y, migration */
    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getint("snap",&snap)) snap=0;   /* interval for snapshots */
    if(!sf_getint("gpz",&gpz)) gpz=0;       /* geophone surface */
    if (!sf_getint("nb",&nb)) nb=30;        /* boundary length */
    if (!sf_getfloat("c",&c)) c = 0.01;     /* decaying parameter */
    if (nb) sf_warning("Absorbing width=%d; decay parameter=%g",nb,c);
    else sf_warning("No ABC applied! \n");

    /* Read/Write axes */
    velocity = sf_input("vel");
    if (!sf_histint(velocity,"n1",&nz)) sf_error("No n1= in vel");
    if (!sf_histfloat(velocity,"d1",&dz)) sf_error("No d1= in vel");
    if (!sf_histfloat(velocity,"o1",&oz)) oz=0.f;

    if (!sf_histint(velocity,"n2",&nx))  sf_error("No n2= in vel");
    if (!sf_histfloat(velocity,"d2",&dx)) sf_error("No d2= in vel");
    if (!sf_histfloat(velocity,"o2",&ox)) ox=0.f;

    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putfloat(image,"o1",oz);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putfloat(image,"o2",ox);
	sf_putstring(image,"label2","Distance");
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");	  /* time samples */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt="); /* rate */

	sf_putint(data,"n1",nt);
	sf_putfloat(data,"d1",dt);
	sf_putfloat(data,"o1",0.);
	sf_putstring(data,"label1","Time");
	sf_putstring(data,"unit1","s");

	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putfloat(data,"o2",ox);
	sf_putstring(data,"label2","Distance");
    }

    if (snap > 0) {
        wfnt = (int)(nt-1)/snap+1;
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	sf_putint(snaps,"n1",nz);
	sf_putfloat(snaps,"d1",dz);
	sf_putfloat(snaps,"o1",oz);
	sf_putstring(snaps,"label1","Depth");
	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putfloat(snaps,"o2",ox);
	sf_putstring(snaps,"label2","Distance");
	sf_putint(snaps,"n3",wfnt);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.f);
	sf_putstring(snaps,"label3","Time");
    } else {
        wfnt = 0;
	snaps = NULL;
    }

    /* for convenience */
    ntx = nt*nx;
    nzx = nz*nx;
    nzb = nz + 2*nb + 4; /* 4 for the fd stencil width */
    nxb = nx + 2*nb + 4;

    /* allocate and initialize arrays */
    vel = sf_floatalloc2(nzb,nxb);
    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nt,nx);
    if (snap > 0) wvfld = sf_floatalloc3(nz,nx,wfnt);
    else wvfld = NULL;
    pars = (par) sf_alloc(1, sizeof(*pars));

    if (adj) {
      sf_floatread(dat[0],ntx,data);
    } else {
      sf_floatread(img[0],nzx,image);
    }

    /*input & extend velocity model*/
    for (ix=nb+2; ix<nx+nb+2; ix++){
        sf_floatread(vel[ix]+nb+2,nz,velocity);
         for (iz=0; iz<nb+2; iz++){
             vel[ix][iz] = vel[ix][nb+2];
             vel[ix][nz+nb+2+iz] = vel[ix][nz+nb+1]; //nb+1=nb+2-1
         }
    }

    for (ix=0; ix<nb+2; ix++){
        for (iz=0; iz<nzb; iz++){
            vel[ix][iz] = vel[nb+2][iz];
            vel[nx+nb+2+ix][iz] = vel[nx+nb+1][iz];
	}
    }

#ifdef _OPENMP
#pragma omp parallel
{   
    if (verb) sf_warning(">>>> Using %d threads <<<<<", omp_get_num_threads());
}
#endif

    /*load constant Par elements*/
    pars->nx  = nx;
    pars->nz  = nz;
    pars->gpz = gpz;
    pars->dx  = dx;
    pars->dz  = dz;
    pars->nt  = nt;
    pars->snap= snap;
    pars->wfnt= wfnt;
    pars->dt  = dt;
    pars->nb  = nb;
    pars->nxb = nxb;
    pars->nzb = nzb;
    pars->c   = c;
    pars->verb= verb;

    fdexp(img, dat, adj, pars, vel, wvfld);

    if (adj) {
	sf_floatwrite(img[0],nzx,image);
    } else {
	sf_floatwrite(dat[0],ntx,data);
    }
    
    if (snap > 0 && NULL != snaps) {
	sf_floatwrite(wvfld[0][0],wfnt*nx*nz,snaps);
    }

    exit(0);
}
