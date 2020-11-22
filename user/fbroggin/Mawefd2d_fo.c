/* Finite-difference time-domain (FDTD) wave propagation modeling in lossless acoustic 2D media.

 This program fdelmodc can be used to model waves conforming the 2D wave equation in different media.
 This program computes a solution of the 2D acoustic wave equation
 defined through the first-order linearized systems of Newton's and Hooke's law.*/

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
/*#include "omputil.h"*/

/* A eighth order accuracy scheme is used
to compute the first derivative along x and z axis.
The derivative operator half-size is 8/2=2 */
#define NOP 2

/* See "Generation of finite difference formulas
on arbitrarly spaced grids" by Bengt Fornberg */
#define C1 +1.196289062	/* +1225/1024 */
#define C2 -0.079752604	/* -1225/(1024*15) */
#define C3 +0.009570312 /* +1225/(1024*125) */
#define C4 -0.000697545 /* -1225/(1024*1715) */
#define D1 +1.125 
#define D2 -0.0416666666666666666666667

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc,recvz; 
    int  jsnap,ntsnap,jdata,srctype;

    /* I/O files */
    sf_file Fwav; /* wavelet   */
    sf_file Fsou; /* sources   */
    sf_file Frec; /* receivers */
    sf_file Fvel; /* velocity  */
    sf_file Fden; /* density   */
    sf_file Fdat; /* pressure data */
    sf_file Fdatvz = NULL; /* vertical particle velocity data */
    sf_file Fwfl = NULL; /* wavefield */

    /* Cube axes */
    sf_axis at,az,ax;
    sf_axis as,ar;

    int     nt,nz,nx,ns,nr,nb,sizem;
    int     it,iz,ix;
    float   dt,dz,dx,idz,idx,ww_avg;

    /* FDM structure */
    fdm2d    fdm;
    sponge   spo = NULL;

    /* I/O arrays */
    float  *ww;           /* wavelet   */
    pt2d   *ss;           /* sources   */
    pt2d   *rr;           /* receivers */
    float  *dd;           /* data      */

    float * ro;            /* density */
    float * vel;           /* velocity */    

    float * vele;		/* velocity with expanded dimensions */
    float * roe;		/* density  with expanded dimensions */
    
    float * roz;		/* density for the FD scheme */
    float * rox;		/* density for the FD scheme */

    float * vel2ro;           /* vel^2*ro = 1/k (k is the compressibility) */

    float * vx,* vz,* p; /* x and y velocities and pressure fields */

    /* Linear interpolation weights/indices */
    lint2d cs,cr;

    /* FD operator size */
    float cax,cbx,caz,cbz;

    /* Wavefield cut params */
    sf_axis   acz = NULL,acx = NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     *uc = NULL;
    
    float scal;
    
    int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;
    
   	float vel2, lamda2mu;
	float bx, bz;
    
    /*------------------------------------------------------------*/
    /* Initialize RSF parameters 								  */
    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* Initialize OMP parameters */
    /*------------------------------------------------------------*/
#ifdef _OPENMP
    omp_init();
#endif
    
	/*------------------------------------------------------------*/
	/* Flags 													  */
	/*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* exploding reflector */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    if(! sf_getbool("recvz",&recvz)) recvz=false; /* vertical particle velocity data */

    /*------------------------------------------------------------*/
    /* I/O files 												  */
    /*------------------------------------------------------------*/
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fdat = sf_output("out"); /* data      */
    if(recvz) {
		Fdatvz = sf_output("datvz"); /* wavefield */
	}
    if(snap) {
		Fwfl = sf_output("wfl"); /* wavefield */
	}

	/*------------------------------------------------------------*/
	/* Axes */
	/*------------------------------------------------------------*/    
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
    /* Source type												  */
    /*------------------------------------------------------------*/
    if(! sf_getint("srctype",&srctype)) srctype=1;

    /*------------------------------------------------------------*/
    /* Time steps for saving data and wavefield					  */
    /*------------------------------------------------------------*/
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
		if(! sf_getint("jsnap",&jsnap)) jsnap=nt;        
    }
    
    /*------------------------------------------------------------*/
    /* Expand domain for FD operators and ABC					  */ 
    /*------------------------------------------------------------*/
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;
    
    /* Initialize the struct fdm */
    fdm=fdutil_init(verb,fsrf,az,ax,nb,1,NOP);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);

    /*------------------------------------------------------------*/
    /* Setup output data and wavefield header					  */
    /*------------------------------------------------------------*/
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);
	if(recvz) {
    	sf_oaxa(Fdatvz,ar,1);
    	sf_oaxa(Fdatvz,at,2);
	}

    if(snap) {
		if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
		if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);
		
		if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
		if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
		
		dqz=sf_d(az);
		dqx=sf_d(ax);
		
		acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
		acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
		
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
    
    /* Allocate memory for wavelet and receivers array			  */    
	/*ww = sf_floatalloc(ns);*/
	dd  = (float *)calloc(nr,sizeof(float));
	ww  = (float *)calloc(nt*ns,sizeof(float));
	sf_floatread(ww,nt*ns,Fwav);
    
	/*------------------------------------------------------------*/
    /* Setup source/receiver coordinates						  */
    /*------------------------------------------------------------*/
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* Setup FD coefficients									  */
    /*------------------------------------------------------------*/
    
    /* dx and dz are always the same in my examples */
    idz = dt/dz;
    idx = dt/dx;

    cax = D1;
    cbx = D2;
    
    caz = D1;
    cbz = D2;
    
    /*------------------------------------------------------------*/
    /* Read density and velocity and initialize rox, rox, and     */
	/* vel2ro (which is the inverse of the compressibility k)      */
    /*------------------------------------------------------------*/
    
    /* In my code nxpad is nx+2*nb*2*NOP-1 */
    sizem = fdm->nxpad*fdm->nzpad;
    
    vel = (float *)calloc(nz*nx,sizeof(float));    
    ro  = (float *)calloc(nz*nx,sizeof(float));
    
    vele = (float *)calloc(sizem,sizeof(float));
    roe =  (float *)calloc(sizem,sizeof(float));
    
    roz = (float *)calloc(sizem,sizeof(float));
    rox = (float *)calloc(sizem,sizeof(float));
    vel2ro  = (float *)calloc(sizem,sizeof(float));

    /* Read density */
    sf_floatread(ro,nz*nx,Fden);

    /* Expand density to fdm->nzpad*fdm->nxpad */
    expand(ro,roe,fdm);
    free(ro); 

    /* Read velocity */
    sf_floatread(vel,nz*nx,Fvel);
    
    /* Expand velocity to fdm->nzpad*fdm->nxpad */
    expand(vel,vele,fdm);
    free(vel);
    
    /* See Jan Thorbecke's code for comparison */
    /* Vx: rox */
	ioXx=NOP;
	ioXz=ioXx-1;
	/* Vz: roz */
	ioZz=NOP;
	ioZx=ioZz-1;
	/* P: vel2ro */
	ioPx=NOP-1;
	ioPz=ioPx;
	
    /*------------------------------------------------------------*/
    /* Initalize rox and roz									  */
    /*------------------------------------------------------------*/	
	iz = fdm->nzpad-NOP-1;
	for (ix=NOP-1; ix<fdm->nxpad-NOP-1; ix++) {
		vel2 = vele[ix*fdm->nzpad+iz]*vele[ix*fdm->nzpad+iz];
		lamda2mu = vel2*roe[ix*fdm->nzpad+iz];

		bx = 0.5*(roe[ix*fdm->nzpad+iz]+roe[(ix+1)*fdm->nzpad+iz]);
		bz = roe[ix*fdm->nzpad+iz];
		rox[(ix+1)*fdm->nzpad+iz] = idx/bx;
		roz[(ix)*fdm->nzpad+iz+1] = idz/bz;
		vel2ro[(ix)*fdm->nzpad+iz] = idx*lamda2mu;
	}
	
	ix = fdm->nxpad-NOP-1;
	for (iz=NOP-1; iz<fdm->nzpad-NOP-1; iz++) {
		vel2 = vele[ix*fdm->nzpad+iz]*vele[ix*fdm->nzpad+iz];
		lamda2mu = vel2*roe[ix*fdm->nzpad+iz];
		
		bx = roe[ix*fdm->nzpad+iz];
		bz = 0.5*(roe[ix*fdm->nzpad+iz]+roe[ix*fdm->nzpad+iz+1]);
		rox[(ix+1)*fdm->nzpad+iz] = idx/bx;
		roz[(ix)*fdm->nzpad+iz+1] = idz/bz;
		vel2ro[(ix)*fdm->nzpad+iz] = idx*lamda2mu;
	}
	
	ix = fdm->nxpad-NOP-1;
	iz = fdm->nzpad-NOP-1;
	vel2 = vele[ix*fdm->nzpad+iz]*vele[ix*fdm->nzpad+iz];
	lamda2mu = vel2*roe[ix*fdm->nzpad+iz];
	bx = roe[ix*fdm->nzpad+iz];
	bz = roe[ix*fdm->nzpad+iz];
	rox[(ix+1)*fdm->nzpad+iz] = idx/bx;
	roz[(ix)*fdm->nzpad+iz+1] = idz/bz;
	vel2ro[(ix)*fdm->nzpad+iz] = idx*lamda2mu;

	for (ix=NOP-1; ix<fdm->nxpad-NOP-1; ix++) {
		for (iz=NOP-1; iz<fdm->nzpad-NOP-1; iz++) {
			vel2 = vele[ix*fdm->nzpad+iz]*vele[ix*fdm->nzpad+iz];
			lamda2mu = vel2*roe[ix*fdm->nzpad+iz];
						
			bx = 0.5*(roe[ix*fdm->nzpad+iz]+roe[(ix+1)*fdm->nzpad+iz]);
			bz = 0.5*(roe[ix*fdm->nzpad+iz]+roe[ix*fdm->nzpad+iz+1]);
			rox[(ix+1)*fdm->nzpad+iz] = idx/bx;
			roz[(ix)*fdm->nzpad+iz+1] = idz/bz;
			vel2ro[(ix)*fdm->nzpad+iz] = idx*lamda2mu;
		}
	}
    
    /*------------------------------------------------------------*/
    /* Free surface												  */ /* Double check this for loop */
    /*------------------------------------------------------------*/
    if(fsrf) {
		for (ix=0; ix<fdm->nxpad; ix++) {
	    	for (iz=0; iz<fdm->nb+NOP; iz++) {
				vel2ro[ix*fdm->nzpad+iz] = 0;
	    	}
		}
    }
    
    /*------------------------------------------------------------*/
	/* Allocate wavefield arrays 								  */
    /*------------------------------------------------------------*/
	vx = (float *)calloc(sizem,sizeof(float));
	vz = (float *)calloc(sizem,sizeof(float));    
	p  = (float *)calloc(sizem,sizeof(float));
	if(snap) {
		uc = (float *)calloc(sf_n(acz)*sf_n(acx),sizeof(float));
	}
	           
	/*------------------------------------------------------------*/
	/* Initialize absorbing boundry condition 					  */
	/*------------------------------------------------------------*/
    if(dabc) {
		/* One-way abc setup */
		/* abc = abcone2d_make(dt,vele,fsrf,fdm); */
		/* Sponge abc setup */
		spo = sponge_make(fdm);
    }
	
	free(vele);
	free(roe);
	
	/*------------------------------------------------------------*/
	/*------------------------------------------------------------*/
	/* TIME LOOP 												  */
	/*------------------------------------------------------------*/
	/*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"Beginning of the for loop over time\n\n");
    
    for (it=0; it<nt; it++) {

#ifdef _OPENMP
#pragma omp parallel default (shared) \
shared (rox, roz, vel2ro, p, vx, vz) \
shared (fdm, expl, verb, dabc, cs, cr, spo, dd, uc, ww)
#endif
{  
		if(verb && it%100==0) fprintf(stderr,"%5d/%5d\b\b\b\b\b\b\b\b\b\b\b",it,nt);
		
		/*-----------------------------------------------------------------*/
		/* Calculate vx for all grid points except on the virtual boundary */
		/*-----------------------------------------------------------------*/
		#ifdef _OPENMP
		#pragma omp for private (ix, iz) schedule(runtime)
		#endif
		for (ix=ioXx; ix<fdm->nxpad-NOP+1; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
			for (iz=ioXz; iz<fdm->nzpad-NOP; iz++) {
				vx[ix*fdm->nzpad+iz] -= rox[ix*fdm->nzpad+iz]*(
						cax*(p[ix*fdm->nzpad+iz]     - p[(ix-1)*fdm->nzpad+iz]) +
						cbx*(p[(ix+1)*fdm->nzpad+iz] - p[(ix-2)*fdm->nzpad+iz]));
			}
		}

		/*-----------------------------------------------------------------*/
		/* Calculate vz for all grid points except on the virtual boundary */
		/*-----------------------------------------------------------------*/
		#ifdef _OPENMP
		#pragma omp for private (ix, iz) schedule(runtime)
		#endif
		for (ix=ioZx; ix<fdm->nxpad-NOP; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
			for (iz=ioZz; iz<fdm->nzpad-NOP+1; iz++) {
				vz[ix*fdm->nzpad+iz] -= roz[ix*fdm->nzpad+iz]*(
						caz*(p[ix*fdm->nzpad+iz]   - p[ix*fdm->nzpad+iz-1]) +
						cbz*(p[ix*fdm->nzpad+iz+1] - p[ix*fdm->nzpad+iz-2]));
			}
		}

		/* Inject vertical force source */
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{	
			if (srctype == 2) {
				/*sf_floatread(ww,ns,Fwav);*/
				
				for (ix=0; ix<cs->n; ix++) {
					/* Correction factor for the souce wavelet - Due to the finite-difference scheme
					See Jan Thorbecke's code for comparison. */
					
					/* IMPORTANT - I had to put a "-1" in the x index, but I'm not sure if that is correct */
					/* time = it*dt;
					   id1 = floor(time/dt);
					   id2 = id1+1; */
					ww_avg = ww[ix+it*ns];
					/*ww_avg = ww[id1]*(id2-time/dt) + ww[id2]*(time/dt-id1);*/
				    scal = roz[(cs->jx[ix]-0)*fdm->nzpad+cs->jz[ix]+0] / dx;
					
					vz[(cs->jx[ix])*fdm->nzpad+cs->jz[ix]+0] += ww_avg * scal;

    			}
			}
		}
 
		/*----------------------------------------------------------------*/
		/* Step forward in time											  */
		/* Calculate p for all grid points except on the virtual boundary */
		/*----------------------------------------------------------------*/
		#ifdef _OPENMP
		#pragma omp	for private (ix, iz) schedule(runtime)
		#endif
		for (ix=ioPx; ix<fdm->nxpad-NOP; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
			for (iz=ioPz; iz<fdm->nzpad-NOP; iz++) {
				p[ix*fdm->nzpad+iz] -= vel2ro[ix*fdm->nzpad+iz]*(
						cax*(vx[(ix+1)*fdm->nzpad+iz] - vx[ix*fdm->nzpad+iz]) +
						cbx*(vx[(ix+2)*fdm->nzpad+iz] - vx[(ix-1)*fdm->nzpad+iz]) +
						caz*(vz[ix*fdm->nzpad+iz+1]   - vz[ix*fdm->nzpad+iz]) +
						cbz*(vz[ix*fdm->nzpad+iz+2]   - vz[ix*fdm->nzpad+iz-1]));
			}
		}

		/* Inject pressure source */
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{	
			if (srctype == 1) {
				/*sf_floatread(ww,ns,Fwav);*/
				
				for (ix=0; ix<cs->n; ix++) {
					/* Correction factor for the souce wavelet - Due to the finite-difference scheme
					See Jan Thorbecke's code for comparison. */
					/* time = it*dt;
					   id1 = floor(time/dt);
					   id2 = id1+1; */
					ww_avg = ww[ix+it*ns];
					/*ww_avg = ww[id1]*(id2-time/dt) + ww[id2]*(time/dt-id1);*/
					scal = vel2ro[cs->jx[ix]*fdm->nzpad+cs->jz[ix]] / dx;
					
					p[cs->jx[ix]*fdm->nzpad+cs->jz[ix]] += ww_avg * scal;
    			}
			}
		}

		/*------------------------------------------------------------*/
		/* Initialize absorbing boundry condition 					  */
		/*------------------------------------------------------------*/		
		if(dabc) {					
			/* one-way abc apply */
			/*abcone2d_apply(vx,vz,abc,fdm);*/
			
			sponge2d_apply(vx,vz,spo,fdm);
		}
		
		#ifdef _OPENMP
		#pragma omp master
		#endif
		{
			/* Extract data at receivers*/
			lint2d_extract(p,dd,cr,fdm);
			if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

			if (recvz) {
				lint2d_extract(vz,dd,cr,fdm);
				if(it%jdata==0) sf_floatwrite(dd,nr,Fdatvz);
			}

			/* Extract wavefield in the "box" */
			if(snap && it%jsnap==0) {
	    		cut2d(p,uc,fdm,acz,acx);
				sf_floatwrite(uc,sf_n(acz)*sf_n(acx),Fwfl);
			}
		}

} /* End of parallel section */
	
    /*------------------------------------------------------------*/    
    } /* End of time loop */
    /*------------------------------------------------------------*/

    if(verb) fprintf(stderr,"\nEnd of loop over time\n");

    /*------------------------------------------------------------*/
    /* Deallocate arrays 										  */
    /*------------------------------------------------------------*/
    free(vx);
    free(vz);
    free(p);
    if(snap) {
		free(uc);
	}

    free(rox);
    free(roz);
    free(vel2ro);

    free(ww);
    free(ss);
    free(rr);
    free(dd);
    
    sf_fileclose(Fwav);
    sf_fileclose(Fvel);
    sf_fileclose(Fden);
    sf_fileclose(Fsou);
    sf_fileclose(Frec);
    sf_fileclose(Fdat);
    if (recvz) {
    	sf_fileclose(Fdatvz);
	}
    if (snap) {
    	sf_fileclose(Fwfl);
    }

    exit(0);
}

