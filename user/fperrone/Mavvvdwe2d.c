/* 2D acoustic variable-velocity variable-density time-domain FD modeling 

The code uses a standard second-order stencil in time.
The coefficients of the spatial stencil are computed 
by matching the transfer function of the discretized 
first-derivative operator to the ideal response. 
The optimized coefficients minimize dispersion 
given that particular size of the stencil.

The term 
	ro div (1/ro grad (u))
where
	ro   : density
	div  : divergence op
	grad : gradient  op
	u    : wavefield
	
is implemented in order to obtain a positive semi-definite operator.

The code implements both the forward (adj=n) and adjoint (adj=y) modeling operator.

============= FILE DESCRIPTIONS   ========================      

Fdat.rsf - An RSF file containing your data in the following format:
			axis 1 - Receiver location
			axis 2 - Time
			
Fwav.rsf - An RSF file containing your VOLUME DENSITY INJECTION RATE 
           wavelet information.  The sampling interval, origin time, 
           and number of time samples will be used as the defaults for the modeling code.
	       i.e. your wavelet needs to have the same length and parameters that you want to model with!
	       The first axis is the number of source locations.
	       The second axis is time.
		   
Fvel.rsf - An N dimensional RSF file that contains the values for the velocity field at every point in the computational domain.
		
Fden.rsf - An N dimensional RSF file that contains the values for density at every point in the computational domain.

Frec.rsf - Coordinates of the receivers
		axis 1 - (x,z) of the receiver
		axis 2 - receiver index

Fsou.rsf - Coordinate of the sources
		axis 1 - (x,y,z) of the source
		axis 2 - source index

Fwfl.rsf - Output wavefield

verb=y/n - verbose flag

snap=y/n - wavefield snapshots flag

free=y/n - free surface on/off

dabc=y/n - absorbing boundary conditions on/off

jdata    - data sampling 

jsnap    - wavefield snapshots sampling


Author: Francesco Perrone
Date: November 2013
*/

/*
	Copyright (C) 2013 Colorado School of Mines

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

#include "fdlib.h"
#include <time.h>
/* check: dt<= 0.2 * min(dx,dz)/vmin */

#define NOP 3 /* derivative operator half-size */


/* LS coefficients */
#define C1 +1.1989919
#define C2 -0.08024696
#define C3 +0.00855954


int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,dabc,adj; 
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
    float **iro=NULL;           /* density */
    float **uat=NULL;          /* 1st derivative of wavefield */
    float **vp=NULL;           /* velocity */
    float **vt=NULL;           /* temporary vp*vp * dt*dt */

	float **fsrfcells=NULL;		/* ghost cells for free surface BC */

    float **um2,**um1,**uo,**ua,**ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint2d cs,cr;

	/* FD coefficients */
	float c1x,c1z,c2x,c2z,c3x,c3z;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL;

	/* for benchmarking */
	clock_t start_t, end_t;
	float total_t;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* Wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* Free surface flag */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
	if(! sf_getbool( "adj",&adj )) adj=false;  /* Adjoint flag*/
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fdat = sf_output("out"); /* data      */

    if(snap) Fwfl = sf_output("wfl"); /* wavefield */


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
	/* absorbing boundary */
	if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;
    /*------------------------------------------------------------*/
    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */


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
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az); /* Saved wfld window nz */
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax); /* Saved wfld window nx */
	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az); /* Saved wfld window oz */
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax); /* Saved wfld window ox */

	dqz=sf_d(az);
	dqx=sf_d(ax);

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
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

	/* wavelet and data vector */
	ww = sf_floatalloc(ns);
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);
	fdbell_init(5);
    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;

	c1x = C1*idx;
	c1z = C1*idz;
	c2x = C2*idx;
	c2z = C2*idz;
	c3x = C3*idx;
	c3z = C3*idz;
    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 

    vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 


    /* input density */
	ro = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	iro = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	sf_floatread(tt[0],nz*nx,Fden); expand(tt,ro ,fdm);

	/* inverse density to avoid computation on the fly */
	iro[0][0] = 1/ro[0][0];
	for 	(ix=1; ix<fdm->nxpad; ix++) {
		for (iz=1; iz<fdm->nzpad; iz++) {
			iro[ix][iz] = 4./(  2*ro[ix  ][iz  ] + 
								ro[ix-1][iz  ] + 
								ro[ix  ][iz-1] );
		}
	}

    /* input velocity */
    sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vp,fdm);
    /* precompute vp^2 * dt^2 */
	for 	(ix=0; ix<fdm->nxpad; ix++) {
		for (iz=0; iz<fdm->nzpad; iz++) {
			vt[ix][iz] = vp[ix][iz] * vp[ix][iz] * dt * dt ;
		}
    }
    if(fsrf) { /* free surface */
		fsrfcells = sf_floatalloc2(4*NOP,fdm->nxpad); 
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
    um2=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    um1=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uat=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    um2[ix][iz]=0;
	    um1[ix][iz]=0;
	    uo[ix][iz]=0;
	    ua[ix][iz]=0;
	}
    }


    /*------------------------------------------------------------*/
    /*                                                            */
    /*  MAIN LOOP                                                 */
    /*                                                            */
    /*------------------------------------------------------------*/
    if (!adj){
	    if(verb) fprintf(stderr,"\nFORWARD ACOUSTIC VARIABLE-DENSITY WAVE EXTRAPOLATION \n");
		// extrapolation			
		start_t=clock();
		for (it=0; it<nt; it++) {
			if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

			#ifdef _OPENMP
			#pragma omp parallel private(ix,iz) 
			#endif
			{

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for    (ix=0; ix<fdm->nxpad; ix++) {
						for(iz=nb; iz<nb+2*NOP; iz++) {
							fsrfcells[ix][2*NOP+(iz-nb)  ] =  um1[ix][iz];
							fsrfcells[ix][2*NOP-(iz-nb)-1] = -um1[ix][iz];
						}
					}
				}


				// spatial derivatives
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		

					// gather
					uat[ix][iz]  = iro[ix][iz]*(
									c3z*(um1[ix  ][iz+2] - um1[ix  ][iz-3]) +
									c2z*(um1[ix  ][iz+1] - um1[ix  ][iz-2]) +
									c1z*(um1[ix  ][iz  ] - um1[ix  ][iz-1])  
									);
					}
				}

				if (fsrf){
					// free surface
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for(iz=nb-NOP; iz<nb+NOP; iz++) {
							uat[ix][iz]  = iro[ix][iz]*(
									c3z*(fsrfcells[ix  ][2*NOP+(iz-nb)+2] - fsrfcells[ix  ][2*NOP+(iz-nb)-3]) +
									c2z*(fsrfcells[ix  ][2*NOP+(iz-nb)+1] - fsrfcells[ix  ][2*NOP+(iz-nb)-2]) +
									c1z*(fsrfcells[ix  ][2*NOP+(iz-nb)  ] - fsrfcells[ix  ][2*NOP+(iz-nb)-1])  
									);
						}
					}
				}
				

				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		
	
					// scatter
					ua[ix][iz  ]  =   	c1z*(uat[ix][iz  ] -
											uat[ix][iz+1]) +
										c2z*(uat[ix][iz-1] -
											uat[ix][iz+2]) +
										c3z*(uat[ix][iz-2] -
											uat[ix][iz+3]);
					}
				}


				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
					// gather
					uat[ix][iz]  = iro[ix][iz]*(
									c3x*(um1[ix+2][iz] - um1[ix-3][iz]) +
									c2x*(um1[ix+1][iz] - um1[ix-2][iz]) +
									c1x*(um1[ix  ][iz] - um1[ix-1][iz])  
									);
					}
				}

				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
					// scatter
					ua[ix  ][iz]  +=    c1x*(uat[ix  ][iz] -
											 uat[ix+1][iz]) +
										c2x*(uat[ix-1][iz] -
											 uat[ix+2][iz]) +
										c3x*(uat[ix-2][iz] -
											 uat[ix+3][iz]);
					}
				}

				/* step forward in time */
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=0; ix<fdm->nxpad; ix++) {
				    for(iz=0; iz<fdm->nzpad; iz++) {
						uo[ix][iz] = 2*um1[ix][iz] 
									-  um2[ix][iz] 
									-  ro[ix][iz]*vt[ix][iz]*ua[ix][iz];
		    
					}
				}				
			}

			// source injection
			/* inject acceleration source */
			sf_floatread(ww,ns,Fwav);
			lint2d_bell(uo,ww,cs);

			/* extract data at receivers */
			lint2d_extract(uo,dd,cr);
			if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

			/* extract wavefield in the "box" */
			if(snap && it%jsnap==0) {
				cut2d(uo,uc,fdm,acz,acx);
				sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
			}
			
			/* one-way abc  + sponge apply */	
			if(dabc) {
				abcone2d_apply(uo,um1,NOP,abc,fdm);
				sponge2d_apply(um1,spo,fdm);
				sponge2d_apply(uo,spo,fdm);
			}
			
			/* circulate wavefield arrays */
			ut=um2;
			um2=um1;
			um1=uo;
			uo=ut;

			
	    } /* end time loop */
	    end_t = clock();
		if(verb) fprintf(stderr,"\r");

		if (verb){	
			total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
			fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
			fprintf(stderr,"Exiting of the program...\n");
		}
    }
    else{
	    if(verb) fprintf(stderr,"\nADJOINT ACOUSTIC VARIABLE-DENSITY WAVE EXTRAPOLATION \n");
		// extrapolation
		start_t = clock();
		for (it=0; it<nt; it++){
			if(verb) fprintf(stderr,"%d/%d  \n",it,nt);

			// source injection
			/* inject acceleration source */
			sf_floatread(ww,ns,Fwav);	
			lint2d_bell(uo,ww,cs);

			/* extract data at receivers */
			if(it%jdata==0) {
				lint2d_extract(uo,dd,cr);
				sf_floatwrite(dd,nr,Fdat);
			}
			/* extract wavefield in the "box" */

			if(snap && it%jsnap==0) {
				cut2d(uo,uc,fdm,acz,acx);
				sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
			}
	
			#ifdef _OPENMP
			#pragma omp parallel private(ix,iz)
			#endif
			{

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for    (ix=0; ix<fdm->nxpad; ix++) {
						for(iz=nb; iz<nb+2*NOP; iz++) {
							fsrfcells[ix][2*NOP+(iz-nb)  ] = uo[ix][iz];
							fsrfcells[ix][2*NOP-(iz-nb)-1] = -uo[ix][iz];
						}
					}
				}
			
				// spatial derivatives		
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=0; ix<fdm->nxpad; ix++) {
					for(iz=0; iz<fdm->nzpad; iz++) {
						ua[ix][iz] = -ro[ix][iz]*vt[ix][iz]*uo[ix][iz];
					}
				}


				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {

					// gather
					uat[ix][iz]  = iro[ix][iz]*(
									c3z*(ua[ix  ][iz+2] - ua[ix  ][iz-3]) +
									c2z*(ua[ix  ][iz+1] - ua[ix  ][iz-2]) +
									c1z*(ua[ix  ][iz  ] - ua[ix  ][iz-1]) 
									);
					}
				}

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for    (ix=0; ix<fdm->nxpad; ix++) {
						for(iz=nb-2*NOP; iz<nb+2*NOP; iz++) {
							ua[ix][iz] = -ro[ix][iz]*vt[ix][iz]*fsrfcells[ix][2*NOP+(iz-nb)];
						}
					}
				}
									
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {					
					// scatter
					um1[ix][iz  ]  +=   c1z*(uat[ix][iz]  - 
											 uat[ix][iz+1]) +
										c2z*(uat[ix][iz-1]  - 
											 uat[ix][iz+2]) +
										c3z*(uat[ix][iz-2]  - 
											 uat[ix][iz+3]);
					}
				}


				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
					// gather
					uat[ix][iz]  = iro[ix][iz]*(
									c3x*(ua[ix+2][iz] - ua[ix-3][iz]) +
									c2x*(ua[ix+1][iz] - ua[ix-2][iz]) +
									c1x*(ua[ix  ][iz] - ua[ix-1][iz])
									);
					}
				}
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
					for	(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
					// scatter
					um1[ix][iz]    +=   c1x*(uat[ix  ][iz] - 
											 uat[ix+1][iz]) + 
										c2x*(uat[ix-1][iz] - 
											 uat[ix+2][iz]) + 
										c3x*(uat[ix-2][iz] - 
											 uat[ix+3][iz]);
					}
				}
			
				/* spray in time */
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for    (ix=0; ix<fdm->nxpad; ix++) {
					for(iz=0; iz<fdm->nzpad; iz++) {
				
						um1[ix][iz] += 2*uo[ix][iz];
						um2[ix][iz] = -uo[ix][iz];
					}
				}
		
			}

			/* one-way abc  + sponge apply */	
			if(dabc) {
				abcone2d_apply(um1,uo,NOP,abc,fdm);
				sponge2d_apply(uo,spo,fdm);
				sponge2d_apply(um1,spo,fdm);
				sponge2d_apply(um2,spo,fdm);				
			}

			/* circulate wavefield arrays */
			ut=uo;
			uo=um1;
			um1=um2;
			um2=ut;

    	} /* end time loop*/
    	end_t = clock();
		if(verb) fprintf(stderr,"\n");

		if (verb){	
			total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
			fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
			fprintf(stderr,"Exiting of the program...\n");
		}

	}
    if(verb) fprintf(stderr,"\n");

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*um2); free(um2);
    free(*um1); free(um1);
    free(*uo); free(uo);
    free(*ua); free(ua);
	free(*uat); free(uat);

    if(snap) {
        free(*uc); free(uc);
    }

	if (fsrf){
		free(*fsrfcells); free(fsrfcells);
	}

    free(*vt);  free(vt);
	free(*ro); free(ro);
	free(*iro); free(iro);

    free(ww);
    free(ss);
    free(rr);
    free(dd);

	if (dabc){
		free_sponge(spo);
		free_abcone2d(abc);	
	}
	free(fdm);
	/* ------------------------------------------------------------------------------------------ */	
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
