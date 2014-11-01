/* 3D acoustic variable-velocity variable-density time-domain FD modeling 

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

/* Muir's derivative operator */
//#define C1 +0.598144  //  1225/ 1024      /2 
//#define C2 -0.039876  // -1225/(1024*  15)/2
//#define C3 +0.004785  //  1225/(1024* 125)/2 
//#define C4 -0.000348  // -1225/(1024*1715)/2 

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
    float ***ro=NULL;           /* density */
    float ***iro=NULL;           /* density */
    float ***uat=NULL;          /* 1st derivative of wavefield */
    float ***vp=NULL;           /* velocity */
    float ***vt=NULL;           /* temporary vp*vp * dt*dt */

	float ***fsrfcells=NULL;		/* ghost cells for free surface BC */

    float ***um2,***um1,***uo,***ua,***ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint3d cs,cr;

	/* FD coefficients */
	float c1x,c1z,c1y,c2x,c2z,c2y,c3x,c3z,c3y;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***uc=NULL;

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
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    ay = sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* y */

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

	dqz=sf_d(az);
	dqx=sf_d(ax);
	dqy=sf_d(ay);

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); sf_raxa(acy);
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

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

	/* wavelet and data vector */
	ww = sf_floatalloc(ns);
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

	c1x = C1*idx;
	c1y = C1*idy;
	c1z = C1*idz;

	c2x = C2*idx;
	c2y = C2*idy;	
	c2z = C2*idz;

	c3x = C3*idx;
	c3y = C3*idy;
	c3z = C3*idz;
    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz, nx, ny); 

    vp = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad); 
    vt = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad); 


    /* input density */
	ro  = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
	iro = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
	sf_floatread(tt[0][0],nz*nx*ny,Fden); expand3d(tt,ro ,fdm);

	/* inverse density to avoid computation on the fly */
	/* 
	there is 1 shell for i1=0 || i2=0 || i3=0 that is zero,
	no big deal but better to fix it
	*/
	for 		(iy=1; iy<fdm->nypad; iy++) {
		for 	(ix=1; ix<fdm->nxpad; ix++) {
			for (iz=1; iz<fdm->nzpad; iz++) {
				iro[iy][ix][iz] = 6./(  3*ro[iy  ][ix  ][iz  ] + 
									  	  ro[iy  ][ix  ][iz-1] +
									  	  ro[iy  ][ix-1][iz  ] + 
										  ro[iy-1][ix  ][iz  ] );
			}
		}
	}

	/* input velocity */
	sf_floatread(tt[0][0],nz*nx*ny,Fvel );    expand3d(tt,vp,fdm);
	/* precompute vp^2 * dt^2 */
	for 		(iy=0; iy<fdm->nypad; iy++) {
		for 	(ix=0; ix<fdm->nxpad; ix++) {
			for (iz=0; iz<fdm->nzpad; iz++) {
				vt[iy][ix][iz] = vp[iy][ix][iz] * vp[iy][ix][iz] * dt * dt ;
			}
		}
	}
	if(fsrf) { /* free surface */
		fsrfcells = sf_floatalloc3(4*NOP,fdm->nxpad,fdm->nypad); 
	}

	free(**tt); free(*tt); free(tt);    

    /*------------------------------------------------------------*/

    if(dabc) {
		/* one-way abc setup */
		abc = abcone3d_make(NOP,dt,vp,fsrf,fdm);
		/* sponge abc setup */
		spo = sponge_make(fdm->nb);
    }

	free(**vp); free(*vp); free(vp);

    /*------------------------------------------------------------*/

    /* allocate wavefield arrays */
    um2 = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
    um1 = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
    uo  = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
    ua  = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
    uat = sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);

	for 		(iy=0; iy<fdm->nypad; iy++) {
		for 	(ix=0; ix<fdm->nxpad; ix++) {
			for (iz=0; iz<fdm->nzpad; iz++) {
				um2[iy][ix][iz]=0;
				um1[iy][ix][iz]=0;
				 uo[iy][ix][iz]=0;
				 ua[iy][ix][iz]=0;
			}
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
			#pragma omp parallel private(iy,ix,iz) 
			#endif
			{

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for 		(iy=0; iy<fdm->nypad; iy++) {
						for 	(ix=0; ix<fdm->nxpad; ix++) {
							for (iz=nb; iz<nb+2*NOP; iz++) {
								fsrfcells[iy][ix][2*NOP+(iz-nb)  ] =  um1[iy][ix][iz];
								fsrfcells[iy][ix][2*NOP-(iz-nb)-1] = -um1[iy][ix][iz];
							}
						}
					}
				}


				// derivatives in z
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {		

							// gather
							uat[iy][ix][iz]  = iro[iy][ix][iz]*(
										c3z*(um1[iy][ix][iz+2] - um1[iy][ix][iz-3]) +
										c2z*(um1[iy][ix][iz+1] - um1[iy][ix][iz-2]) +
										c1z*(um1[iy][ix][iz  ] - um1[iy][ix][iz-1])  
										);
						}
					}
				}

				if (fsrf){
					// free surface
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
						for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
							for (iz=nb-NOP; iz<nb+NOP; iz++) {
								uat[iy][ix][iz]  = iro[iy][ix][iz]*(
									c3z*(fsrfcells[iy][ix][2*NOP+(iz-nb)+2] - fsrfcells[iy][ix][2*NOP+(iz-nb)-3]) +
									c2z*(fsrfcells[iy][ix][2*NOP+(iz-nb)+1] - fsrfcells[iy][ix][2*NOP+(iz-nb)-2]) +
									c1z*(fsrfcells[iy][ix][2*NOP+(iz-nb)  ] - fsrfcells[iy][ix][2*NOP+(iz-nb)-1])  
									);
							}
						}
					}
				}


				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {		
	
						// scatter
						ua[iy][ix][iz] = c1z*(	uat[iy][ix][iz  ] -
												uat[iy][ix][iz+1]) +
										c2z*(	uat[iy][ix][iz-1] -
												uat[iy][ix][iz+2]) +
										c3z*(	uat[iy][ix][iz-2] -
												uat[iy][ix][iz+3]);
						}
					}
				}

				// derivatives in x
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
						// gather
						uat[iy][ix][iz]  = iro[iy][ix][iz]*(
										c3x*(um1[iy][ix+2][iz] - um1[iy][ix-3][iz]) +
										c2x*(um1[iy][ix+1][iz] - um1[iy][ix-2][iz]) +
										c1x*(um1[iy][ix  ][iz] - um1[iy][ix-1][iz])  
										);
						}
					}
				}

				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
					// scatter
					ua[iy][ix][iz] += c1x*(	uat[iy][ix  ][iz] -
											uat[iy][ix+1][iz]) +
									c2x*(	uat[iy][ix-1][iz] -
											uat[iy][ix+2][iz]) +
									c3x*(	uat[iy][ix-2][iz] -
											uat[iy][ix+3][iz]);
						}
					}
				}
		
				// derivatives in y
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
						// gather
						uat[iy][ix][iz]  = iro[iy][ix][iz]*(
										c3y*(um1[iy+2][ix][iz] - um1[iy-3][ix][iz]) +
										c2y*(um1[iy+1][ix][iz] - um1[iy-2][ix][iz]) +
										c1y*(um1[iy  ][ix][iz] - um1[iy-1][ix][iz])  
										);
						}
					}
				}

				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
					// scatter
					ua[iy][ix][iz] += c1y*(	uat[iy  ][ix][iz] -
											uat[iy+1][ix][iz]) +
									c2y*(	uat[iy-1][ix][iz] -
											uat[iy+2][ix][iz]) +
									c3y*(	uat[iy-2][ix][iz] -
											uat[iy+3][ix][iz]);
						}
					}
				}


				/* step forward in time */
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {	
						uo[iy][ix][iz] = 2*um1[iy][ix][iz] 
									-  um2[iy][ix][iz] 
									-  ro[iy][ix][iz]*vt[iy][ix][iz]*ua[iy][ix][iz];
						}
					}
				}		
			} /* end parallel section */

			// source injection
			/* inject acceleration source */
			sf_floatread(ww,ns,Fwav);
			lint3d_bell(uo,ww,cs);

			/* extract data at receivers */
			lint3d_extract(uo,dd,cr);
			if(it%jdata==0) sf_floatwrite(dd,nr,Fdat);

			/* extract wavefield in the "box" */
			if(snap && it%jsnap==0) {
				cut3d(uo,uc,fdm,acz,acx,acy);
				sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
			}

			/* one-way abc  + sponge apply */	
			if(dabc) {
				abcone3d_apply(uo,um1,NOP,abc,fdm);
				sponge3d_apply(uo,spo,fdm);
				sponge3d_apply(um1,spo,fdm);
			}

			
			/* circulate wavefield arrays */
			ut=um2;
			um2=um1;
			um1=uo;
			uo=ut;
			
	    } /* end time loop */
	    end_t = clock();
		if(verb) fprintf(stderr,"\n");

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
			if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

			// source injection
			/* inject acceleration source */
			sf_floatread(ww,ns,Fwav);	
			lint3d_bell(uo,ww,cs);

			/* extract data at receivers */
			if(it%jdata==0) {
				lint3d_extract(uo,dd,cr);
				sf_floatwrite(dd,nr,Fdat);
			}
			/* extract wavefield in the "box" */

			if(snap && it%jsnap==0) {
				cut3d(uo,uc,fdm,acz,acx,acy);
				sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
			}
	
			#ifdef _OPENMP
			#pragma omp parallel private(iy,ix,iz)
			#endif
			{

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for 		(iy=0; iy<fdm->nypad; iy++) {
						for 	(ix=0; ix<fdm->nxpad; ix++) {
							for (iz=nb; iz<nb+2*NOP; iz++) {
								fsrfcells[iy][ix][2*NOP+(iz-nb)  ] =  uo[iy][ix][iz];
								fsrfcells[iy][ix][2*NOP-(iz-nb)-1] = -uo[iy][ix][iz];
							}
						}
					}
				}

				// scaling		
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
							ua[iy][ix][iz] = -ro[iy][ix][iz]*vt[iy][ix][iz]*uo[iy][ix][iz];
						}
					}
				}


				// derivative in z
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=NOP; iy<fdm->nypad-NOP; iy++) {
					for 	(ix=NOP; ix<fdm->nxpad-NOP; ix++) {
						for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {

						// gather
						uat[iy][ix][iz]  = iro[iy][ix][iz]*(
										c3z*(ua[iy][ix][iz+2] - ua[iy][ix][iz-3]) +
										c2z*(ua[iy][ix][iz+1] - ua[iy][ix][iz-2]) +
										c1z*(ua[iy][ix][iz  ] - ua[iy][ix][iz-1]) 
										);
						}
					}
				}

				if (fsrf){
					/* free surface */
					#ifdef _OPENMP
					#pragma omp for schedule(dynamic,fdm->ompchunk)
					#endif
					for 		(iy=0; iy<fdm->nypad; iy++) {					
						for 	(ix=0; ix<fdm->nxpad; ix++) {
							for (iz=nb-2*NOP; iz<nb+2*NOP; iz++) {
								ua[iy][ix][iz] = -ro[iy][ix][iz]*vt[iy][ix][iz]*
													fsrfcells[iy][ix][2*NOP+(iz-nb)];
							}
						}
					}
				}
									
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
					// scatter
					um1[iy][ix][iz]  +=   c1z*(uat[iy][ix][iz]  - 
											 uat[iy][ix][iz+1]) +
										c2z*(uat[iy][ix][iz-1]  - 
											 uat[iy][ix][iz+2]) +
										c3z*(uat[iy][ix][iz-2]  - 
											 uat[iy][ix][iz+3]);
						}
					}
				}

				// derivative in x
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
					// gather
					uat[iy][ix][iz]  = iro[iy][ix][iz]*(
									c3x*(ua[iy][ix+2][iz] - ua[iy][ix-3][iz]) +
									c2x*(ua[iy][ix+1][iz] - ua[iy][ix-2][iz]) +
									c1x*(ua[iy][ix  ][iz] - ua[iy][ix-1][iz])
									);
						}
					}
				}
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
						// scatter
						um1[iy][ix][iz] += c1x*(uat[iy][ix  ][iz] - 
												uat[iy][ix+1][iz]) + 
										c2x*(	uat[iy][ix-1][iz] - 
												uat[iy][ix+2][iz]) + 
										c3x*(	uat[iy][ix-2][iz] - 
												uat[iy][ix+3][iz]);
						}
					}
				}

				// derivative in y
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
					// gather
					uat[iy][ix][iz]  = iro[iy][ix][iz]*(
									c3y*(ua[iy+2][ix][iz] - ua[iy-3][ix][iz]) +
									c2y*(ua[iy+1][ix][iz] - ua[iy-2][ix][iz]) +
									c1y*(ua[iy  ][ix][iz] - ua[iy-1][ix][iz])
									);
						}
					}
				}
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
						// scatter
						um1[iy][ix][iz] += c1y*(uat[iy  ][ix][iz] - 
												uat[iy+1][ix][iz]) + 
										c2y*(	uat[iy-1][ix][iz] - 
												uat[iy+2][ix][iz]) + 
										c3y*(	uat[iy-2][ix][iz] - 
												uat[iy+3][ix][iz]);
						}
					}
				}

				/* spray in time */
				#ifdef _OPENMP
				#pragma omp for schedule(dynamic,fdm->ompchunk)
				#endif
				for 		(iy=0; iy<fdm->nypad; iy++) {
					for 	(ix=0; ix<fdm->nxpad; ix++) {
						for (iz=0; iz<fdm->nzpad; iz++) {
				
							um1[iy][ix][iz] += 2*uo[iy][ix][iz];
							um2[iy][ix][iz] = -uo[iy][ix][iz];
						}
					}
				}
		
			}

			/* one-way abc  + sponge apply */	
			if(dabc) {
				abcone3d_apply(um1,uo,NOP,abc,fdm);
				sponge3d_apply(uo,spo,fdm);
				sponge3d_apply(um1,spo,fdm);
				sponge3d_apply(um2,spo,fdm);				
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
    free(**um2); free(*um2); free(um2);
    free(**um1); free(*um1); free(um1);
    free(**uo);  free(*uo);  free(uo);
    free(**ua);  free(*ua);  free(ua);
	free(**uat); free(*uat); free(uat);

    if(snap) {
        free(**uc); free(*uc); free(uc);
    }

	if (fsrf){
		free(**fsrfcells); free(*fsrfcells); free(fsrfcells);
	}

    free(**vt);  free(*vt);  free(vt);
	free(**ro);  free(*ro);  free(ro);
	free(**iro); free(*iro); free(iro);

    free(ww);
    free(ss);
    free(rr);
    free(dd);

	if (dabc){
		free(spo);
		free(abc);	
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
