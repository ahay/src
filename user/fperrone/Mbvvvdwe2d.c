/* Born variable-density variable-velocity acoustic 2D time-domain FD modeling */
/*
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

The "reflectivity" that is used in the code is intended to be function of the 
change in VELOCITY. In particular, it is supposed to be equal to the product between 
the background and the perturbation in the velocity field, that is, the linear term in
the perturbation when you expand the square of the perturbed velocity
	
	v^2 = (vo + vp)^2 ~ vo^2 + 2*vo*vp
	
by assuming the perturbation is small compared to the background, the term vp^2 
can be neglected. The factor 2 is included in the source injection term in the code.

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

Fref.rsf - Reflectivity (same dimensions of the velocity model)

Frec.rsf - Coordinates of the receivers
		axis 1 - (x,z) of the receiver
		axis 2 - receiver index

Fsou.rsf - Coordinate of the sources
		axis 1 - (x,y,z) of the source
		axis 2 - source index

Fwfl.rsf - Output wavefield

Fliw.rsf - linearized scattered wavefield

Flid.rsf - linearized scattered data

verb=y/n - verbose flag

snap=y/n - wavefield snapshots flag

free=y/n - free surface on/off



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
    bool verb,fsrf,snap,dabc; 
    int  jsnap,ntsnap;
    int  jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */

    sf_file Fvel=NULL; /* velocity  */
    sf_file Fref=NULL; /* reflectivity */
    sf_file Fden=NULL; /* density   */

    sf_file Fdat=NULL; /* data (background)      */
    sf_file Fwfl=NULL; /* wavefield (background) */

    sf_file Flid=NULL; /* data (scattered)      */
    sf_file Fliw=NULL; /* wavefield (scattered) */

    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */

    float **vpin=NULL;         /* velocity  */
    float **roin=NULL;         /* density   */
    float **rfin=NULL;         /* reflectivity */

    float **vp=NULL;           /* velocity     in expanded domain */
    float **ro=NULL;           /* density      in expanded domain */
    float **iro=NULL;          /* inverse density (buoyancy) */

    float **rf=NULL;           /* reflectivity in expanded domain */

    float  *bdd=NULL;          /* data (background) */
    float  *sdd=NULL;          /* data (scattered)  */

    float **vt=NULL;           /* temporary vp*vp * dt*dt */

	float **fsrfbck=NULL;		/* ghost cells for free surface BC */
	float **fsrfsct=NULL;		/* ghost cells for free surface BC */
	
    float **bum,**buo,**bup,**bua,**buat,**but; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **sum,**suo,**sup,**sua,**suat,**sut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* cube axes */
    sf_axis at,a1,a2,as,ar;
    int     nt,n1,n2,ns,nr,nb;
    int     it,i1,i2;
    float   dt,d1,d2,id1,id2,dt2;

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    fdm2d    fdm;
    abcone2d abc;     /* abc */
    sponge spo;

	/* FD coefficients */
	float c1x,c1z,c2x,c2z,c3x,c3z;

    int ompchunk; 
#ifdef _OPENMP
    int ompnth,ompath;
#endif

    sf_axis   ac1=NULL,ac2=NULL;
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
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  
    /* OpenMP data chunk size */
#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;  
    /* OpenMP available threads */

#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* Absorbing BC */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */

    Fwav = sf_input ("in" ); /* wavelet   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */

    Fvel = sf_input ("vel"); /* velocity  */
    Fden = sf_input ("den"); /* density   */
    Fref = sf_input ("ref"); /* reflectivity */

    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    Fliw = sf_output("liw"); /* wavefield (scattered) */
    Flid = sf_output("lid"); /* data (scattered) */

    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    a1 = sf_iaxa(Fvel,1); sf_setlabel(a1,"z"); if(verb) sf_raxa(a1); /* depth */
    a2 = sf_iaxa(Fvel,2); sf_setlabel(a2,"x"); if(verb) sf_raxa(a2); /* space */

    nt = sf_n(at); dt = sf_d(at);
    ns = sf_n(as);
    nr = sf_n(ar);
    n1 = sf_n(a1); d1 = sf_d(a1);
    n2 = sf_n(a2); d2 = sf_d(a2);

    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,a1,a2,nb,ompchunk);

    sf_setn(a1,fdm->nzpad); sf_seto(a1,fdm->ozpad); if(verb) sf_raxa(a1);
    sf_setn(a2,fdm->nxpad); sf_seto(a2,fdm->oxpad); if(verb) sf_raxa(a2);
    /*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Flid,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);
    sf_oaxa(Flid,at,2);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(a1);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(a2);
	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(a1);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(a2);
	dqz=sf_d(a1);
	dqx=sf_d(a2);

	ac1 = sf_maxa(nqz,oqz,dqz);
	ac2 = sf_maxa(nqx,oqx,dqx);

	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(ac1),sf_n(ac2));

	ntsnap=0;
        for(it=0; it<nt; it++) {
            if(it%jsnap==0) ntsnap++;
        }
        sf_setn(at,  ntsnap);
        sf_setd(at,dt*jsnap);
        if(verb) sf_raxa(at);

/*	sf_setn(at,nt/jsnap);
	sf_setd(at,dt*jsnap); */

	sf_oaxa(Fwfl,ac1,1);
	sf_oaxa(Fwfl,ac2,2);
	sf_oaxa(Fwfl,at, 3);

	sf_oaxa(Fliw,ac1,1);
	sf_oaxa(Fliw,ac2,2);
	sf_oaxa(Fliw,at, 3);
    }

	ww = sf_floatalloc(ns);
	
	
    bdd =sf_floatalloc(nr);
    sdd =sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);
	
	bell2d bell;
	bell = fdbell_init(3);
    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    dt2 =    dt*dt;
    id1 = 1/d1;
    id2 = 1/d2;

	c1x = C1*id2;
	c1z = C1*id1;
	c2x = C2*id2;
	c2z = C2*id1;
	c3x = C3*id2;
	c3z = C3*id1;

    /*------------------------------------------------------------*/ 
    /* input density */
    roin=sf_floatalloc2(n1,   n2   ); 
    ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    iro =sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    sf_floatread(roin[0],n1*n2,Fden); 
    expand(roin,ro,fdm);

	/* inverse density to avoid computation on the fly */
	iro[0][0] = 1/ro[0][0];
	for 	(i2=1; i2<fdm->nxpad; i2++) {
		for (i1=1; i1<fdm->nzpad; i1++) {
			iro[i2][i1] = 4./(  2*ro[i2  ][i1  ] + 
								  ro[i2-1][i1  ] + 
								  ro[i2  ][i1-1] );
		}
	}
	
    free(*roin); free(roin);

    /*------------------------------------------------------------*/
    /* input velocity */
    vpin=sf_floatalloc2(n1,   n2   ); 
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(vpin[0],n1*n2,Fvel);
    expand(vpin,vp,fdm);
    free(*vpin); free(vpin);

    /*------------------------------------------------------------*/
    /* input reflectivity */
    rfin=sf_floatalloc2(n1,   n2   ); 
    rf  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(rfin[0],n1*n2,Fref); 
    expand(rfin,rf,fdm);
    free(*rfin); free(rfin);

    for    (i2=0; i2<fdm->nxpad; i2++) {
	for(i1=0; i1<fdm->nzpad; i1++) {
	    vt[i2][i1] = vp[i2][i1] * vp[i2][i1] * dt2;
	}
    }

	/* free surface */
    if(fsrf){
    	fsrfbck = sf_floatalloc2(4*NOP,fdm->nxpad);
    	fsrfsct = sf_floatalloc2(4*NOP,fdm->nxpad); 
    }
    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    bum=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    buo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    bup=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    bua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    buat=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    
    sum=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    suo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    sup=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    sua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    suat=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    
    for    (i2=0; i2<fdm->nxpad; i2++) {
	for(i1=0; i1<fdm->nzpad; i1++) {
	    bum[i2][i1]=0;
	    buo[i2][i1]=0;
	    bup[i2][i1]=0;
	    bua[i2][i1]=0;

	    sum[i2][i1]=0;
	    suo[i2][i1]=0;
	    sup[i2][i1]=0;
	    sua[i2][i1]=0;
	}
    }

    /*------------------------------------------------------------*/
    /* one-way abc setup */
    abc = abcone2d_make(NOP,dt,vp,fsrf,fdm);
    /* sponge abc setup */
    spo = sponge_make(fdm->nb);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
	if(verb) fprintf(stderr,"\nFORWARD BORN ACOUSTIC VARIABLE-DENSITY WAVE EXTRAPOLATION \n");
	/* extrapolation */
	start_t=clock();
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d  \r",it,nt);
	
	#ifdef _OPENMP
	#pragma omp parallel private(i2,i1)
	#endif
	{

	if (fsrf){
		/* free surface */
		#ifdef _OPENMP
		#pragma omp for schedule(dynamic,fdm->ompchunk)
		#endif
		for    (i2=0; i2<fdm->nxpad; i2++) {
			for(i1=nb; i1<nb+2*NOP; i1++) {
				fsrfbck[i2][2*NOP+(i1-nb)  ] =  buo[i2][i1];
				fsrfbck[i2][2*NOP-(i1-nb)-1] = -buo[i2][i1];
				
				fsrfsct[i2][2*NOP+(i1-nb)  ] =  suo[i2][i1];
				fsrfsct[i2][2*NOP-(i1-nb)-1] = -suo[i2][i1];
			}
		}
	}

	// spatial derivatives
	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,fdm->ompchunk)
	#endif
	for    (i2=NOP; i2<fdm->nxpad-NOP; i2++) {
		for(i1=NOP; i1<fdm->nzpad-NOP; i1++) {		

			// gather
			buat[i2][i1]  = iro[i2][i1]*(
							c3z*(buo[i2  ][i1+2] - buo[i2  ][i1-3]) +
							c2z*(buo[i2  ][i1+1] - buo[i2  ][i1-2]) +
							c1z*(buo[i2  ][i1  ] - buo[i2  ][i1-1])  
							);
							
			suat[i2][i1]  = iro[i2][i1]*(
							c3z*(suo[i2  ][i1+2] - suo[i2  ][i1-3]) +
							c2z*(suo[i2  ][i1+1] - suo[i2  ][i1-2]) +
							c1z*(suo[i2  ][i1  ] - suo[i2  ][i1-1])  
							);
		}
	}

	if (fsrf){
		// free surface
		#ifdef _OPENMP
		#pragma omp for schedule(dynamic,fdm->ompchunk)
		#endif
		for    (i2=NOP; i2<fdm->nxpad-NOP; i2++) {
			for(i1=nb-NOP; i1<nb+NOP; i1++) {
			
				buat[i2][i1]  = iro[i2][i1]*(
				c3z*(fsrfbck[i2  ][2*NOP+(i1-nb)+2] - fsrfbck[i2  ][2*NOP+(i1-nb)-3]) +
				c2z*(fsrfbck[i2  ][2*NOP+(i1-nb)+1] - fsrfbck[i2  ][2*NOP+(i1-nb)-2]) +
				c1z*(fsrfbck[i2  ][2*NOP+(i1-nb)  ] - fsrfbck[i2  ][2*NOP+(i1-nb)-1])
				);
				
				suat[i2][i1]  = iro[i2][i1]*(
				c3z*(fsrfsct[i2  ][2*NOP+(i1-nb)+2] - fsrfsct[i2  ][2*NOP+(i1-nb)-3]) +
				c2z*(fsrfsct[i2  ][2*NOP+(i1-nb)+1] - fsrfsct[i2  ][2*NOP+(i1-nb)-2]) +
				c1z*(fsrfsct[i2  ][2*NOP+(i1-nb)  ] - fsrfsct[i2  ][2*NOP+(i1-nb)-1])
				);
			}
		}
	}

	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,fdm->ompchunk)
	#endif
	for    (i2=NOP; i2<fdm->nxpad-NOP; i2++) {
		for(i1=NOP; i1<fdm->nzpad-NOP; i1++) {		
	
			// scatter
			bua[i2][i1  ] = c1z*(buat[i2][i1  ] -
								 buat[i2][i1+1]) +
							c2z*(buat[i2][i1-1] -
								 buat[i2][i1+2]) +
							c3z*(buat[i2][i1-2] -
								 buat[i2][i1+3]);
								 
			sua[i2][i1  ] = c1z*(suat[i2][i1  ] -
								 suat[i2][i1+1]) +
							c2z*(suat[i2][i1-1] -
								 suat[i2][i1+2]) +
							c3z*(suat[i2][i1-2] -
								 suat[i2][i1+3]);
		}
	}


	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,fdm->ompchunk)
	#endif
	for 	(i2=NOP; i2<fdm->nxpad-NOP; i2++) {
		for	(i1=NOP; i1<fdm->nzpad-NOP; i1++) {
			// gather
			buat[i2][i1]  = iro[i2][i1]*(
						c3x*(buo[i2+2][i1] - buo[i2-3][i1]) +
						c2x*(buo[i2+1][i1] - buo[i2-2][i1]) +
						c1x*(buo[i2  ][i1] - buo[i2-1][i1])  
						);

			suat[i2][i1]  = iro[i2][i1]*(
						c3x*(suo[i2+2][i1] - suo[i2-3][i1]) +
						c2x*(suo[i2+1][i1] - suo[i2-2][i1]) +
						c1x*(suo[i2  ][i1] - suo[i2-1][i1])  
						);


		}
	}

	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,fdm->ompchunk)
	#endif
	for 	(i2=NOP; i2<fdm->nxpad-NOP; i2++) {
		for	(i1=NOP; i1<fdm->nzpad-NOP; i1++) {
				// scatter
				bua[i2  ][i1] += c1x*(buat[i2  ][i1] -
									  buat[i2+1][i1]) +
								 c2x*(buat[i2-1][i1] -
									  buat[i2+2][i1]) +
								 c3x*(buat[i2-2][i1] -
									  buat[i2+3][i1]);
									  
				sua[i2  ][i1] += c1x*(suat[i2  ][i1] -
									  suat[i2+1][i1]) +
								 c2x*(suat[i2-1][i1] -
									  suat[i2+2][i1]) +
								 c3x*(suat[i2-2][i1] -
									  suat[i2+3][i1]);									  
									  
		}
	}


	/* step forward in time */
	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,fdm->ompchunk)
	#endif
	for    (i2=0; i2<fdm->nxpad; i2++) {
	    for(i1=0; i1<fdm->nzpad; i1++) {
		bup[i2][i1] = 2*buo[i2][i1] 
		    -           bum[i2][i1] 
		    -           ro[i2][i1]*vt[i2][i1]*bua[i2][i1];

		sup[i2][i1] = 2*suo[i2][i1] 
		    -           sum[i2][i1] 
		    -           ro[i2][i1]*vt[i2][i1]*sua[i2][i1];

	    }
	}

	/* single scattering */
	#ifdef _OPENMP
	#pragma omp for schedule(dynamic,ompchunk)
	#endif
	for     (i2=0; i2<fdm->nxpad; i2++) {
	    for (i1=0; i1<fdm->nzpad; i1++) {
			sup[i2][i1] -= 2*rf[i2][i1]*ro[i2][i1]*bua[i2][i1]*dt2;
	    }
	}


	}	/* end of the parallel section */
	
	/* inject acceleration source */
	sf_floatread(ww,ns,Fwav);	
	lint2d_bell(bup,ww,bell,cs);

	/* extract data */
	lint2d_extract(bup,bdd,cr);
	lint2d_extract(sup,sdd,cr);

	if(snap && it%jsnap==0) {
	    cut2d(bup,uc,fdm,ac1,ac2);
	    sf_floatwrite(uc[0],sf_n(ac1)*sf_n(ac2),Fwfl);

	    cut2d(sup,uc,fdm,ac1,ac2);
	    sf_floatwrite(uc[0],sf_n(ac1)*sf_n(ac2),Fliw);
	}
	if(        it%jdata==0) {
	    sf_floatwrite(bdd,nr,Fdat);
	    sf_floatwrite(sdd,nr,Flid);
	}

	/* one-way abc  + sponge apply */	
	if(dabc) {
		abcone2d_apply(bup,buo,NOP,abc,fdm);
		sponge2d_apply(bup,        spo,fdm);
		sponge2d_apply(buo,        spo,fdm);

		abcone2d_apply(sup,suo,NOP,abc,fdm);
		sponge2d_apply(sup,        spo,fdm);
		sponge2d_apply(suo,        spo,fdm);
	}

	/* circulate wavefield arrays */
	but=bum;
	bum=buo;
	buo=bup;
	bup=but;

	sut=sum;
	sum=suo;
	suo=sup;
	sup=sut;


    } /* end time loop */
	end_t = clock();
	if(verb) fprintf(stderr,"\n");

	if (verb){	
		total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
		fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
		fprintf(stderr,"Exiting of the program...\n");
	}   


    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*bum); free(bum);
    free(*buo); free(buo);
    free(*bup); free(bup);
    free(*bua); free(bua);
    free(*buat); free(buat);

    free(*sum); free(sum);
    free(*suo); free(suo);
    free(*sup); free(sup);
    free(*sua); free(sua);
    free(*suat); free(suat);

    if(snap) {
        free(*uc); free(uc);
    }

	if (fsrf){
		free(*fsrfbck); free(fsrfbck);
		free(*fsrfsct); free(fsrfsct);
	}

	free(*vp);	free(vp);
    free(*vt);  free(vt);

	free(*ro); free(ro);
	free(*iro); free(iro);

	free(*rf);	free(rf);

    free(ww);
    free(ss);
    free(rr);

    free(bdd);
    free(sdd);

	if (dabc){
		free(spo);
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

	if (Fref!=NULL) sf_fileclose(Fref);


    if (Fdat!=NULL) sf_fileclose(Fdat);

    if (Fwfl!=NULL) sf_fileclose(Fwfl);


	if (Fliw!=NULL) sf_fileclose(Fliw);
	if (Flid!=NULL) sf_fileclose(Flid);

    exit (0);
}
