/* Acoustic staggered-gridded time-domain FD modeling,
   automatically determines whether or not to use 3D or 2D.

Acoustic wave equation finite difference modeling in both 2D and 3D, using an explicit time-domain solver.

*** Please see the SConstruct in book/tutorial/ewe for a SConstruct that demonstrates how to use 
predefined functions for using this program. ***

This program solves a system of first-order PDE's for pressure and particle velocity using a staggered-grid approach.
The model parameters are compressibility (1/K K: bulk modulus) and density.

The program is parallelized using OpenMP, so be sure to use a compatible compiler to take
advantage of the performance boost

============= STAGGERED-GRID   ========================

		o -- x -- o -- x -- o -- x -- o
        |    |    |    |    |    |    |
		x -- x -- x -- x -- x -- x -- x
		|    |    |    |    |    |    |
		o -- x -- o -- x -- o -- x -- o

The "o"'s are the points where the pressures at computed (integer grid). The "x"'s 
are the points where the particle velocities are computed (half grid).

============= FILE DESCRIPTIONS   ========================      

Fdat.rsf - An RSF file containing your data in the following format:
            axis 1 - source location
            axis 2 - wavefield component (z,x,y) order
            axis 3 - Time
			
Fwav.rsf - An RSF file containing your VOLUME DENSITY INJECTION RATE AND DENSITY OF FORCE 
           wavelet information.  The sampling interval, origin time, 
           and number of time samples will be used as the defaults for the modeling code.
	       i.e. your wavelet needs to have the same length and parameters that you want to model with!
	       The first axis is the number of source locations.
	       The second axis contains [fz, fx, (fy,) q],respectively. If the file is 1D then the source is assumed
	       to be a isotropic pressure source.
	       The third axis is time.
	       The code check the dimensions of the model and the dimensions of the wavelt file; for 2D modeling, the wavelet
	       file may have n2=1 or n2=3, for 3D modeling, n2=1 or n2=4.  An error is returned if the dimensions don't match.
		   
Fcom.rsf - An N dimensional RSF file that contains the values for the compressibility (inverse of the bulk modulus K) at every point
        in the computational domain.
		
Fden.rsf - An N dimensional RSF file that contains the values for density at every point in the computational domain.

Fsou.rsf, Frec.rsf - The source and receiver RSF files respectively.  
                   The 1st axis contains the locations for the points like so:
				   [x,y,z]
                   The second axis is a concatenated list of all points in the list.
				   So, for an array of receivers, it would look like:
                   [x1,y1,z1]
                   [x2,y2,z2]
                   [x3,y3,z3]
                   [x4,y4,z4]
				   
Fwfl.rsf     - The name of the file to save the PRESSURE wavefield snapshots to.  This will be an N+2 
          dimensional file.  The file will be organized as follows:
              1-2(3) axes, spatial coordinates
              3(4) axis, wavefield value
              4(5) axis, time, sequential snapshots
              ***The parentheses indicate what the axes will be for 3D models.

Fdat.rsf     - The name of the file to save the receiver data to.  The data has the format of:
	      spatial coordinates, then value of the wavefield.  Lastly, time.
		  
======= PARAMETERS ========

free = y/[n]   - Free surface boundary condition (the free surface is for PRESSURE).

abc  = y/[n]   - Absorbing Boundary Conditions (PML).

nb             - thickness of the absorbing boundary  

verb = y/[n]   - verbosity flag


		  		  
======= TIPS ========

If the simulation seems to slow down as it's running, its a pretty
good indication that the simulation has become unstable and is overflowing
with NaNs.

*/
/*
  Copyright (C) 2008 Colorado School of Mines
  
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
#include "omputil.h"
#endif

#include "fdutil.h"

#include <math.h>

/*------------------------------------------------------------*/
/* stencils and stencil parameters*/
#define NOP 4 /* derivative operator half-size */
#define alfa 3e-1 /*PML sigmoid std dev*/
#define mPML 100 /* max value of the PML*/
/* Muir's derivative operator */
#define C1 +0.598144  //  1225/ 1024      /2 
#define C2 -0.039876  // -1225/(1024*  15)/2
#define C3 +0.004785  //  1225/(1024* 125)/2 
#define C4 -0.000348  // -1225/(1024*1715)/2 

/***************/
/* 2D STENCILS */
/***************/
/*  2D forward FD derivative stencils */
#define Fz(a,ix,iz,s) (2*C4*(a[ix  ][iz+4] - a[ix  ][iz-3]) +     \
                       2*C3*(a[ix  ][iz+3] - a[ix  ][iz-2]) +     \
                       2*C2*(a[ix  ][iz+2] - a[ix  ][iz-1]) +     \
                       2*C1*(a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Fx(a,ix,iz,s) (2*C4*(a[ix+4][iz  ] - a[ix-3][iz  ]) +     \
                       2*C3*(a[ix+3][iz  ] - a[ix-2][iz  ]) +     \
                       2*C2*(a[ix+2][iz  ] - a[ix-1][iz  ]) +     \
                       2*C1*(a[ix+1][iz  ] - a[ix  ][iz  ])  )*s

/* 2D backward FD derivative stencils */
#define Bz(a,ix,iz,s) (2*C4*(a[ix  ][iz+3] - a[ix  ][iz-4]) +     \
                       2*C3*(a[ix  ][iz+2] - a[ix  ][iz-3]) +     \
                       2*C2*(a[ix  ][iz+1] - a[ix  ][iz-2]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix  ][iz-1])  )*s
#define Bx(a,ix,iz,s) (2*C4*(a[ix+3][iz  ] - a[ix-4][iz  ]) +     \
                       2*C3*(a[ix+2][iz  ] - a[ix-3][iz  ]) +     \
                       2*C2*(a[ix+1][iz  ] - a[ix-2][iz  ]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix-1][iz  ])  )*s

/***************/
/* 3D STENCILS */
/***************/
/*  3D forward FD derivative stencils */
#define Fz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+4] - a[iy ][ix  ][iz-3]) +     \
                       2*C3*( a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-2]) +     \
                       2*C2*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-1]) +     \
                       2*C1*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz  ])  )*s
#define Fx3(a,iy,ix,iz,s) (2*C4*(a[iy][ix+4][iz  ] - a[iy ][ix-3][iz  ]) +     \
                       2*C3*( a[iy][ix+3][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       2*C2*( a[iy][ix+2][iz  ] - a[iy ][ix-1][iz  ]) +     \
                       2*C1*( a[iy][ix+1][iz  ] - a[iy ][ix  ][iz  ])  )*s
#define Fy3(a,iy,ix,iz,s) (2*C4*(a[iy+4][ix][iz  ] - a[iy-3][ix ][iz  ]) +     \
                       2*C3*( a[iy+3][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       2*C2*( a[iy+2][ix][iz  ] - a[iy-1][ix ][iz  ]) +     \
                       2*C1*( a[iy+1][ix][iz  ] - a[iy  ][ix ][iz  ])  )*s

/* 3D backward FD derivative stencils */
#define Bz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-4]) +     \
                       2*C3*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-3]) +     \
                       2*C2*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz-2]) +     \
                       2*C1*( a[iy][ix  ][iz  ] - a[iy ][ix  ][iz-1])  )*s
#define Bx3(a,iy,ix,iz,s) (2*C4*(a[iy][ix+3][iz  ] - a[iy ][ix-4][iz  ]) +     \
                       2*C3*( a[iy][ix+2][iz  ] - a[iy ][ix-3][iz  ]) +     \
                       2*C2*( a[iy][ix+1][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       2*C1*( a[iy][ix  ][iz  ] - a[iy ][ix-1][iz  ])  )*s
#define By3(a,iy,ix,iz,s) (2*C4*(a[iy+3][ix][iz  ] - a[iy-4][ix ][iz  ]) +     \
                       2*C3*( a[iy+2][ix][iz  ] - a[iy-3][ix ][iz  ]) +     \
                       2*C2*( a[iy+1][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       2*C1*( a[iy  ][ix][iz  ] - a[iy-1][ix ][iz  ])  )*s
/*------------------------------------------------------------*/

int main(int argc, char* argv[])
{

    /* Declare RSF params */
    bool verb,fsrf,snap,abc,is2D,debug;
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* source term file */
    sf_file Fsou=NULL;  /* sources   */
    sf_file Frec=NULL;  /* receivers */
    sf_file Fcom=NULL;    /* compressibility */
    sf_file Fden=NULL;  /* density */
    sf_file Fdat=NULL;  /* PRESSURE data */
    sf_file Fwfl=NULL;  /* wavefield */

/*set all y variables to be either zero or null to avoid compiler warnings
about being uninitialized */
    /* cube axes */
    sf_axis af,at,az,ax,ay=NULL; 
    sf_axis as,ar;
	
    int     nt,nz,nx,ny=0,nf,ns,nr,nb;
    int     it,iz,ix,iy=0;
    float   dt,dz,dx,dy=0,idz,idx,idy=0;

    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy=0;
    float     oqz,oqx,oqy=0;
    float     dqz,dqx,dqy=0;
	
    /* I/O arrays */
    float   *ww=NULL;          /* source term array*/
    float  *dd=NULL;           /* data      */	


    /*------------------------------------------------------------*/
	/* init RSF */
    sf_init(argc,argv);
    
    /*------------------------------------------------------------*/
    /* OpenMP parameters*/
    
	int ompchunk; 
	#ifdef _OPENMP
		int ompnth,ompath;
	#endif

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
    

	/* parse the command line */	
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("abc",&abc)) abc=false; /* "PML ABC" */
	if(! sf_getbool("debug",&debug)) debug=false; /* debug */ 
    
    /*------------------------------------------------------------*/    
    /* I/O files*/
    Fwav = sf_input ("in" ); /* wavelet   */
    Fcom = sf_input ("com"); /* compressibiliy  */
    Fden = sf_input ("den"); /* density   */	
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */


	/* Determine dimensionality, if 2D then axis 3 has n size of 1 */
	sf_axis test = sf_iaxa(Fcom,3);
	if(sf_n(test) == 1) is2D = true;
	else is2D = false;

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    af = sf_iaxa(Fwav,2); sf_setlabel(af,"component"); if(verb) sf_raxa(af); /* time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    az = sf_iaxa(Fcom,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fcom,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    nt = sf_n(at); dt = sf_d(at);
    nf = sf_n(af); /* number of force terms*/
    ns = sf_n(as);
    nr = sf_n(ar);
    nz = sf_n(az); dz = sf_d(az); idz=1/dz;
    nx = sf_n(ax); dx = sf_d(ax); idx=1/dx;


    if(!is2D){ /*If 3D*/
		ay=sf_iaxa(Fcom,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /*space*/
		ny=sf_n(ay); dy=sf_d(ay); idy=1/dy;
	}

	if (debug) {
		fprintf(stderr,"idx = %g\n",idx);
		fprintf(stderr,"idz = %g\n",idz);

		if(!is2D)
			fprintf(stderr,"idy = %g\n",idy);
	}



    /*------------------------------------------------------------*/
    /* wavefield extraction parameters */
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) 
	{  
	   if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
       /* save wavefield every *jsnap* time steps */
	}
    /*------------------------------------------------------------*/	

if (is2D){
    /*****************/
    /* Begin 2D code */
    /*****************/
    
    /* FDM structure */
    fdm2d    fdm;
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
    

    float **tt=NULL;           /* temporary  */
    float **com=NULL;           /* velocity in expanded domain */
    float **ro=NULL;           /* density  in expanded domain */


    /*------------------------------------------------------------*/
    /* pressure: u = U @ t */
    float **u; 
    /* velocity: v = U @ t */
    float **vx;
    float **vz;

    /* PML structures */
    PML2D pml=NULL;

    /* linear interpolation weights/indices */
    lint2d cs,cr;
    /* wavefield snapshot array */
    float     **uc=NULL;	
	
	   
    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,ompchunk);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    
    float *sigma=NULL;
    sigma=sf_floatalloc(nb);
    
    /*------------------------------------------------------------*/

    /* setup output data header */
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);
    sf_setn(af,1);
    sf_setd(af,1);
    sf_oaxa(Fdat,af,3);    

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

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,at,3);
    }


    /* Check the source file dimensions for defining the type of source*/
    if(nf==1 || nf==3) ww = sf_floatalloc(ns); /*source term array allocation */
    if(nf==2 || nf>3) sf_error("The source file does not have the correct dimensions");
    
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    if (debug) fprintf(stderr,"xs = %g\n",ss[0].x);
    if (debug) fprintf(stderr,"zs = %g\n",ss[0].z);


    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/ 
    /* input density */
    tt=sf_floatalloc2(nz,   nx   ); 
    ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    sf_floatread(tt[0],nz*nx,Fden); 
    expand(tt,ro,fdm);

	/*inverse density, we don't like to divide */
	for (ix=0; ix<fdm->nxpad; ix++){
		for(iz=0; iz<fdm->nzpad; iz++){
			ro[ix][iz] = 1/ro[ix][iz];
		}
	}

    /*------------------------------------------------------------*/
    /* input velocity */
    com  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(tt[0],nz*nx,Fcom);
    expand(tt,com,fdm);

    free(*tt); free(tt);

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    u=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    
    /********************************/
    /*PML structures initialization */
    /********************************/
    if (abc){
        if (verb) fprintf(stderr,"initializing pml.. ");    
        pml = pml2d_init(fdm);
        
        
        for (ix=0;ix<nb;ix++){
			//sigma[ix]=mPML*(1-exp(-alfa*(nb-1-ix)*(nb-1-ix)));
			sigma[ix]=mPML*(nb-1-ix)/nb;
			//sigma[ix]=mPML*alfa*(1-tanh(alfa*(nb-1 - (ix + nb/2)))*tanh(alfa*(nb-1 - (ix + nb/2))));
		}
        
        
        if (verb) fprintf(stderr,"DONE\n");
    }
    
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(ix,iz) \
       shared(fdm,u,vx,vz)
    #endif
    
    for(ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++){
	        u[ix][iz]=0;
			vx[ix][iz]=0;
			vz[ix][iz]=0;
		}
	}

    /*------------------------------------------------------------*/
    /*                                                            */
    /*  2D MAIN LOOP                                              */
    /*                                                            */ 
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    /* TIME LOOPING */
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

    /* VELOCITY SOURCE INJECTION */
    if (nf==3)
    {
       /* z component*/
       sf_floatread(ww,ns,Fwav);
       lint2d_inject(vz,ww,cs);
       /* x component*/
       sf_floatread(ww,ns,Fwav);
       lint2d_inject(vx,ww,cs);       
    }

    /**************************/
	/* dv/dt = 1/ro * grad(u) */
	/**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(ix,iz) \
       shared(fdm,vz,vx,u,ro,idx,idz)
    #endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			vx[ix][iz] += (ro[ix][iz])*Bx(u,ix,iz,idx)*dt;
			vz[ix][iz] += (ro[ix][iz])*Bz(u,ix,iz,idz)*dt;			
	    }
	}	


    if (debug) fprintf(stderr,"Just before computing the auxiliary contribution for velocity\n");

    /******************************/
    /* PML FOR THE VELOCITY FIELD */
    /******************************/
    if (abc){
       if (debug) fprintf(stderr,"pml to velocity.. ");
       pml2d_velApply(vz,vx,dt,sigma,fdm);
       if (debug) fprintf(stderr,"DONE\n");
    }
        
    if (debug) fprintf(stderr,"Just before computing the auxiliary contribution for pressure\n");

    /**************************/
    /* du/dt = 1/com * div(v) */
    /**************************/

	/* SOURCE INJECTION */
    sf_floatread(ww,ns,Fwav);
    lint2d_inject(u,ww,cs);
    
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(ix,iz) \
       shared(fdm,u,vx,vz,com,idx,idz,dt)
    #endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
			u[ix][iz] += (com[ix][iz])*( Fx(vx,ix,iz,idx) + Fz(vz,ix,iz,idz) )*dt;
		}
	}

    /*******************************/
    /* PML FOR THE PRESSURE FIELD  */
    /*******************************/
    if (abc){
    /************************************************************************/
       if (debug) fprintf(stderr,"pml to pressure.. ");
       pml2d_presApply(u,vx,vz,dt,pml,com,sigma,fdm);
       if (debug) fprintf(stderr,"DONE\n");
    /************************************************************************/
    }
	/* BOUNDARY CONDITIONS*/
	if (fsrf){
    
       for (ix=NOP; ix<fdm->nxpad-NOP; ix++){
          for (iz=NOP; iz<nb; iz++){
             u[ix][iz]=0;
          }
       }
       
	} /* end free surface boundary conditions */

	/* extract data */
	lint2d_extract(u,dd,cr);

	if(snap && it%jsnap==0) {
	    cut2d(u,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
	if(        it%jdata==0) sf_floatwrite(dd,nr,Fdat);


    /*------------------------------------------------------------*/
    }
    if (verb) fprintf(stderr,"Deallocating arrays.. ");
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(  *u);   free( u);
    free( *vx);   free(vx);
    free( *vz);   free(vz);

    if (abc) {pml2d_free(pml); free(pml); }
    
    if(snap) { free(*uc); free(uc); }
    
    free(*ro); free(ro);
    free(*com); free(com);


    free(ww);
    free(ss);
    free(rr);
    free(dd);    
    if (verb) fprintf(stderr,"DONE\n");
    
    exit (0);
    /**************/
    /* end 2D code*/
    /**************/
}
else{
    /*****************/
    /* Begin 3D code */
    /*****************/

    /* FDM structure */
    fdm3d    fdm;
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    

    float ***tt=NULL;           /* temporary array  */
    float ***com=NULL;           /* velocity in expanded domain */
    float ***ro=NULL;           /* density  in expanded domain */

    /*------------------------------------------------------------*/
    /* pressure: u = U @ t */
    float ***u; 
    /* velocity: v = U @ t */
    float ***vx,***vy,***vz;
    
	/* PML */
    PML3D pml=NULL;

    /* linear interpolation weights/indices */
    lint3d cs,cr;
    /* wavefield snapshot array */
    float     ***uc=NULL;	

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,ompchunk);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);

    float *sigma=NULL;
    sigma=sf_floatalloc(nb);
    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);
    sf_setn(af,1);
    sf_setd(af,1);
    sf_oaxa(Fdat,af,3);   

    /* setup output wavefield header */
    if(snap) {
		if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
		if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);
		if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);

		if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
		if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
		if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);

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

    /* Check the source file dimensions for defining the type of source*/
    if(nf==1 || nf==4) ww = sf_floatalloc(ns); /*source term array allocation */
    if(nf==2 || nf==3 || nf>4) sf_error("The source file does not have the correct dimensions");

    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    if (debug) fprintf(stderr,"xs = %g\n",ss[0].x);
    if (debug) fprintf(stderr,"ys = %g\n",ss[0].y);
    if (debug) fprintf(stderr,"zs = %g\n",ss[0].z);


    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    if (debug) fprintf(stderr,"ns = %d\n",cs->n);
    if (debug) fprintf(stderr,"cr = %d\n",cr->n);

    /*------------------------------------------------------------*/ 
    /* input density */
    tt = sf_floatalloc3(nz, nx, ny); 
    ro = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    sf_floatread(tt[0][0],nz*nx*ny,Fden); 
    expand3d(tt,ro,fdm);

	/*inverse density, we don't like to divide */
	for (iy=0; iy<fdm->nypad; iy++){
		for (ix=0; ix<fdm->nxpad; ix++){
			for(iz=0; iz<fdm->nzpad; iz++){
				ro[iy][ix][iz] = 1/ro[iy][ix][iz];
			}
		}
	}
    /*------------------------------------------------------------*/
    /* input velocity */
    com = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    sf_floatread(tt[0][0],nz*nx*ny,Fcom);
    expand3d(tt,com,fdm);

    free(**tt); free(*tt); free(tt);    

    /********************************/
    /*PML structures initialization */
    /********************************/
    if (abc){
        if (verb) fprintf(stderr,"initializing pml.. ");    
        pml = pml3d_init(fdm);

		/* compute the vector sigma once and for all */
		/* I use ix but it doesn't matter, just a dummy index*/

		for (ix=0;ix<nb;ix++){
			sigma[ix]=mPML*(1-exp(-alfa*(nb-ix)*(nb-ix)));
		}

        if (verb) fprintf(stderr,"DONE\n");
    }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
     u =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vx =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vy =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vz =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

	#ifdef _OPENMP
	#pragma omp parallel for \
		schedule(dynamic,fdm->ompchunk) \
		private(iy,ix,iz) \
		shared(fdm,u,vy,vx,vz)
	#endif
	for(iy=0; iy<fdm->nypad; iy++){
	   for(ix=0; ix<fdm->nxpad; ix++) {
		   for(iz=0; iz<fdm->nzpad; iz++){
					u[iy][ix][iz]=0;
				   vx[iy][ix][iz]=0;
				   vy[iy][ix][iz]=0;
				   vz[iy][ix][iz]=0;
			}
		}
	} 
	
	if (debug) fprintf(stderr,"nxpad %d\n", fdm->nxpad);
	if (debug) fprintf(stderr,"nypad %d\n", fdm->nypad);
	if (debug) fprintf(stderr,"nzpad %d\n", fdm->nzpad);	
	
	if (debug) fprintf(stderr,"sf_n(acz) %d sf_n(acx) %d sf_n(acy) %d\n", sf_n(acz),sf_n(acx),sf_n(acy));
	if (debug) fprintf(stderr,"sf_o(acz) %g sf_o(acx) %g sf_o(acy) %g\n", sf_o(acz),sf_o(acx),sf_o(acy));	
	if (debug) fprintf(stderr,"ozpad %g oxpad %g oypad %g\n", fdm->ozpad,fdm->oxpad,fdm->oypad);		
    /*------------------------------------------------------------*/
    /*                                                            */
	/*  3D MAIN LOOP                                              */
    /*                                                            */ 
    /*------------------------------------------------------------*/

    if(verb) fprintf(stderr,"\n");
    /* TIME LOOPING */
    for (it=0; it<nt; it++) {
	    if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

    /*****************************/
    /* VELOCITY SOURCE INJECTION */
    /*****************************/    
    if (nf==4)
    {
       /* z component*/
       sf_floatread(ww, ns, Fwav);
       lint3d_inject(vz, ww, cs);
       /* x component*/
       sf_floatread(ww, ns, Fwav);
       lint3d_inject(vx, ww, cs);       
       /* y component*/
       sf_floatread(ww, ns, Fwav);
       lint3d_inject(vy, ww, cs);
    }

    /**************************/
	/* dv/dt = 1/ro * grad(u) */
	/**************************/
	
	/**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(iy,ix,iz) \
       shared(fdm,vz,vy,vx,u,ro,idy,idx,idz,dt)
    #endif
    /******************************/	
	for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	  for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
           vx[iy][ix][iz] +=  ro[iy][ix][iz]*Bx3(u,iy,ix,iz,idx)*dt;
           vy[iy][ix][iz] +=  ro[iy][ix][iz]*By3(u,iy,ix,iz,idy)*dt;           
           vz[iy][ix][iz] +=  ro[iy][ix][iz]*Bz3(u,iy,ix,iz,idz)*dt;
	    }
	  }
    }
    
    /******************************/
    /* PML FOR THE VELOCITY FIELD */
    /******************************/
    if (abc){
       if (debug) fprintf(stderr,"pml to velocity.. ");
       pml3d_velApply(vx, vy, vz, dt, sigma, fdm);
       if (debug) fprintf(stderr,"DONE\n");
    }
    
    

    /**************************/
    /* du/dt = 1/com * div(v) */
    /**************************/    

    /********************/
	/* SOURCE INJECTION */
    /********************/	
    sf_floatread(ww,ns,Fwav);
    lint3d_inject(u,ww,cs);    


    /**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(iy,ix,iz) \
       shared(fdm,u,vx,vy,vz,com,idx,idy,idz,dt)
    #endif
   /**************************/
   
   	for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	  for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
			u[iy][ix][iz] += com[iy][ix][iz]*( Fx3(vx,iy,ix,iz,idx) + 
			                                     Fy3(vy,iy,ix,iz,idy) + 
			                                     Fz3(vz,iy,ix,iz,idz) )*dt;
			
			}
		}
    }

    /*******************************/
    /* PML FOR THE PRESSURE FIELD  */
    /*******************************/
    if (abc){
    /************************************************************************/
       if (debug) fprintf(stderr,"pml to pressure.. ");
       pml3d_presApply(u,vx,vy,vz,dt,pml,com,sigma,fdm);
       if (debug) fprintf(stderr,"DONE\n");
    /************************************************************************/
    }
    
    /****************************/
    /* FREE BOUNDARY CONDITIONS */
    /****************************/
    
	if (fsrf){
	   /****************************/
       #ifdef _OPENMP
       #pragma omp parallel for \
          schedule(dynamic,fdm->ompchunk) \
          private(iy,ix,iz) \
          shared(fdm,nb,u)
       #endif
       /****************************/

       for     (iy=0; iy<fdm->nypad; iy++){
         for   (ix=0; ix<fdm->nxpad; ix++){
           for (iz=0; iz<nb; iz++){
               u[iy][ix][iz]=0;
           }
         }
       }

	}


	if(snap && it%jsnap==0) {
	    cut3d(u,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
	}
	if(        it%jdata==0){
			/* extract data */
			lint3d_extract(u,dd,cr);
			sf_floatwrite(dd,nr,Fdat);
	} 

    /*------------------------------------------------------------*/ 
    /* end 3D main loop*/
    }


    /*------------------------------------------------------------*/
    /* deallocate arrays */
    if (verb) fprintf(stderr,"Deallocating arrays.. ");
    
    free(   **u);   free(  *u);   free( u);
    free(  **vx);   free( *vx);   free(vx);
    free(  **vy);   free( *vy);   free(vy);
    free(  **vz);   free( *vz);   free(vz);

    if (abc) {pml3d_free(pml); free(pml); }

    if (snap) { free(**uc); free(*uc); free(uc); }

    free( **ro);  free( *ro);  free( ro);
    free(**com);  free(*com);  free(com);

	free(sigma);

    free(ww);
    free(ss);
    free(rr);
    free(dd);

    if (verb) fprintf(stderr,"DONE\n ");
    
    exit(0);
    /**************/
    /* end 3D code*/
    /**************/
}
/***************/
/* END PROGRAM */
/***************/
}
