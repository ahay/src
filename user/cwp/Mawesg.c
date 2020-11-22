/* Acoustic staggered-gridded time-domain FD modeling,
   automatically determines whether or not to use 3D or 2D.

Acoustic wave equation finite difference modeling in both 2D and 3D, using an explicit time-domain solver.

*** Please see the SConstruct in book/tutorial/ewe for a SConstruct that demonstrates how to use 
predefined functions for using this program. ***

This program solves a system of first-order PDE's for pressure and particle velocity using a staggered-grid approach.
The model parameters are incompressibility (K: bulk modulus) and density.

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
		   
Fbulk.rsf - An N dimensional RSF file that contains the values for the incompressibility (bulk modulus K) at every point in the computational domain.
		
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
#include <time.h>
#include <math.h>

/*------------------------------------------------------------*/
/* stencils and stencil parameters*/
#define NOP 3 /* derivative operator half-size */
/* PML parameters*/
#define mPML 100 /* max value of the PML*/
/* Muir's derivative operator */
//#define C1 +0.598144  //  1225/ 1024      /2 
//#define C2 -0.039876  // -1225/(1024*  15)/2
//#define C3 +0.004785  //  1225/(1024* 125)/2 
//#define C4 -0.000348  // -1225/(1024*1715)/2 
/* Optimized 6-point 1st derivative stencil */
#define C1 +1.1989919
#define C2 -0.08024696
#define C3 +0.00855954

/***************/
/* 2D STENCILS */
/***************/
// Muir's operators
/* 1D forward FD derivative stencils */
/*#define Fz(a,ix,iz,s) (2*C4*(a[ix  ][iz+4] - a[ix  ][iz-3]) +     \
                       2*C3*(a[ix  ][iz+3] - a[ix  ][iz-2]) +     \
                       2*C2*(a[ix  ][iz+2] - a[ix  ][iz-1]) +     \
                       2*C1*(a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Fx(a,ix,iz,s) (2*C4*(a[ix+4][iz  ] - a[ix-3][iz  ]) +     \
                       2*C3*(a[ix+3][iz  ] - a[ix-2][iz  ]) +     \
                       2*C2*(a[ix+2][iz  ] - a[ix-1][iz  ]) +     \
                       2*C1*(a[ix+1][iz  ] - a[ix  ][iz  ])  )*s
*/

// Optimized LS
/* 1D forward FD derivative stencils */
#define Fz(a,ix,iz,s) (C3*(a[ix  ][iz+3] - a[ix  ][iz-2]) +     \
                       C2*(a[ix  ][iz+2] - a[ix  ][iz-1]) +     \
                       C1*(a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Fx(a,ix,iz,s) (C3*(a[ix+3][iz  ] - a[ix-2][iz  ]) +     \
                       C2*(a[ix+2][iz  ] - a[ix-1][iz  ]) +     \
                       C1*(a[ix+1][iz  ] - a[ix  ][iz  ])  )*s
                       
// Muir's operators
/* 1D backward FD derivative stencils */
/*#define Bz(a,ix,iz,s) (2*C4*(a[ix  ][iz+3] - a[ix  ][iz-4]) +     \
                       2*C3*(a[ix  ][iz+2] - a[ix  ][iz-3]) +     \
                       2*C2*(a[ix  ][iz+1] - a[ix  ][iz-2]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix  ][iz-1])  )*s
#define Bx(a,ix,iz,s) (2*C4*(a[ix+3][iz  ] - a[ix-4][iz  ]) +     \
                       2*C3*(a[ix+2][iz  ] - a[ix-3][iz  ]) +     \
                       2*C2*(a[ix+1][iz  ] - a[ix-2][iz  ]) +     \
                       2*C1*(a[ix  ][iz  ] - a[ix-1][iz  ])  )*s
*/

// Optimized LS
/* 1D forward FD derivative stencils */
#define Bz(a,ix,iz,s) (C3*(a[ix  ][iz+2] - a[ix  ][iz-3]) +     \
                       C2*(a[ix  ][iz+1] - a[ix  ][iz-2]) +     \
                       C1*(a[ix  ][iz  ] - a[ix  ][iz-1])  )*s
#define Bx(a,ix,iz,s) (C3*(a[ix+2][iz  ] - a[ix-3][iz  ]) +     \
                       C2*(a[ix+1][iz  ] - a[ix-2][iz  ]) +     \
                       C1*(a[ix  ][iz  ] - a[ix-1][iz  ])  )*s
/***************/
/* 3D STENCILS */
/***************/
// Muir's operators
/* 1D forward FD derivative stencils */
/*#define Fz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+4] - a[iy ][ix  ][iz-3]) +     \
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
*/

// Optimized LS
/* 1D forward FD derivative stencils */
#define Fz3(a,iy,ix,iz,s) (C3*( a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-2]) +     \
                       C2*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-1]) +     \
                       C1*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz  ])  )*s
#define Fx3(a,iy,ix,iz,s) (C3*( a[iy][ix+3][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       C2*( a[iy][ix+2][iz  ] - a[iy ][ix-1][iz  ]) +     \
                       C1*( a[iy][ix+1][iz  ] - a[iy ][ix  ][iz  ])  )*s
#define Fy3(a,iy,ix,iz,s) (C3*( a[iy+3][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       C2*( a[iy+2][ix][iz  ] - a[iy-1][ix ][iz  ]) +     \
                       C1*( a[iy+1][ix][iz  ] - a[iy  ][ix ][iz  ])  )*s
                       
// Muir's operators
/* 1D backward FD derivative stencils */
/*#define Bz3(a,iy,ix,iz,s) (2*C4*(a[iy][ix  ][iz+3] - a[iy ][ix  ][iz-4]) +     \
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
*/
                       
// LS optimized
/* 1D backward FD derivative stencils */
#define Bz3(a,iy,ix,iz,s) (C3*( a[iy][ix  ][iz+2] - a[iy ][ix  ][iz-3]) +     \
                       C2*( a[iy][ix  ][iz+1] - a[iy ][ix  ][iz-2]) +     \
                       C1*( a[iy][ix  ][iz  ] - a[iy ][ix  ][iz-1])  )*s
#define Bx3(a,iy,ix,iz,s) (C3*( a[iy][ix+2][iz  ] - a[iy ][ix-3][iz  ]) +     \
                       C2*( a[iy][ix+1][iz  ] - a[iy ][ix-2][iz  ]) +     \
                       C1*( a[iy][ix  ][iz  ] - a[iy ][ix-1][iz  ])  )*s
#define By3(a,iy,ix,iz,s) (C3*( a[iy+2][ix][iz  ] - a[iy-3][ix ][iz  ]) +     \
                       C2*( a[iy+1][ix][iz  ] - a[iy-2][ix ][iz  ]) +     \
                       C1*( a[iy  ][ix][iz  ] - a[iy-1][ix ][iz  ])  )*s
/*------------------------------------------------------------*/

int main(int argc, char* argv[])
{

    /* Declare RSF params */
    bool verb,fsrf,snap,abc,abcpml,is2D,debug;
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav =NULL; /* source term file */
    sf_file Fsou =NULL; /* sources   */
    sf_file Frec =NULL; /* receivers */
    sf_file Fbulk=NULL; /* incompressibility */
    sf_file Fden =NULL; /* density */
    sf_file Fdat =NULL; /* PRESSURE data */
    sf_file Fwfl =NULL; /* wavefield */

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
    float  *ww=NULL; /* source term array*/
    float  *dd=NULL; /* data      */	

	/* for benchmarking */
	clock_t start_t, end_t;
	float total_t;

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
    if(! sf_getbool("abc",&abc)) abc=false; /* ABC if the abcpml=n: spongeABC */
    if(! sf_getbool("pml",&abcpml)) abcpml=false; /* "PML ABC" */
	if(! sf_getbool("debug",&debug)) debug=false; /* debug */ 
    
    /*------------------------------------------------------------*/    
    /* I/O files*/
    Fwav  = sf_input  ("in" );  /* wavelet   */
    Fbulk = sf_input  ("bulk"); /* incompressibiliy  */
    Fden  = sf_input  ("den");  /* density   */	
    Fsou  = sf_input  ("sou");  /* sources   */
    Frec  = sf_input  ("rec");  /* receivers */
    Fwfl  = sf_output ("wfl");  /* wavefield */
    Fdat  = sf_output ("out");  /* data      */


	/* Determine dimensionality, if 2D then axis 3 has n size of 1 */
	sf_axis test = sf_iaxa(Fbulk,3);
	if(sf_n(test) == 1) is2D = true;
	else is2D = false;

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    af = sf_iaxa(Fwav,2); sf_setlabel(af,"component"); if(verb) sf_raxa(af); /*time */
    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */
    az = sf_iaxa(Fbulk,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fbulk,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    nt = sf_n(at); dt = sf_d(at);
    nf = sf_n(af); /* number of force terms*/
    ns = sf_n(as);
    nr = sf_n(ar);
    nz = sf_n(az); dz = sf_d(az); idz=1/dz;
    nx = sf_n(ax); dx = sf_d(ax); idx=1/dx;


    if(!is2D){ /*If 3D*/
		ay=sf_iaxa(Fbulk,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /*space*/
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
    pt2d   *ss=NULL;   /* sources   */
    pt2d   *rr=NULL;   /* receivers */
    

    float **tt  =NULL;	/* temporary  */
    float **bulk=NULL;	/* velocity in expanded domain */
	float **iro =NULL;	/* inverse of the density */
	float **ro =NULL;	/* inverse of the density */	
    float **vp  =NULL;	/* auxiliary wavespeed model  in expanded domain*/

    /*------------------------------------------------------------*/
    /* pressure: u = U @ t */
    float **u; 
    /* velocity: v = U @ t */
    float **vx;
    float **vz;

    /* pressure: u = U @ t-1 */
    float **um = NULL; 
    /* velocity: v = U @ t-1 */
    float **vmx = NULL;
    float **vmz = NULL;

	/* PML structures */
	PML2D pml=NULL;
	/* Sponge structures*/
	abcone2d abcp=NULL;
	sponge   spo=NULL;


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
		sf_oaxa(Fwfl,at, 3);
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
	fdbell_init(5);

    /*------------------------------------------------------------*/ 
    /* input density */
    tt =sf_floatalloc2(     nz,        nx   ); 
    iro =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ro =sf_floatalloc2(fdm->nzpad,fdm->nxpad);	
	
	/* read density and expand */
    sf_floatread(tt[0],nz*nx,Fden); 
    expand(tt,ro,fdm);

	/*inverse density, to avoid division in the extrapolation */
	iro[0][0] = 1/ro[0][0];
	for (ix=1; ix<fdm->nxpad; ix++){
		for(iz=1; iz<fdm->nzpad; iz++){
			iro[ix][iz] = 4/( 2*ro[ix  ][iz] + 
								ro[ix-1][iz] + 
								ro[ix][iz-1]);
		}
	}

	free(*ro);free(ro);

    /*------------------------------------------------------------*/
    /* input velocity */
    bulk  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    sf_floatread(tt[0],nz*nx,Fbulk);
    expand(tt,bulk,fdm);

    free(*tt); free(tt);

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    u=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

	if (abc ){
       um=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
       vmx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
       vmz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
	}
    
    if (abc){
    /********************************/
    /*PML structures initialization */
    /********************************/
    	if (abcpml){
	        if (verb) fprintf(stderr,"initializing pml.. ");    
            pml = pml2d_init(fdm);
        
            for (ix=0;ix<nb;ix++)
               sigma[ix]=(mPML*(nb-1-ix)/nb);

            if (verb) fprintf(stderr,"DONE\n");    	
    	}	
    /***************************************************/
    /*One-way ABC and Sponge structures initialization */
    /***************************************************/
		else{
            /* one-way abc setup   */
    	    if (verb) fprintf(stderr,"setup one-way abc.. ");
            vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad);  
            for    (ix=0; ix<fdm->nxpad; ix++) {
                for(iz=0; iz<fdm->nzpad; iz++) {
                   vp[ix][iz] = pow(bulk[ix][iz]*iro[ix][iz],.5);
                }
            }
            abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
            free(*vp); free(vp);
    	
            if (verb) fprintf(stderr,"initializing sponge... ");    
            /* sponge abc setup */
            spo = sponge_make(fdm->nb);
		}
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
			if (abc){
				um[ix][iz]=0;
				vmx[ix][iz]=0;
				vmz[ix][iz]=0;
			}
		}
	}

    /*------------------------------------------------------------*/
    /*                                                            */
    /*  2D MAIN LOOP                                              */
    /*                                                            */ 
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    /* TIME LOOPING */
    /*
    SEVERAL CASES TO REMOVE THE IF-STATEMENTS FROM THE EXTRAPOLATION LOOP
    */
    start_t = clock();
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

    /**************************/
	/* dv/dt = 1/ro * grad(u) */
	/**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(ix,iz) \
       shared(fdm,vz,vx,u,iro,idx,idz)
    #endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
	    	if (abc){
				vmx[ix][iz] = vx[ix][iz];
				vmz[ix][iz] = vz[ix][iz];
	    	}
			vx[ix][iz] += (iro[ix][iz])*Bx(u,ix,iz,idx)*dt;
			vz[ix][iz] += (iro[ix][iz])*Bz(u,ix,iz,idz)*dt;			
	    }
	}	
	
    /* VELOCITY SOURCE INJECTION */
    if (nf==3)
    {
       /* z component*/
       sf_floatread(ww,ns,Fwav);
       //lint2d_inject(vz,ww,cs);
       lint2d_bell(vz,ww,cs);
       /* x component*/
       sf_floatread(ww,ns,Fwav);
       //lint2d_inject(vx,ww,cs);
       lint2d_bell(vx,ww,cs);
    }

    if (debug) fprintf(stderr,"Just before computing the auxiliary contribution for velocity\n");

    if (abc){
	    /******************************/
		/* PML FOR THE VELOCITY FIELD */
    	/******************************/
		if (abcpml){
       		if (debug) fprintf(stderr,"pml to velocity.. ");
       		pml2d_velApply(vz,vx,dt,sigma,fdm);
       		if (debug) fprintf(stderr,"DONE\n");
		}
		/******************************/
    	/* SPONGE FOR THE VELOCITY FIELD */
    	/******************************/
    	else{
    	
	    	/* one-way ABC */
	    	if (debug) fprintf(stderr,"apply one-way ABC to velocity.. ");
		    abcone2d_apply(vz,vmz,NOP,abcp,fdm);
	        abcone2d_apply(vx,vmx,NOP,abcp,fdm);
    	
	       /* sponge setup*/
	       if (debug) fprintf(stderr,"apply sponge to velocity.. ");
	       sponge2d_apply(vz,spo,fdm);
	       sponge2d_apply(vx,spo,fdm);
	       if (debug) fprintf(stderr,"DONE\n");
	    }

    }
        
    if (debug) fprintf(stderr,"Just before computing the auxiliary contribution for pressure\n");

    /**************************/
    /* du/dt = K * div(v) */
    /**************************/

    
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(ix,iz) \
       shared(fdm,u,vx,vz,bulk,idx,idz,dt)
    #endif
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
	    	if (abc){
		    	um[ix][iz] = u[ix][iz];
		    }
		    
			u[ix][iz] += (bulk[ix][iz])*( Fx(vx,ix,iz,idx) + Fz(vz,ix,iz,idz) )*dt;
		}
	}

	/* SOURCE INJECTION */
    sf_floatread(ww,ns,Fwav);
    //lint2d_inject(u,ww,cs);
	lint2d_bell(u,ww,cs);

	
    if (abc){
		/*******************************/
    	/* PML FOR THE PRESSURE FIELD  */
    	/*******************************/
		if (abcpml){
       		if (debug) fprintf(stderr,"pml to pressure.. ");
       		pml2d_presApply(u,vx,vz,dt,pml,bulk,sigma,fdm);
       		if (debug) fprintf(stderr,"DONE\n");
		}

		/******************************/
		/* SPONGE FOR THE PRESSURE FIELD */
		/******************************/
		else{   	
	       /* one-way ABC */
	       if (debug) fprintf(stderr,"apply one-way ABC to pressure.. ");
		   abcone2d_apply(u,um,NOP,abcp,fdm);

			if (debug) fprintf(stderr,"apply sponge to pressure.. ");
			sponge2d_apply(u,spo,fdm);
			if (debug) fprintf(stderr,"DONE\n");
		}
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
	end_t = clock();
	if(verb) fprintf(stderr,"\n");

	if (verb){	
		total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
		fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
		fprintf(stderr,"Exiting of the program...\n");
	}
    if (verb) fprintf(stderr,"Deallocating arrays.. ");
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(  *u);   free( u);
    free( *vx);   free(vx);
    free( *vz);   free(vz);



    if (abc){
		/*free memory PML */
    	if (abcpml){
			pml2d_free(pml); free(pml);     	
    	}
		/*free memory for the sponge and one-way ABC*/
	    free(  *um);   free( um);
    	free( *vmx);   free(vmx);
    	free( *vmz);   free(vmz);

    }
    
    if(snap) { free(*uc); free(uc); }
    
    free(*iro); free(iro);
    free(*bulk); free(bulk);

	free(sigma);

    free(ww);
    free(ss);
    free(rr);
    free(dd);

    free(spo);
	/* --------------------------------------------------------------- */	
	/* CLOSE FILES AND EXIT */
    if (Fwav!=NULL) sf_fileclose(Fwav); 

    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Frec!=NULL) sf_fileclose(Frec);

    if (Fbulk!=NULL) sf_fileclose(Fbulk);
    if (Fden!=NULL) sf_fileclose(Fden);

    if (Fdat!=NULL) sf_fileclose(Fdat);

    if (Fwfl!=NULL) sf_fileclose(Fwfl);

    if (verb) fprintf(stderr,"DONE\n ");

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
    float ***bulk=NULL;           /* velocity in expanded domain */
    float ***ro=NULL;           /* density  in expanded domain */
    float ***iro=NULL;           /* inverse density  in expanded domain */

	float ***vp=NULL;			/* auxiliary velocity for ONE-WAY ABC*/
    /*------------------------------------------------------------*/
    /* pressure: u = U @ t */
    float ***u; 
    /* velocity: v = U @ t */
    float ***vx,***vy,***vz;

    /* pressure: u = U @ t-1 */
    float ***um; 
    /* velocity: v = U @ t-1 */
    float ***vmx,***vmy,***vmz;    

	/* PML */
    PML3D pml=NULL;

    /* One-way ABC and Sponge */
    abcone3d abcp=NULL;
    sponge   spo=NULL;

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
    tt  = sf_floatalloc3(     nz,        nx,        ny   ); 
    iro = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ro = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
	/* read density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden); 
    expand3d(tt,ro,fdm);

	/*inverse density, to avoid division in the extrapolation */
	iro[0][0][0]= 1/ro[0][0][0];
	#ifdef _OPENMP
	#pragma omp parallel for \
		schedule(dynamic,fdm->ompchunk) \
		private(iy,ix,iz)
	#endif	
	for (iy=1; iy<fdm->nypad; iy++){
		for (ix=1; ix<fdm->nxpad; ix++){
			for(iz=1; iz<fdm->nzpad; iz++){
				iro[iy][ix][iz] = 6./(	3*ro[iy  ][ix  ][iz  ] +
										ro[iy  ][ix  ][iz-1] +
										ro[iy  ][ix-1][iz  ] +
										ro[iy-1][ix  ][iz  ]);
			}
		}
	}

	free(**ro);free(*ro);free(ro);
	
    /*------------------------------------------------------------*/
    /* input velocity */
    bulk = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    sf_floatread(tt[0][0],nz*nx*ny,Fbulk);
    expand3d(tt,bulk,fdm);

    free(**tt); free(*tt); free(tt);


    if (abc){
    /********************************/
    /*PML structures initialization */
    /********************************/
        if (abcpml){
            if (verb) fprintf(stderr,"initializing pml.. ");    
            pml = pml3d_init(fdm);

            /* compute the vector sigma once and for all */
            /* I use ix but it doesn't matter, just a dummy index*/

            for (ix=0;ix<nb;ix++){
                sigma[ix]=mPML*(nb-1-ix)/nb;
            }
            if (verb) fprintf(stderr,"DONE\n");
        }

    /************************************/
    /*Sponge structures initialization */
    /***********************************/
        else{
            /* one-way abc setup   */
            if (verb) fprintf(stderr,"initializing one-way ABC... ");    
            vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
			#ifdef _OPENMP
			#pragma omp parallel for \
			schedule(dynamic,fdm->ompchunk) \
			private(iy,ix,iz)
			#endif	
            for        (iy=0; iy<fdm->nypad; iy++) {
                for    (ix=0; ix<fdm->nxpad; ix++) {
                    for(iz=0; iz<fdm->nzpad; iz++) {
                        vp[iy][ix][iz] =
                        1/sqrt(bulk[iy][ix][iz]*(iro[iy][ix][iz]) );
                    }
                }
            }
            abcp = abcone3d_make(NOP,dt,vp,fsrf,fdm);
            free(**vp); free(*vp); free(vp);
    	
            if (verb) fprintf(stderr,"initializing sponge... ");    
            /* sponge abc setup */
            spo = sponge_make(fdm->nb);
        }
    }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
     u =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vx =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vy =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vz =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    
    um =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vmx =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vmy =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vmz =  sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

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
				   
				   um[iy][ix][iz]=0;
				   vmx[iy][ix][iz]=0;
				   vmy[iy][ix][iz]=0;
				   vmz[iy][ix][iz]=0;
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
    start_t = clock();
    for (it=0; it<nt; it++) {
	    if(verb) fprintf(stderr,"%d/%d  \r",it,nt);

    /**************************/
	/* dv/dt = 1/ro * grad(u) */
	/**************************/
	
	/**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(iy,ix,iz) \
       shared(fdm,vz,vy,vx,u,iro,idy,idx,idz,dt)
    #endif
    /******************************/	
	for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	  for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
           vmx[iy][ix][iz] = vx[iy][ix][iz];
           vmy[iy][ix][iz] = vy[iy][ix][iz];
           vmz[iy][ix][iz] = vz[iy][ix][iz];
	    	
	       vx[iy][ix][iz] +=  iro[iy][ix][iz]*Bx3(u,iy,ix,iz,idx)*dt;
           vy[iy][ix][iz] +=  iro[iy][ix][iz]*By3(u,iy,ix,iz,idy)*dt;           
           vz[iy][ix][iz] +=  iro[iy][ix][iz]*Bz3(u,iy,ix,iz,idz)*dt;	    }
	  }
    }

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

    
    /******************************/
    /* PML FOR THE VELOCITY FIELD */
    /******************************/
    if (abc){
    
    if (abcpml){
       if (debug) fprintf(stderr,"pml to velocity.. ");
       pml3d_velApply(vx, vy, vz, dt, sigma, fdm);
       if (debug) fprintf(stderr,"DONE\n");
    }
    /*******************************/
    /* SPONGE FOR THE VELOCITY FIELD  */
    /*******************************/
    else{
    	
    	/* one-way ABC */
    	if (debug) fprintf(stderr,"apply one-way ABC to velocity.. ");
	    abcone3d_apply(vz,vmz,NOP,abcp,fdm);
	    abcone3d_apply(vy,vmy,NOP,abcp,fdm);
        abcone3d_apply(vx,vmx,NOP,abcp,fdm);
    	
       if (debug) fprintf(stderr,"sponge pressure.. ");
       sponge3d_apply(vx,spo,fdm);
       sponge3d_apply(vy,spo,fdm);
       sponge3d_apply(vz,spo,fdm);
       if (debug) fprintf(stderr,"DONE\n");
    }


    }
    
    

    /**************************/
    /* du/dt = K * div(v) */
    /**************************/

    /**************************/
    #ifdef _OPENMP
    #pragma omp parallel for \
       schedule(dynamic,fdm->ompchunk) \
       private(iy,ix,iz) \
       shared(fdm,u,vx,vy,vz,bulk,idx,idy,idz,dt)
    #endif
   /**************************/
   
   	for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	  for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
	    	um[iy][ix][iz] = u[iy][ix][iz];
		
			u[iy][ix][iz] += bulk[iy][ix][iz]*( Fx3(vx,iy,ix,iz,idx) + 
			                                     Fy3(vy,iy,ix,iz,idy) + 
			                                     Fz3(vz,iy,ix,iz,idz) )*dt;
			
			}
		}
    }

    /********************/
	/* SOURCE INJECTION */
    /********************/	
    sf_floatread(ww,ns,Fwav);
    lint3d_inject(u,ww,cs); 

    /*******************************/
    /* PML FOR THE PRESSURE FIELD  */
    /*******************************/
    if (abc){
    if (abcpml){
       /********************************************************************/
       if (debug) fprintf(stderr,"pml to pressure.. ");
       pml3d_presApply(u,vx,vy,vz,dt,pml,bulk,sigma,fdm);
       if (debug) fprintf(stderr,"DONE\n");
       /********************************************************************/
    }
    /*******************************/
    /* SPONGE FOR THE PRESSURE FIELD  */
    /*******************************/
    else{
    	
    	/* one-way ABC */
    	if (debug) fprintf(stderr,"apply one-way ABC to pressure.. ");
	    abcone3d_apply(u,um,NOP,abcp,fdm);
    	
       if (debug) fprintf(stderr,"sponge pressure.. ");
       sponge3d_apply(u,spo,fdm);
       if (debug) fprintf(stderr,"DONE\n");
    }

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

	end_t = clock();
	if(verb) fprintf(stderr,"\n");

	if (verb){	
		total_t = (float)(end_t - start_t) / CLOCKS_PER_SEC;
		fprintf(stderr,"Total time taken by CPU: %g\n", total_t  );
		fprintf(stderr,"Exiting of the program...\n");
	}
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    if (verb) fprintf(stderr,"Deallocating arrays.. ");
    
    free(   **u);   free(  *u);   free( u);
    free(  **vx);   free( *vx);   free(vx);
    free(  **vy);   free( *vy);   free(vy);
    free(  **vz);   free( *vz);   free(vz);

    free(   **um);   free(  *um);   free( um);
    free(  **vmx);   free( *vmx);   free(vmx);
    free(  **vmy);   free( *vmy);   free(vmy);
    free(  **vmz);   free( *vmz);   free(vmz);

    if (abc){
    	/* free memory PML */
    	if (abcpml){    
    		pml3d_free(pml); free(pml);
    	}
    	/* free memory sponge and ONE-WAY ABC */
    }

    if (snap) { free(**uc); free(*uc); free(uc); }

    free( **iro);  free( *iro);  free( iro);
    free(**bulk);  free(*bulk);  free(bulk);

	free(sigma);

    free(ww);
    free(ss);
    free(rr);
    free(dd);

	free(spo);

	/* --------------------------------------------------------------- */	
	/* CLOSE FILES AND EXIT */
    if (Fwav!=NULL) sf_fileclose(Fwav); 

    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Frec!=NULL) sf_fileclose(Frec);

    if (Fbulk!=NULL) sf_fileclose(Fbulk);
    if (Fden!=NULL) sf_fileclose(Fden);

    if (Fdat!=NULL) sf_fileclose(Fdat);

    if (Fwfl!=NULL) sf_fileclose(Fwfl);

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
