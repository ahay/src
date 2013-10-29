/* Elastic time-domain FD modeling, automatically determines whether or not to use 3D or 2D, supports different types of elastic.

Elastic wave equation finite difference modeling in both 2D and 3D, using an explicit time-domain solver.

*** Please see the SConstruct in book/tutorial/ewe for a SConstruct that demonstrates how to use 
predefined functions for using this program. ***

This program is designed to be as generic as possible, and allows you to use files
with arbitrary models, and arbitrary source and receiver geometries.  Source types are
as generic as possible.  Supports arbitrary types of anisotropy as well.  

The downside to the generality, is that the program is not as performant as dedicated solvers
that are less flexible.  The program is parallelized using OpenMP, so be sure to use a compatible compiler to take
advantage of the performance boost.
=========== Rotated Staggered grid ==========================
        Ux,Uz=====================Ux,Uz
        ||            |             || 
        ||                          ||
        ||            |             ||
        ||             tij          ||
        ||- - - - - - |- - - - - - -|| 
        ||             Cij          ||
        ||            |             ||
        ||                          ||
        ||            |             ||
        Ux,Uz=====================Ux,Uz
===========  OPTIONS  ======================================= 
ani - The type of anisotropy for this simulation.  Valid options:
        For 2D:
        ISO/HTI/VTI = 0
        TTI = 1

        For 3D:
        ISO/HTI/VTI = 0
        TTI    = 1

        VTI, HTI, and Isotropic media are special cases of ISO/HTI/VTI media. 
        TTI media can be represented using TTI media.

cfl   - Execute the CFL check.  If the CFL check fails, then it will cause the program to fail. 
        The CFL check will check both the stability and accuracy conditions for both p-waves and
        s-waves. Depending on the type of anisotropy that you specify, the CFL condition will
        use a safety factor (that you can override if necessary).  
        
        NOTE: the CFL condition will return both minimum and maximum
        constraints on the grid given your velocity model, desired frequency content, and other
        parameters.  IT IS POSSIBLE TO HAVE NO STABLE, AND ACCURATE SOLUTIONS FOR A GIVEN 
        MODEL WITH GIVEN PARAMETERS. THE CFL CONDITION WILL WARN YOU IF THIS IS THE CASE.

        YOU MUST SPECIFY fmax Parameter as well!

         ----- STABILITY ------
        The stability condition is related to the maximum wave speed and minimum grid sampling
        as follows:

        dt < min(dx,dy,dz) / (sqrt(2)*vmax)

        Given a time sampling dt, it is possible to determine the minimum dx,dy,dz for stability.
        vmax is the MAXIMUM velocity for all waves in the model (usually P-wave velocity).

        For elastic FD, the P-wave most greatly influences the stability, as it moves fastest
        on the grid.  

        The stability condition gives us a LOWER bound on the grid sampling for a given dt.

        ------ ACCURACY -------
        The accuracy condition is related to the number of gridpoints per wavelength.  Thus,

        safety*vmin / fmax > N * sqrt(dx^2+dy^2+dz^2) 

        where vmin is the minimum wave velocity in the model (usually S-wave), fmax is some
        relative measure of the maximum frequency of your wavelet (usually 1.5*peak for Ricker), 
        N is the number of points desired per wavelength (5), and safety is a safety factor that 
        is dependent on the type of anisotropy.  

        For elastic FD, the S-wave most greatly impacts the accuracy of the solution, as the S-wave
        is typically much higher frequency and travels at slower wave speeds, meaning shorter 
        wavelengths.  

        The accuracy condition places an UPPER bound on the grid sampling.

        ---- SAFETY FACTOR -----
        The safety factor depends on the type of anisotropy specified, and attempts to place a lower
        bound on the slowest S-wave velocity (guess):

        ISO/HTI/VTI - (3/4)
        TTI    - (1/2)

        You can also override the safety factor using the safety parameter.
safety- Override the safety factor for the CFL condition.  This should be a floating point (0-1.0).
fmax  - An estimate of the highest frequency content in your wavelet (for Ricker use 1.5*peak)

fsrf  - Use a free surface at the top boundary (z=0).  
        WARNING: The free surface condition appears to introduce numerical artifacts into the simulation.  
        USE AT YOUR OWN RISK.

snap  - Save snapshots of the wavefield every n iterations of the modeling program. 

jsnap - Number of iterations between snapshots of the wavefield.  
	    i.e. jsnap=10, means save a snapshot every 10 iterations. 
	    If you had 1000 total iterations, then you would have 100 snapshots total.
	    The default, will output no snapshots.

jdata - Number of time imterations between exporting the data at the receivers.
	    i.e. jdata=1, means save a snapshot every iteration, which should be the default.
	    This can be used to change the sampling of the data to something different from 
        the wavelet/wavefield.

verb  - Print useful information
debug - Print debugging information.  This is more detailed than verbose.

srctype - An integer which determines where the source wavelet is injected
        in the simulation.  Valid options are:  
            0 - Acceleration source
            1 - Displacement source
            2 - Stress source
            3 - Tensor source
        The default option is 2: Acceleration source.
        For Stress, Displacement and Acceleration sources, your wavelet
        needs to have only 3 components (z,x,y).
        For a Tensor source, you must specify wavelet components for 
        all 3 (2D) or 6 (3D) tensor components in the following order:
        2D: tzz, txx, tzx
        3D: tzz, txx, tyy, tyz, tzx, txy

        Hint:  To inject an acoustic source, use a stress source,
            with equal components on all three components.

dabc  - Use a sponge layer to attenuate waves off the edge of the grid.  Use this in 
        combination with the nb parameter.
abcone- In addition to the sponge layer, using a severe ramp at the very edge of the expanded 
        sponge layer to severely attenuate zero-incidence waves at the boundaries. 
        It's not clear if this condition actually affects most computations.
opot  - True: output is second spatial derivative of potentials; False: output wavefield.
nbell - Size of gaussian used to linearly interpolate curves.  A value of 5 seems to work well.  
nb    - Not listed, but is an important parameter.  Allows you to control the size of the sponge 
        layer for the absorbing boundary condition.  If you are getting reflections off the sides, 
        with dabc=y, then make this number larger (int).  This pads the grid by this amount on all sides.  
        For example:
                   |--------------------------|
                   |            ramp layer    |
                   |r |--------------------|  |
                   |a |        nb          |r |
                   |m |      |~~~~~~~~|    |a |
                   |p |      |  MODEL |    |m |
                   |  |  nb  |  SPACE | nb |p |
                   |  |      |~~~~~~~~|    |  |
                   |  |         nb         |  |
                   |  |--------------------|  |
                   |         ramp layer       |
                   |--------------------------| 
nqz, nqx, oqz, oqx, nqy, oqy, - Allows you to set the parameters for the axes.  Leave as defaults.

=============BOUNDARY CONDITIONS ========================

This code enforces a fixed reflecting boundary condition at the 
edge of the computational domain.  The absorbing sponge is used
IN ADDITION to this condition.

=============FILE DESCRIPTIONS   ========================      

Fdat.rsf - An RSF file containing your data in the following format:
            axis 1 - source location
            axis 2 - wavefield component (z,x,y) order
            axis 3 - Time

Fwav.rsf - An RSF file containing your wavelet information.  For elastic modeling, the wavelet needs 
           to have 3 samples on N1 one for each component Z-X-Y (or just Z-X for 2D).  The second 
           axis describes the component as a function of time.  The sampling interval, origin time, 
           and number of time samples will be used as the defaults for the modeling code.
	       i.e. your wavelet needs to have the same length and parameters that you want to model with!
	   Ex:
	   1st axis    index
	   Z component  0     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
	   X component  1     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
	   Y component  2     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
				    2nd axis
       NOTE: For tensor sources, you must have an appropriate number of components.  See srctype for more information.

cccc  - An N+1 dimensional RSF file that contains the values for the stiffness coefficients to be used 
           as the model for the modeling code.  So, for 2D, this would be a 3D array of values concatenated 
           together in the order as described in the anisotropy section.  Each coefficient file contains 
           the value of that coefficient for every point in space. 
           The axes for this file are: Axis 1: Z; Axis 2: X; Axis 3: Y;
    
        The stiffness tensor coefficients are defined uniformly as follows, where 
        --x---y---z--(y)-----(y) describes how the coefficients depend on space.
        |C11 C12 C13 C14 C15 C16|
        |    C22 C23 C24 C25 C26|
        |        C33 C34 C35 C36|
        |            C44 C45 C46|
        |                C55 C56|
        |                    C66|

        The tensor is assumed to be symmetric.  

        Order of the coefficients in the N+1 dimensional file...
        (First coefficient is the first 2D array in the 3D array).
        2D Anisotropy Modes:

        ISO/HTI/VTI: C11, C33, C55, C13
        "TTI:" C11, C13, C15, C33, C35, C55 
        ***TTI basically allows access to all coefs in 2D, but is not really triclinic media
        ------------------------------------------------------------
        (First coefficient is the first 3D array in the 4D array).
        3D Anisotropy Modes:

        ISO/HTI/VTI: C11, C22, C33, C44, C55, C66, C12, C13, C23
        TTI: C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, 
            C35, C36, C44, C45, C46, C55, C56, C66

   
den      - An N dimensional RSF file that contains the valuese for the density to be used for the model.  
           For 2D, this would be a 2D array.  

sou, rec -The source and receiver RSF files respectively.  
          The 1st axis contains the locations for the points like so:
          [x,y,z]
          The second axis is a concatenated list of all points in the list.
          So, for an array of receivers, it would look like:
              [x1,y1,z1]
              [x2,y2,z2]
              [x3,y3,z3]
              [x4,y4,z4]

wfl     - The name of the file to save the wavefield snapshots to.  This will be an N+2 
          dimensional file.  The file will be organized as follows:
              1-2(3) axes, spatial coordinates
              3(4) axis, wavefield components, in the Z,X,(Y) order
              4(5) axis, time, sequential snapshots
              ***The parentheses indicate what the axes will be for 3D models.

dat     - The name of the file to save the receiver data to.  The data has the format of:
	      spatial coordinates, then the data components of the elastic wavefield in the 
	      same order as the wavefield.  Lastly, time.

========== USEFUL COMMANDS  ============================= 	  

To view the wavefield snapshots (2D case):
sfwindow < Fwfl.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

To view the data (2D case):
sfwindow < Fdat.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

========== TROUBLESHOOTING ===============================

If you aren't getting output, or your output is full of Nans, make sure
that you have the proper dimensions for your wavelet files, and that
your input files make sense.

Make sure your source and receiver points are located inside the 
model space, otherwise you will get all NaNs and the simulation will
run forever.

======= TIPS ========

If the simulation seems to slow down as it's running, its a pretty
good indication that the simulation has become unstable and is overflowing
with NaNs.


Modified by Yuting Duan, Colorado School of Mines,2013-10-22. 

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
#define NOP 4 /* derivative operator half-size */

/* Coefficients are:  (-1)^(k+1)/k * (n!)^2 / (n+k)!(n-k)!
for n = order of approximation i.e. eighth order
and k = coefficient number */
/* CENTERED derivatives 8th order */
#define C1 +0.598144530000  /*1.196289060000 */  
#define C2 -0.039876302100  /*0.079752604200 */
#define C3 +0.004785156250  /*0.009570312500 */
#define C4 -0.0003487723215/*0.000697544643*/ 
#define C5 +0.299072265000  /*1.196289060000/4 */  
#define C6 -0.019938151050  /*0.079752604200/4 */
#define C7 +0.002392578125  /*0.009570312500/4 */
#define C8 -0.00017438616075/*0.000697544643/4 */  

#define Dx(a,ix,iz,s) (C4*(a[ix+4][iz+4] + a[ix+4][iz-3] - \
                           a[ix-3][iz+4] - a[ix-3][iz-3]) + \
		                   C3*(a[ix+3][iz+3] + a[ix+3][iz-2] -\
                           a[ix-2][iz+3] - a[ix-2][iz-2]) +	\
		                   C2*(a[ix+2][iz+2] + a[ix+2][iz-1] - \
                           a[ix-1][iz+2] - a[ix-1][iz-1]) +	\
		                   C1*(a[ix+1][iz+1] + a[ix+1][iz  ] - \
                           a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Dz(a,ix,iz,s) (C4*(a[ix+4][iz+4] - a[ix+4][iz-3] + \
                           a[ix-3][iz+4] - a[ix-3][iz-3]) + \
		                   C3*(a[ix+3][iz+3] - a[ix+3][iz-2] + \
                           a[ix-2][iz+3] - a[ix-2][iz-2]) +	\
		                   C2*(a[ix+2][iz+2] - a[ix+2][iz-1] + \
                           a[ix-1][iz+2] - a[ix-1][iz-1]) +	\
		                   C1*(a[ix+1][iz+1] - a[ix+1][iz  ] + \
                           a[ix  ][iz+1] - a[ix  ][iz  ])  )*s

#define Dx3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy-3][ix+4][iz+4] + \
                   			       a[iy+4][ix+4][iz-3] + a[iy-3][ix+4][iz-3] - \
                   			       a[iy-3][ix-3][iz-3] - a[iy+4][ix-3][iz-3] - \
                   				     a[iy-3][ix-3][iz+4] - a[iy+4][ix-3][iz+4]) + \
                   			   C7*(a[iy+3][ix+3][iz+3] + a[iy-2][ix+3][iz+3] + \
                   			       a[iy+3][ix+3][iz-2] + a[iy-2][ix+3][iz-2] - \
                   			       a[iy-2][ix-2][iz-2] - a[iy+3][ix-2][iz-2] - \
                   				     a[iy-2][ix-2][iz+3] - a[iy+3][ix-2][iz+3]) +   \
                   			   C6*(a[iy+2][ix+2][iz+2] + a[iy-1][ix+2][iz+2] + \
                   			       a[iy+2][ix+2][iz-1] + a[iy-1][ix+2][iz-1] - \
                   			       a[iy-1][ix-1][iz-1] - a[iy+2][ix-1][iz-1] - \
                   				     a[iy-1][ix-1][iz+2] - a[iy+2][ix-1][iz+2]) +	\
                   			   C5*(a[iy+1][ix+1][iz+1] + a[iy  ][ix+1][iz+1] + \
                   			       a[iy+1][ix+1][iz  ] + a[iy  ][ix+1][iz  ] - \
                   			       a[iy  ][ix  ][iz  ] - a[iy+1][ix  ][iz  ] - \
                   				     a[iy  ][ix  ][iz+1] - a[iy+1][ix  ][iz+1])  )*s
#define Dy3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy+4][ix-3][iz-3] + \
			                         a[iy+4][ix+4][iz-3] + a[iy+4][ix-3][iz+4] - \
			                         a[iy-3][ix-3][iz-3] - a[iy-3][ix+4][iz+4] - \
	                    			   a[iy-3][ix-3][iz+4] - a[iy-3][ix+4][iz-3]) + \
			                     C7*(a[iy+3][ix+3][iz+3] + a[iy+3][ix-2][iz-2] + \
			                         a[iy+3][ix+3][iz-2] + a[iy+3][ix-2][iz+3] - \
			                         a[iy-2][ix-2][iz-2] - a[iy-2][ix+3][iz+3] - \
			                      	 a[iy-2][ix-2][iz+3] - a[iy-2][ix+3][iz-2]) +   \
			                     C6*(a[iy+2][ix+2][iz+2] + a[iy+2][ix-1][iz-1] + \
			                         a[iy+2][ix+2][iz-1] + a[iy+2][ix-1][iz+2] - \
			                         a[iy-1][ix-1][iz-1] - a[iy-1][ix+2][iz+2] - \
				                       a[iy-1][ix-1][iz+2] - a[iy-1][ix+2][iz-1]) +	\
			                     C5*(a[iy+1][ix+1][iz+1] + a[iy+1][ix  ][iz  ] + \
			                         a[iy+1][ix+1][iz  ] + a[iy+1][ix  ][iz+1] - \
			                         a[iy  ][ix  ][iz  ] - a[iy  ][ix+1][iz+1] - \
				                       a[iy  ][ix  ][iz+1] - a[iy  ][ix+1][iz  ])  )*s
#define Dz3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy-3][ix+4][iz+4] + \
			                         a[iy-3][ix-3][iz+4] + a[iy+4][ix-3][iz+4] - \
			                         a[iy-3][ix-3][iz-3] - a[iy+4][ix-3][iz-3] - \
				                       a[iy+4][ix+4][iz-3] - a[iy-3][ix+4][iz-3]) + \
			                     C7*(a[iy+3][ix+3][iz+3] + a[iy-2][ix+3][iz+3] + \
			                         a[iy-2][ix-2][iz+3] + a[iy+3][ix-2][iz+3] - \
			                         a[iy-2][ix-2][iz-2] - a[iy+3][ix-2][iz-2] - \
				                       a[iy+3][ix+3][iz-2] - a[iy-2][ix+3][iz-2]) +   \
			                     C6*(a[iy+2][ix+2][iz+2] + a[iy-1][ix+2][iz+2] + \
			                         a[iy-1][ix-1][iz+2] + a[iy+2][ix-1][iz+2] - \
			                         a[iy-1][ix-1][iz-1] - a[iy+2][ix-1][iz-1] - \
				                       a[iy+2][ix+2][iz-1] - a[iy-1][ix+2][iz-1]) +	\
			                     C5*(a[iy+1][ix+1][iz+1] + a[iy  ][ix+1][iz+1] + \
			                         a[iy  ][ix  ][iz+1] + a[iy+1][ix  ][iz+1] - \
			                         a[iy  ][ix  ][iz  ] - a[iy+1][ix  ][iz  ] - \
				                       a[iy+1][ix+1][iz  ] - a[iy  ][ix+1][iz  ])  )*s


enum AniType {ORTHORHOMBIC=0,TRICLINIC=1};
enum SourceType {STRESS=2, DISPLACEMENT=1, ACCELERATION=0, TENSOR=3};

int main(int argc, char* argv[])
{
    /* Declare RSF params */
    bool verb,fsrf,snap,dabc,debug,abcone,cfl,is2d,opot;
    int  jsnap,ntsnap,jdata,nbell;
    enum AniType type;
    enum SourceType srctype;
    float fmax,safety;
    

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    int     nt,nz,nx,ny=0,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy=0,idz,idx,idy;

    sf_axis at,ax,ay=NULL,az;
    sf_axis as,ar,ac,acq;

    /* init RSF */
    sf_init(argc,argv);
   
   
    int tsrc;
    if(! sf_getint("srctype",&tsrc)) tsrc=0; /* source type, see comments */
    
    srctype = tsrc;

    /* MAKE SURE USER HAS DECLARED ANISOTROPY TYPE*/
    int ttype;
    if(! sf_getint("ani", &ttype)) ttype=-1;/* Anisotropy type, see comments */
    if(ttype == -1 || ttype >2 || ttype < 0) sf_error("Invalid Anisotropy type");
    type = ttype;

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    float **dd=NULL;           /* data      */
    /*------------------------------------------------------------*/
    /* Temp values for strain */
    float    szz,   sxx,   syy,   sxy,   syz,   szx;
    /*------------------------------------------------------------*/
    /* execution flags */
    if(! sf_getbool("verb",&verb)) verb=true; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("dabc",&dabc)) dabc=true; /* use sponge absorbing BC */
    if(! sf_getbool("abcone",&abcone)) abcone=false; /* use sharp brake at end of boundary layer */
    if(! sf_getbool("debug",&debug)) debug=true; /* print debugging info */
    if(! sf_getbool("cfl",&cfl)) cfl=true; /* use CFL check, will cause program to fail if not satisfied */
    if(! sf_getbool("opot",&opot)) opot=false; /* output potentials */
	/*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    /* print debugging information */
    /*------------------------------------------------------------*/
    /* Determine if 3D or 2D, test for existence of 4th axis, n=1 if none */
    sf_axis test = sf_iaxa(Fccc,4);
    if (sf_n(test) == 1) is2d = true;
    else is2d = false;
    
    
    if(!dabc){
        fprintf(stderr,"*********************\n");
        fprintf(stderr,"NOT USING ABSORBING BOUNDARY CONDITIONS!!!\n");
        fprintf(stderr,"THIS WILL CAUSE LARGE REFLECTIONS AT EDGES\n");
        fprintf(stderr,"TO TURN ON ABSORBING, SET DABC=y\n");
        fprintf(stderr,"*********************\n");
    }

    /* Make sure that we have the proper number of coefs in our stiffness file */
    sf_axis temp; int tempn;

    if (is2d) {
        temp = sf_iaxa(Fccc,3);
    } else {
        temp = sf_iaxa(Fccc,4);
    }
    tempn = sf_n(temp);
    sf_warning("Stiffness coefficient file has %d coefficients\n",tempn);

    switch(type) {
        case ORTHORHOMBIC:
            if (is2d && tempn != 4) 
              sf_error("ISO/HTI/VTI requires 4 coefficients in 2D");
            else if ( !is2d && tempn != 9) 
              sf_error("ISO/HTI/VTI requires 9 coefficients in 3D");
            break;
        case TRICLINIC:
            if (is2d && tempn != 6) 
              sf_error("TTI requires 6 coefficients in 2D");
            else if ( !is2d && tempn != 21) 
              sf_error("TTI requires 21 coefficients in 3D");
            break;
    }
    /*-------------------------------------------------------------*/
    /* Wavefield cut parameters */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;

    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"time"); 
    	sf_setunit(at,"s"); if(verb) sf_raxa(at); /* time */

    az = sf_iaxa(Fccc,1); sf_setlabel(az,"space z"); 
   	 sf_setunit(az,"km");if(verb) sf_raxa(az); /* depth */

    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"space x");  
    	sf_setunit(ax,"km");if(verb) sf_raxa(ax); /* space x */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"sources");  
    	sf_setunit(as,"km");if(verb) sf_raxa(as); /* sources */

    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"receivers");
    	sf_setunit(ar,"km");if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    
    if(!is2d){ 
    	ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"space y");  
	  	sf_setunit(ay,"km");if(verb) sf_raxa(ay); /* space y */
    	ny = sf_n(ay); dy = sf_d(ay);
    }

    ns = sf_n(as);
    nr = sf_n(ar);

    if(snap){
      if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
      if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

      if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
      if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

      if(!is2d){
        if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);
        if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);
      }
    }
   /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;
    
    if( !sf_getint("nb",&nb)) nb=100; /* padding size for absorbing boundary */
    if (nb < NOP) nb = NOP;
    
    if(snap) {  /* save wavefield every *jsnap* time steps */
	    if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }
    /* CFL check parameters */
    if (cfl){
      if(! sf_getfloat("fmax",&fmax)) sf_error("CFL: you did not specify fmax");

      if (! sf_getfloat("safety",&safety) || safety < 0) {
        switch(type){
          case ORTHORHOMBIC:
            safety = 0.75f;
            break;
          case TRICLINIC:
            safety = 0.50f;
            break;
        }
      }
      sf_warning("CFL: safety margin %f",safety);
    }

    /* --------------2d code---------------- */
if (is2d){
    /* IO Arrays */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
    /* FDM structure */
    fdm2d    fdm=NULL;
    abcone2d abcp=NULL,abcs=NULL;
    sponge   spo=NULL;

    /*------------------------------------------------------------*/
    float **tt=NULL;
    float **ro=NULL;           /* density   */
    float **c11=NULL;
    float **c33=NULL;
    float **c55=NULL;
    float **c13=NULL;
    float **c15=NULL;
    float **c35=NULL;
    float **vp,**vs;
    float **qp=NULL,**qs=NULL;
	
    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **umz,**uoz,**upz,**uaz,**utz; 
    float **umx,**uox,**upx,**uax,**utx;
    /* stress/strain tensor */ 
    float **tzz,**tzx,**txx; 
    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices */
    lint2d cs,cr;
    /* wavefield cut params */
    float     **uc=NULL;
    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);
    fdbell_init(nbell);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/
    /* 2D vector components */
    nc=2;
    ac=sf_maxa(nc,0,1);
    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);
    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);
    /* setup output wavefield header */
    if(snap) {
      dqz=sf_d(az);
      dqx=sf_d(ax);
 
      acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
      acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
      /* TODO: check if the imaging window fits in the wavefield domain */
 
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
      sf_oaxa(Fwfl,ac, 3);
      sf_oaxa(Fwfl,at, 4);
    }

    /*------------------------------------------------------------*/
    /* source array */
    if(srctype == TENSOR) {
        ww=sf_floatalloc3(ns,3,nt); 
        sf_floatread(ww[0][0],nt*3*ns,Fwav);
    } else {
        ww=sf_floatalloc3(ns,nc,nt); 
        sf_floatread(ww[0][0],nt*nc*ns,Fwav);
    }

    /* data array */
    dd=sf_floatalloc2(nr,nc);
    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm); /* setup linear interpolation */
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */

    idz = 1/dz;
    idx = 1/dx;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 

    ro =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input density */
    sf_floatread(tt[0],nz*nx,Fden);     expand(tt,ro ,fdm);

    switch(type){
      case ORTHORHOMBIC:
        c11=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c33=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c55=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c13=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        /* input stiffness */
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,fdm);
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,fdm);
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,fdm);
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,fdm);    
      break;
      case TRICLINIC:
        c11=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c33=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c55=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c13=sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
        c15=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        c35=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
        /* input stiffness */
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,fdm);
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,fdm); 
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c15,fdm); 
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,fdm);
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c35,fdm);  
        sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,fdm);              
      break;
    }
    free(*tt); free(tt);
    /*------------------------------------------------------------*/
	/* one-way abc setup   */
   	vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
	  vs = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    float vpmax = 0.0f;
    float vpmin = 100000.0f;
    float vsmax = 0.0f;
    float vsmin = 100000.0f;
    float c11f;
    float c55f;
   	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for  (iz=0; iz<fdm->nzpad; iz++) {
        c11f = c11[ix][iz]; 
        c55f = c55[ix][iz];
        if(c11f < 0 ) {
          vp[ix][iz] = sqrt( -c11f/ro[ix][iz] );
        } else {
          vp[ix][iz] = sqrt( c11f/ro[ix][iz] );
        }
        if (vp[ix][iz] > vpmax) vpmax = vp[ix][iz];
        if (vp[ix][iz] < vpmin) vpmin = vp[ix][iz];
        if(c55f < 0 ) {
          vs[ix][iz] = sqrt(-c55f/ro[ix][iz] );
        } else {
          vs[ix][iz] = sqrt( c55f/ro[ix][iz] );
        }
        if( vs[ix][iz] > vsmax) vsmax = vs[ix][iz];
        if( vs[ix][iz] < vsmin) vsmin = vs[ix][iz];
	    }
    }
    if (cfl) {
      cfl_elastic( vpmin,vpmax,vsmin,vsmax,
                   dx,-1.0,dz,dt,fmax,safety,4);
    }
    /* Absorbing boundary conditions setup */
    if (abcone) {
      abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
      abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);
    }
	  free(*vp); free(vp);
	  free(*vs); free(vs);
	  /* sponge abc setup */
	  spo = sponge_make(fdm->nb);

    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */
    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
	        ro [ix][iz] = dt*dt/ro[ix][iz];   
	    }
    }
    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    umz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uoz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    upz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uaz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    umx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    upx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    tzz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    tzx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    txx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    /* Zero out values */

    for(  ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
        umz[ix][iz]=0; umx[ix][iz]=0;
        uoz[ix][iz]=0; uox[ix][iz]=0;
        upz[ix][iz]=0; upx[ix][iz]=0;
        uaz[ix][iz]=0; uax[ix][iz]=0;
        tzz[ix][iz]=0; tzx[ix][iz]=0; 
        txx[ix][iz]=0;
	    }
    }

    if(opot) {
    	qp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    	qs=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
    	if(verb) fprintf(stderr,"%d/%d \r",it,nt);

    	/*------------------------------------------------------------*/
    	/* from displacement to strain                                */
    	/*------------------------------------------------------------*/
    	/* 
    	 * exx = Fx(ux)
    	 * ezz = Fz(uz)
    	 * ezx = Bx(uz) + Bz(ux)
    	 */
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk)		\
    private(iz,ix)				\
    shared(fdm,tzz,tzx,txx,uoz,uox,idz,idx)
#endif
    	for  (ix=NOP; ix<fdm->nxpad-NOP-1; ix++) {
    	  for(iz=NOP; iz<fdm->nzpad-NOP-1; iz++) {
          txx[ix][iz] = Dx(uox,ix,iz,idx);
          tzz[ix][iz] = Dz(uoz,ix,iz,idz);
          tzx[ix][iz] = Dx(uoz,ix,iz,idx) + Dz(uox,ix,iz,idz);
    	  }
    	}
			
    	/*------------------------------------------------------------*/
    	/* from strain to stress                                      */
    	/*------------------------------------------------------------*/
      switch(type){
        case ORTHORHOMBIC:
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk)		\
    private(iz,ix,szz,szx,sxx)			\
    shared(fdm,tzz,tzx,txx,c11,c33,c55,c13)
#endif
        for  (ix=0; ix<fdm->nxpad; ix++) {
          for(iz=0; iz<fdm->nzpad; iz++) {

            sxx = c11[ix][iz] * txx[ix][iz]
                + c13[ix][iz] * tzz[ix][iz];
                    
            szz = c13[ix][iz] * txx[ix][iz]
                + c33[ix][iz] * tzz[ix][iz]; 
            
            szx = c55[ix][iz] * tzx[ix][iz];

            txx[ix][iz] = sxx;
            tzz[ix][iz] = szz;
            tzx[ix][iz] = szx;
          }
        }    
        break;
        case TRICLINIC:
#ifdef _OPENMP
#pragma omp parallel for	    \
    schedule(dynamic,fdm->ompchunk)		\
    private(iz,ix,szz,szx,sxx)			\
    shared(fdm,tzz,tzx,txx,c11,c33,c55,c13,c15,c35)
#endif
        for    (ix=0; ix<fdm->nxpad; ix++) {
          for(iz=0; iz<fdm->nzpad; iz++) {
            sxx = c11[ix][iz] * txx[ix][iz]
                + c13[ix][iz] * tzz[ix][iz]
                + c15[ix][iz] * tzx[ix][iz];
                    
            szz = c13[ix][iz] * txx[ix][iz]
                + c33[ix][iz] * tzz[ix][iz] 
                + c35[ix][iz] * tzx[ix][iz]; 
            
            szx = c15[ix][iz] * txx[ix][iz]
                + c35[ix][iz] * tzz[ix][iz]
                + c55[ix][iz] * tzx[ix][iz];

            txx[ix][iz] = sxx;
            tzz[ix][iz] = szz;
            tzx[ix][iz] = szx;
          }
        }    
        break;
      }

	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	/* Francesco: the z component of the traction must be zero at the free surface */
	
    	if(fsrf) {
    	  for(ix=0; ix<fdm->nxpad; ix++) {
    	    for(iz=0; iz<fdm->nb; iz++) {
    	      txx[ix][iz]=0;
    	      tzz[ix][iz]=0;
    	      tzx[ix][iz]=0;
    	    }
    	  }
    	 
    	  for(ix=0; ix<fdm->nxpad; ix++) {
        tzz[ix][iz]=0.0;  
    		} 
    	}

    	/*------------------------------------------------------------*/
    	/* inject stress source                                       */
    	/*------------------------------------------------------------*/
    	if(srctype == STRESS || srctype == TENSOR) {
  	    lint2d_bell(tzz,ww[it][0],cs);
  	    lint2d_bell(txx,ww[it][1],cs);
    	}
      if (srctype == TENSOR){
        lint2d_bell(tzx,ww[it][2],cs);
      }
      if(dabc){
        sponge2d_apply(txx,spo,fdm);
        sponge2d_apply(tzx,spo,fdm);
        sponge2d_apply(tzz,spo,fdm);    
    	}

    	/*------------------------------------------------------------*/
    	/* from stress to acceleration                                */
    	/*------------------------------------------------------------*/
    	/* 
    	 * ax = Bx(txx) + Fz(txz)
    	 * az = Fx(txz) + Bz(tzz)
    	 */
#ifdef _OPENMP
#pragma omp parallel for			\
    schedule(dynamic,fdm->ompchunk)		\
    private(iz,ix)				\
    shared(fdm,tzz,tzx,txx,uaz,uax,idz,idx)
#endif
    	for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
    	  for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
          uax[ix][iz] = Dx( txx,ix-1,iz-1,idx ) + Dz( tzx,ix-1,iz-1,idz );
          uaz[ix][iz] = Dx( tzx,ix-1,iz-1,idx ) + Dz( tzz,ix-1,iz-1,idz );
    	  }
    	}


    	/*------------------------------------------------------------*/
    	/* inject acceleration source                                 */
    	/*------------------------------------------------------------*/
    	if(srctype == ACCELERATION) {
        lint2d_bell(uaz,ww[it][0],cs); 
        lint2d_bell(uax,ww[it][1],cs);
    	}
  
    	/*------------------------------------------------------------*/
    	/* step forward in time                                       */
    	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for				\
    schedule(dynamic,fdm->ompchunk)			\
    private(iz,ix)					\
    shared(fdm,uoz,uox,umz,umx,upz,upx,uaz,uax,ro)
#endif
    	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
  	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
          upz[ix][iz] = 2*uoz[ix][iz] 
                      -   umz[ix][iz] 
                      +   uaz[ix][iz] * ro[ix][iz]; 

          upx[ix][iz] = 2*uox[ix][iz] 
                      -   umx[ix][iz] 
                      +   uax[ix][iz] * ro[ix][iz]; 
  	    }
    	}


    	if(srctype == DISPLACEMENT) {
        lint2d_bell(upz,ww[it][0],cs); 
        lint2d_bell(upx,ww[it][1],cs);
    	}
    	/*------------------------------------------------------------*/
    	/* apply the boundary condition                               */
    	/*------------------------------------------------------------*/
      for(ix = 0; ix < fdm->nxpad; ix++){
        for(iz = 0; iz < NOP; iz++){
          upz[ix][iz] = 0.0f;
          upx[ix][iz] = 0.0f;
        }
      }
      for(ix = 0; ix < fdm->nxpad; ix++){
        for(iz = fdm->nzpad-NOP; iz < fdm->nzpad; iz++){
          upz[ix][iz] = 0.0f;
          upx[ix][iz] = 0.0f;
        }
      }
      for(ix = 0; ix < NOP; ix++){
        for(iz = 0; iz < fdm->nzpad; iz++){
          upz[ix][iz] = 0.0f;
          upx[ix][iz] = 0.0f;
        }
      }
      for(ix = fdm->nxpad-NOP; ix < fdm->nxpad; ix++){
        for(iz = 0; iz < fdm->nzpad; iz++){
          upz[ix][iz] = 0.0f;
          upx[ix][iz] = 0.0f;
        }
      }
    	/* circulate wavefield arrays */
      /* Change pointers around */
    	utz=umz; utx=umx;
    	umz=uoz; umx=uox;
    	uoz=upz; uox=upx;
    	upz=utz; upx=utx;
  
  
      /* Apply the zero-incidence boundary condition */
      if(abcone){
        abcone2d_apply(uoz,umz,NOP,abcp,fdm);
        abcone2d_apply(uox,umx,NOP,abcp,fdm);
        
        abcone2d_apply(uoz,umz,NOP,abcs,fdm);
        abcone2d_apply(uox,umx,NOP,abcs,fdm);
      }
      /* sponge ABC */
      if(dabc){
        sponge2d_apply(umz,spo,fdm);
        sponge2d_apply(uoz,spo,fdm);
        sponge2d_apply(upz,spo,fdm);
        
        sponge2d_apply(umx,spo,fdm);
        sponge2d_apply(uox,spo,fdm);
        sponge2d_apply(upx,spo,fdm);
      }
    	/*------------------------------------------------------------*/
    	/* cut wavefield and save */
    	/*------------------------------------------------------------*/
    	lint2d_extract(uoz,dd[0],cr);
    	lint2d_extract(uox,dd[1],cr);

    	if (opot){
    	  if(snap && it%jsnap==0) {
        #ifdef _OPENMP
        #pragma omp parallel for			\
            schedule(dynamic,fdm->ompchunk)		\
            private(iz,ix)				\
            shared(fdm,uoz,uox,qp,qs,idz,idx)
        #endif
   
          for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
        		 for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
          
              qp[ix][iz] = Dz( uoz,ix,iz,idz ) 
                         + Dx( uox,ix,iz,idx );
          
              qs[ix][iz] = Dz( uox,ix,iz,idz ) 
                         - Dx( uoz,ix,iz,idx );
            }
          }
    	  	cut2d(qp,uc,fdm,acz,acx);
    	  	sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
    	  	
    	  	cut2d(qs,uc,fdm,acz,acx);
    	  	sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
        }	
    	  if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
    	    
    	}else{
    	  if(snap && it%jsnap==0) {
    	  	cut2d(uoz,uc,fdm,acz,acx);
    	  	sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
    	  	
    	  	cut2d(uox,uc,fdm,acz,acx);
    	  	sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
    	  }
    	  if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
    	}
    	/* Save snapshots of wavefield if selected */
    	/* Save data snapshots*/
    }
        
    if(verb) fprintf(stderr,"\n");      
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);
   
    free(*ro);  free(ro);
    /* Only free double pointers that are not null... */
    switch (type){
      case TRICLINIC:
        free(*c15); free(c15);
        free(*c35); free(c35); //flow through
      case ORTHORHOMBIC:
        free(*c11); free(c11);
        free(*c33); free(c33);
        free(*c55); free(c55);
        free(*c13); free(c13);
      break;
    }
/* Free all outstanding parameters */
    free(*umz); free(umz);
    free(*uoz); free(uoz);
    free(*upz); free(upz);
    free(*uaz); free(uaz);

    free(*umx); free(umx);
    free(*uox); free(uox);
    free(*upx); free(upx);
    free(*uax); free(uax);

    free(*tzz); free(tzz);
    free(*txx); free(txx);
    free(*tzx); free(tzx);

    if(snap){free(*uc);  free(uc);    }
  	if(opot) {
    	free(*qp);  free(qp);    
    	free(*qs);  free(qs);    
    }

    exit (0);
/* ------------end 2d ------------------ */
  } else {

/* ------------ 3d code ---------------*/  
    /* IO Arrays */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abcp=NULL,abcs=NULL;
    sponge   spo=NULL;
    /*------------------------------------------------------------*/
    float ***tt=NULL;
    float ***ro=NULL;           /* density */

    float ***c11=NULL,***c12=NULL,***c13=NULL;
    float ***c14=NULL,***c15=NULL,***c16=NULL;
    float ***c22=NULL,***c23=NULL,***c24=NULL,***c25=NULL,***c26=NULL;
    float ***c33=NULL,***c34=NULL,***c35=NULL,***c36=NULL;
    float ***c44=NULL,***c45=NULL,***c46=NULL;
    float ***c55=NULL,***c56=NULL;
    float ***c66=NULL;
  
    float ***vp,***vs;
    float ***qp=NULL;

    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float ***umz,***uoz,***upz,***uaz,***utz; 
    float ***umx,***uox,***upx,***uax,***utx;
    float ***umy,***uoy,***upy,***uay,***uty;
    /*------------------------------------------------------------*/

    /* stress/strain tensor */ 
    float ***tzz,***txx,***tyy,***txy,***tyz,***tzx;       

    /* linear interpolation weights/indices */
    lint3d cs,cr;

    /* wavefield cut params */
    float     ***uc=NULL;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    fdbell3d_init(nbell);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); 
    sf_setlabel(az,"expanded z");sf_setunit(az,"km");if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); 
    sf_setlabel(ax,"expanded x");sf_setunit(ax,"km");if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); 
    sf_setlabel(ay,"expanded y");sf_setunit(ay,"km");if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

    /* 3D vector components */
    nc=3;
    ac=sf_maxa(nc,0,1);

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header */
    if(debug) fprintf(stderr, "Setting up output wavefield header\n");
    if(snap) {
    	dqz=sf_d(az);
    	dqx=sf_d(ax);
    	dqy=sf_d(ay);
    
    	acz = sf_maxa(nqz,oqz,dqz);  sf_setunit(acz,"km");
      sf_setlabel(acz,"snapshot z");sf_raxa(acz);
    	acx = sf_maxa(nqx,oqx,dqx);  sf_setunit(acx,"km");
      sf_setlabel(acx,"snapshot x");sf_raxa(acx);
    	acy = sf_maxa(nqy,oqy,dqy);  sf_setunit(acy,"km");
      sf_setlabel(acy,"snapshot y");sf_raxa(acy);
    	/* check if the imaging window fits in the wavefield domain */
    
    	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));
    
    	ntsnap=0;
    	for(it=0; it<nt; it++) {
    	  if(it%jsnap==0) ntsnap++;
    	}
    	sf_setn(at,  ntsnap);
    	sf_setd(at,dt*jsnap);
    	sf_setlabel(at,"snapshot frames");
    	if(verb)  sf_raxa(at);
      
      if(opot) acq=sf_maxa(4,0,1);
      else acq = sf_maxa(nc,0,1);
    	sf_oaxa(Fwfl,acz,1);
    	sf_oaxa(Fwfl,acx,2);
    	sf_oaxa(Fwfl,acy,3);
    	sf_oaxa(Fwfl,acq,4);
    	sf_oaxa(Fwfl,at, 5);
    }

    /*------------------------------------------------------------*/
    /* source array */
    if (srctype == TENSOR) {
      ww = sf_floatalloc3(ns,6,nt);
      sf_floatread(ww[0][0],nt*6*ns,Fwav);
    }
    else {
      ww=sf_floatalloc3(ns,nc,nt); 
      sf_floatread(ww[0][0],nt*nc*ns,Fwav);
    }

    /* data array */
    dd=sf_floatalloc2(nr,nc);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */
    
    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 
    ro =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* input density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden);     expand3d(tt,ro ,fdm);

    if(debug) fprintf(stderr, "Beginning to read into stiffness arrays\n");

	/* allocate and read into arrays for stiffnesses */
    switch(type){
    	case ORTHORHOMBIC:
	  	c11=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c22=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c33=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c44=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c55=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c66=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c12=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c13=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c23=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);  

	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,fdm);

	    break;

    	case TRICLINIC:
	  	c11=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c12=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c13=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c14=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c15=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c16=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c22=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c23=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c24=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);    
	  	c25=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c26=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c33=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c34=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c35=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c36=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c44=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c45=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c46=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);    
	  	c55=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c56=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
	  	c66=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);  

	  	 /* input stiffness */
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c14,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c15,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c16,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c24,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c25,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c26,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c34,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c35,fdm);    
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c36,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c45,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c46,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c56,fdm);
	  	sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,fdm);
	    break;
    }
    free(**tt); free(*tt); free(tt);
    if(debug) 
      fprintf(stderr,"Successfully allocated and read in coefficients...\n");

    /*------------------------------------------------------------*/
    if(dabc) {
    	/* one-way abc setup   */
    	vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    	vs = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    
      float vpmin = 100000000;
      float vpmax = 0.0f;
      float vsmin = 100000000;
      float vsmax = 0.0f;
    	for (iy=0; iy<fdm->nypad; iy++) {
  	    for (ix=0; ix<fdm->nxpad; ix++) {
  		    for(iz=0; iz<fdm->nzpad; iz++) {
            if (c11[iy][ix][iz] < 0) {
              vp[iy][ix][iz] = sqrt(-c11[iy][ix][iz]/ro[iy][ix][iz]);
            } else {
              vp[iy][ix][iz] = sqrt(c11[iy][ix][iz]/ro[iy][ix][iz]);
            } 
            float vpt = vp[iy][ix][iz];
            if (vpt > vpmax) vpmax = vpt;
            if (vpt < vpmin) vpmin = vpt;
            if (c55[iy][ix][iz] < 0){
              vs[iy][ix][iz] = sqrt(-c55[iy][ix][iz]/ro[iy][ix][iz]);
            } else {
              vs[iy][ix][iz] = sqrt(c55[iy][ix][iz]/ro[iy][ix][iz]);
            }
            float vst = vs[iy][ix][iz];
            if (vst > vsmax) vsmax = vst;
            if (vst < vsmin) vsmin = vst;
  		    }
  	    }
    	}

      if (cfl) {
        cfl_elastic(vpmin,vpmax,vsmin,vsmax,
            dx,dy,dz,dt,fmax,safety,4);
      }
      if (abcone){
        abcp = abcone3d_make(NOP,dt,vp,fsrf,fdm);
        abcs = abcone3d_make(NOP,dt,vs,fsrf,fdm);
      }
    	free(**vp); free(*vp); free(vp);
    	free(**vs); free(*vs); free(vs);
    	/* sponge abc setup */
    	spo = sponge_make(fdm->nb);
    }
    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */    
    for (iy=0; iy<fdm->nypad; iy++) {
	    for  (ix=0; ix<fdm->nxpad; ix++) {
	      for(iz=0; iz<fdm->nzpad; iz++) {
	  	    ro[iy][ix][iz] = dt*dt/ro[iy][ix][iz];
	      }
	    }
    }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    if(debug) fprintf(stderr,"Allocating wavefield arrays\n");
    umz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uaz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uax=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uay=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    tzz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    txy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tyz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    tzx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    for (iy=0; iy<fdm->nypad; iy++) {
    	for (ix=0; ix<fdm->nxpad; ix++) {
	      for(iz=0; iz<fdm->nzpad; iz++) {
          umz[iy][ix][iz]=0; umx[iy][ix][iz]=0; umy[iy][ix][iz]=0;
          uoz[iy][ix][iz]=0; uox[iy][ix][iz]=0; uoy[iy][ix][iz]=0;
          upz[iy][ix][iz]=0; upx[iy][ix][iz]=0; upy[iy][ix][iz]=0;
          uaz[iy][ix][iz]=0; uax[iy][ix][iz]=0; uay[iy][ix][iz]=0;
	      }
	    }
    }
    if(opot) {
    	qp=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    }
    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");

    for (it=0; it<nt; it++) {
	    if(verb) fprintf(stderr,"%d/%d \r",it,nt);

    	/*------------------------------------------------------------*/
    	/* from displacement to strain                                */
    	/*------------------------------------------------------------*/	
    	/* 
    	 * exx = Dx(ux)
    	 * eyy = Dy(uy)
    	 * ezz = Dz(uz)
    	 * exy = Dy(ux) + Dx(uy)
    	 * eyz = Dz(uy) + Dy(uz)
    	 * ezx = Dx(uz) + Dz(ux)
    	 */
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,idx,idy,idz)
#endif
    	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
   	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
      		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) { 
      		  txx[iy][ix][iz] = Dx3(uox,ix,iy,iz,idx);
      		  tyy[iy][ix][iz] = Dy3(uoy,ix,iy,iz,idy);
      		  tzz[iy][ix][iz] = Dz3(uoz,ix,iy,iz,idz);
      		  
      		  txy[iy][ix][iz] = Dy3(uox,ix,iy,iz,idy) + Dx3(uoy,ix,iy,iz,idx);
      		  tyz[iy][ix][iz] = Dz3(uoy,ix,iy,iz,idz) + Dy3(uoz,ix,iy,iz,idy);
      		  tzx[iy][ix][iz] = Dx3(uoz,ix,iy,iz,idx) + Dz3(uox,ix,iy,iz,idz);
      		}
        }
    	}
	
    	/*------------------------------------------------------------*/
    	/* from strain to stress                                      */
    	/*------------------------------------------------------------*/
      if(debug) fprintf(stderr,"Going from strain to stress\n");
      switch (type){
      
        case TRICLINIC:
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
#endif
      	for        (iy=0; iy<fdm->nypad; iy++) {
     	    for    (ix=0; ix<fdm->nxpad; ix++) {
        		for(iz=0; iz<fdm->nzpad; iz++) {
        		    
        	    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c12[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c13[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c14[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c15[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c16[iy][ix][iz] * txy[iy][ix][iz];
        
        	    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c22[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c23[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c24[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c25[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c26[iy][ix][iz] * txy[iy][ix][iz];
        
        	    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c23[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c33[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c34[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c35[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c36[iy][ix][iz] * txy[iy][ix][iz];
        
        	    syz = c14[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c24[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c34[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c44[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c45[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c46[iy][ix][iz] * txy[iy][ix][iz];
        	    
        	    szx = c15[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c25[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c35[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c45[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c55[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c56[iy][ix][iz] * txy[iy][ix][iz];
        
        	    sxy = c16[iy][ix][iz] * txx[iy][ix][iz]
        	    		+ c26[iy][ix][iz] * tyy[iy][ix][iz]
        	    		+ c36[iy][ix][iz] * tzz[iy][ix][iz]
        	    		+ c46[iy][ix][iz] * tyz[iy][ix][iz]
        	    		+ c56[iy][ix][iz] * tzx[iy][ix][iz]
        	    		+ c66[iy][ix][iz] * txy[iy][ix][iz];
        
        	    txx[iy][ix][iz] = sxx;
        	    tyy[iy][ix][iz] = syy;
        	    tzz[iy][ix][iz] = szz;
        	    txy[iy][ix][iz] = sxy;
        	    tyz[iy][ix][iz] = syz;
        	    tzx[iy][ix][iz] = szx;
        	  }
      	  }
      	}
        break;
  
        case ORTHORHOMBIC:
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,c11,c22,c33,c44,c55,c66,c12,c13,c23)
#endif
      	for        (iy=0; iy<fdm->nypad; iy++) {
     	    for    (ix=0; ix<fdm->nxpad; ix++) {
        		for(iz=0; iz<fdm->nzpad; iz++) {
        		    
       		    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
            			+ c12[iy][ix][iz] * tyy[iy][ix][iz]
            			+ c13[iy][ix][iz] * tzz[iy][ix][iz];
            		    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
            			+ c22[iy][ix][iz] * tyy[iy][ix][iz]
            			+ c23[iy][ix][iz] * tzz[iy][ix][iz];
            		    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
            			+ c23[iy][ix][iz] * tyy[iy][ix][iz]
            			+ c33[iy][ix][iz] * tzz[iy][ix][iz];
        		    
      		    sxy = c66[iy][ix][iz] * txy[iy][ix][iz];
      		    syz = c44[iy][ix][iz] * tyz[iy][ix][iz];
      		    szx = c55[iy][ix][iz] * tzx[iy][ix][iz];
      		    
      		    txx[iy][ix][iz] = sxx;
      		    tyy[iy][ix][iz] = syy;
      		    tzz[iy][ix][iz] = szz;
      
      		    txy[iy][ix][iz] = sxy;
      		    tyz[iy][ix][iz] = syz;
      		    tzx[iy][ix][iz] = szx;
        		}
     	    }
      	}
        break; 
      }
    /* End Anisotropy LOGIC */
      if(dabc) {
  	    /* sponge ABC */
  	    sponge3d_apply(txx,spo,fdm);
  	    sponge3d_apply(tzx,spo,fdm);
  	    sponge3d_apply(tyz,spo,fdm);
  	    
  	    sponge3d_apply(tzz,spo,fdm);
  	    sponge3d_apply(tyy,spo,fdm);
  	    sponge3d_apply(txy,spo,fdm);
    	}
    	/*------------------------------------------------------------*/
    	/* free surface */
    	/*------------------------------------------------------------*/
    	/* Francesco: the z component of the traction must be
         zero at the free surface */
    	if(fsrf) {
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx)
        for(iy=0; iy<fdm->nypad; iy++) {
    	    for(ix=0; ix<fdm->nxpad; ix++) {
            for(iz=0; iz < fdm->nb; iz++) {
              txx[iy][ix][iz]=0;
              tyy[iy][ix][iz]=0;
              tzz[iy][ix][iz]=0;
              tyz[iy][ix][iz]=0;
              tzx[iy][ix][iz]=0;
              txy[iy][ix][iz]=0;
            }
    	    }
        }
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx)
    	  for(iy=0; iy<fdm->nypad; iy++) {
    		  for(ix=0; ix<fdm->nxpad; ix++) {
    		    for(iz=fdm->nb; iz < fdm->nb+1; iz++) {
              tzz[iy][ix][iz]=0;
              tyz[iy][ix][iz]=0;
            }
          }
        }
    	}

    	/*------------------------------------------------------------*/
    	/* inject stress source                                       */
    	/*------------------------------------------------------------*/
    	if(srctype == STRESS || srctype == TENSOR) {
    	  lint3d_bell(txx,ww[it][0],cs);
    	  lint3d_bell(tyy,ww[it][1],cs);
    	  lint3d_bell(tzz,ww[it][2],cs);
    	}
      if (srctype == TENSOR) {
        lint3d_bell(tyz,ww[it][3],cs);
        lint3d_bell(tzx,ww[it][4],cs);
        lint3d_bell(txy,ww[it][5],cs);
      }
  
    	/*------------------------------------------------------------*/
    	/* from stress to acceleration                                */
    	/*------------------------------------------------------------*/
    	/* 
    	 * ax = Dx(txx) + Dy(txy) + Dz(txz)
    	 * ay = Dx(txy) + Dy(tyy) + Dz(tyz)
    	 * az = Dx(txz) + Dy(tyz) + Dz(tzz)
    	 */	
      if(debug) fprintf(stderr,"Going from stress to acceleration \n");
#ifdef _OPENMP
#pragma omp parallel for					\
    schedule(dynamic,fdm->ompchunk)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,idx,idy,idz)
#endif
    	for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
    	  for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
    		  for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		    
    		    uax[iy][ix][iz] = Dx3( txx,ix-1,iy-1,iz-1,idx) 
                            + Dy3( txy,ix-1,iy-1,iz-1,idy) 
                            + Dz3( tzx,ix-1,iy-1,iz-1,idz) ;
    		    uay[iy][ix][iz] = Dx3( txy,ix-1,iy-1,iz-1,idx) 
                            + Dy3( tyy,ix-1,iy-1,iz-1,idy) 
                            + Dz3( tyz,ix-1,iy-1,iz-1,idz) ;
    		    uaz[iy][ix][iz] = Dx3( tzx,ix-1,iy-1,iz-1,idx) 
                            + Dy3( tyz,ix-1,iy-1,iz-1,idy) 
                            + Dz3( tzz,ix-1,iy-1,iz-1,idz) ;		    
    		  }
    	  }
    	}
    	/*------------------------------------------------------------*/
    	/* inject acceleration source                                 */
    	/*------------------------------------------------------------*/
    	if(srctype == ACCELERATION) {
    	  lint3d_bell(uaz,ww[it][0],cs);
    	  lint3d_bell(uax,ww[it][1],cs);
    	  lint3d_bell(uay,ww[it][2],cs);
    	}
    
    	/*------------------------------------------------------------*/
    	/* step forward in time                                       */
    	/*------------------------------------------------------------*/
      if(debug) fprintf(stderr,"Trying to step forward in time\n");
#ifdef _OPENMP
#pragma omp parallel for						\
    schedule(dynamic,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(fdm,uox,uoy,uoz,umx,umy,umz,upx,upy,upz,uax,uay,uaz,ro)
#endif
    	for        (iy=0; iy<fdm->nypad; iy++) {
    	  for    (ix=0; ix<fdm->nxpad; ix++) {
    		  for(iz=0; iz<fdm->nzpad; iz++) {
    		    upx[iy][ix][iz] = 2*uox[iy][ix][iz] 
    		                  	-   umx[iy][ix][iz] 
    		                  	+   uax[iy][ix][iz] * ro[iy][ix][iz]; 
    
    		    upy[iy][ix][iz] = 2*uoy[iy][ix][iz] 
    		                  	-   umy[iy][ix][iz] 
    		                  	+   uay[iy][ix][iz] * ro[iy][ix][iz]; 
    
    		    upz[iy][ix][iz] = 2*uoz[iy][ix][iz] 
    	                   		-   umz[iy][ix][iz] 
    	                   		+   uaz[iy][ix][iz] * ro[iy][ix][iz]; 
    		    
    		  }
    	  }
    	}

    	if(srctype == DISPLACEMENT) {
    	  lint3d_bell(upz,ww[it][0],cs);
    	  lint3d_bell(upx,ww[it][1],cs);
    	  lint3d_bell(upy,ww[it][2],cs);
    	}
  

      if(debug) fprintf(stderr,"Done stepping forward\n");
    	/* circulate wavefield arrays */
    	utz=umz; uty=umy; utx=umx;
    	umz=uoz; umy=uoy; umx=uox;
    	uoz=upz; uoy=upy; uox=upx;
    	upz=utz; upy=uty; upx=utx;
    	
    	/* one-way ABC */
      if(abcone){
        abcone3d_apply(uoz,umz,NOP,abcp,fdm);
        abcone3d_apply(uox,umx,NOP,abcp,fdm);
        abcone3d_apply(uoy,umy,NOP,abcp,fdm);
        
        abcone3d_apply(uoz,umz,NOP,abcs,fdm);
        abcone3d_apply(uox,umx,NOP,abcs,fdm);
        abcone3d_apply(uoy,umy,NOP,abcs,fdm);
      }
  
    	if(dabc) {
    	  /* sponge ABC */
    	  sponge3d_apply(umz,spo,fdm);
    	  sponge3d_apply(uoz,spo,fdm);
    	  sponge3d_apply(upz,spo,fdm);
    	  
    	  sponge3d_apply(umx,spo,fdm);
    	  sponge3d_apply(uox,spo,fdm);
    	  sponge3d_apply(upx,spo,fdm);
    
    	  sponge3d_apply(umy,spo,fdm);
    	  sponge3d_apply(uoy,spo,fdm);
    	  sponge3d_apply(upy,spo,fdm);
    	}	  
  
    	/*------------------------------------------------------------*/
    	/* cut wavefield and save */
    	/*------------------------------------------------------------*/
      if(debug) fprintf(stderr,"Trying to extract data\n");
      lint3d_extract(uoz,dd[0],cr);
      lint3d_extract(uox,dd[1],cr);
      lint3d_extract(uoy,dd[2],cr);
      if(debug) fprintf(stderr,"Trying to save snapshots\n");
      
      if(opot){
        if(snap && it%jsnap==0) {
          for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
            for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
         	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
    	          qp[iy][ix][iz] = Dz3( uoz,ix,iy,iz,idz ) 
    	                         + Dy3( uoy,ix,iy,iz,idy ) 
    	                         + Dx3( uox,ix,iz,iz,idx );
              }
    	    	}
    	    }
          cut3d(qp,uc,fdm,acz,acx,acy);
          sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
      
          for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
            for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
         		  for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
    	          qp[iy][ix][iz] = Dy3( uoz,ix,iy,iz,idy ) 
    	                         - Dz3( uoy,ix,iy,iz,idz );
    	    	  }
            }
    	    }
          cut3d(qp,uc,fdm,acz,acx,acy);
          sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
      
          for    (iy=NOP; iy<fdm->nypad-NOP; iy++) {
            for  (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
         		  for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
    	          qp[iy][ix][iz] = Dz3( uox,ix,iy,iz,idz ) 
    	                         - Dx3( uoz,ix,iy,iz,idx );
              }
    	    	}
    	    }
          cut3d(qp,uc,fdm,acz,acx,acy);
          sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
          for     (iy=NOP; iy<fdm->nypad-NOP; iy++) {
            for   (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
         		  for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {
    	          qp[iy][ix][iz] = Dx3( uoy,ix,iy,iz,idx ) 
    	                         - Dy3( uox,ix,iy,iz,idy );
              }
    	    	}
    	    }
          cut3d(qp,uc,fdm,acz,acx,acy);
          sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
        }
        if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
      }else{
        if(snap && it%jsnap==0) {
    	    cut3d(uoz,uc,fdm,acz,acx,acy);
    	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
    
    	    cut3d(uox,uc,fdm,acz,acx,acy);
    	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
    
    	    cut3d(uoy,uc,fdm,acz,acx,acy);
    	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
        }
        if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
      }
    }
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    if(debug) fprintf(stderr,"Finished loop, trying to deallocate...\n"); 
    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);
    if(debug) fprintf(stderr,"Deallocating coefficient matrices...\n");
  
    free(**ro);  free(*ro);  free(ro);
    switch(type){
      case TRICLINIC:
        free(**c14); free(*c14); free(c14);
        free(**c15); free(*c15); free(c15);
        free(**c16); free(*c16); free(c16);
        free(**c24); free(*c24); free(c24);
        free(**c25); free(*c25); free(c25);
        free(**c26); free(*c26); free(c26);
        free(**c34); free(*c34); free(c34);
        free(**c35); free(*c35); free(c35);
        free(**c36); free(*c36); free(c36);
        free(**c45); free(*c45); free(c45);
        free(**c46); free(*c46); free(c46);
        free(**c56); free(*c56); free(c56);
      case ORTHORHOMBIC:
        free(**c11); free(*c11); free(c11);
        free(**c22); free(*c22); free(c22);
        free(**c33); free(*c33); free(c33);
        free(**c44); free(*c44); free(c44);
        free(**c55); free(*c55); free(c55);
        free(**c66); free(*c66); free(c66);
        free(**c12); free(*c12); free(c12);
        free(**c13); free(*c13); free(c13);
        free(**c23); free(*c23); free(c23);
      break;
    }
  
    free(**umz); free(*umz); free(umz);
    free(**uoz); free(*uoz); free(uoz);
    free(**upz); free(*upz); free(upz);
    free(**uaz); free(*uaz); free(uaz);

    free(**umx); free(*umx); free(umx);
    free(**uox); free(*uox); free(uox);
    free(**upx); free(*upx); free(upx);
    free(**uax); free(*uax); free(uax);

    free(**umy); free(*umy); free(umy);
    free(**uoy); free(*uoy); free(uoy);
    free(**upy); free(*upy); free(upy);
    free(**uay); free(*uay); free(uay);

    free(**tzz); free(*tzz); free(tzz);
    free(**txx); free(*txx); free(txx);
    free(**tyy); free(*tyy); free(tyy);
    free(**txy); free(*txy); free(txy);
    free(**tyz); free(*tyz); free(tyz);
    free(**tzx); free(*tzx); free(tzx);
    if(snap) {free(**uc);  free(*uc);  free(uc);}
    if(opot) {free(**qp);  free(*qp);  free(qp);}
    exit (0);
  }
}
/* ---------------- end 3d code --------------*/
