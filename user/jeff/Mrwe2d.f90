!  Riemannian Wavefield Extrapolation for SHOT-PROFILE migration in  2D generalized coordinates.
!
!		The key concept is that wave-equation migration can be performed
!		on almost any type of grid called Riemannian coordintes
!		This program allows for user to input 
!		generalized coordinate and source and receiver wavefields that are recorded on the 
!		zeroth extrapolation step, and then propagate them into the user-specified
!		computational mesh.


!!$  Copyright (C) 2007 Stanford University
!!$  
!!$  This program is free software; you can redistribute it and/or modify
!!$  it under the terms of the GNU General Public License as published by
!!$  the Free Software Foundation; either version 2 of the License, or
!!$  (at your option) any later version.
!!$  
!!$  This program is distributed in the hope that it will be useful,
!!$  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!$  GNU General Public License for more details.
!!$  
!!$  You should have received a copy of the GNU General Public License
!!$  along with this program; if not, write to the Free Software
!!$  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
program Mrwe2D   
  !! V3VX-9BHSTSNC 
  !! . . Riemannian Wavefield Extrapolation program for 
  !! . . SHOT-PROFILE migration 
  !! . . in 2D generalized coordinates. 
  !!
  !! . . General interal processing steps
  !! 		1) Input wavefields, coordinate system, Velocity model (CARTESIAN)
  !!		2) Calculate geometrical coefficients for coordinate system and velocity model
  !!		3) Propagate wavefields and apply imaging condition
  !!		4) Output Image as: a) interpolated to Cartesian; and b) In original coordinates
  !!
  !!
  !!		Requisite files:
  !!
  !!		1) "rwf" tag is the input receiver wavefield in nx-nw-ns format
  !!		2) "swf" tag is the input source wavefield in nx-nw-ns format
  !! 		3) "vel" tag is the input velocity model in nz-nx format
  !!		4) "rays" tag is the input coordinate system in nz-nx-2 format
  !!				where nz-nx-1 is the x-coordinate field and nz-
  !!                              nx-2 the z-coordinate field
  !!		5) "image" tag is the Cartesian output imagein cnz-cnx-ns format
  !!		6) "Rimage" tag is the Ray based ouput image in nz-nx-ns formate
  !!
  !!		Required parameters:
  !!
  !!		1) nref=256 is the starting number of points for calculating reference velocities
  !!		2) verbose=1 is the level of verbosity for the RWE2D migration code
  !!		3) min_region=.1 controls how close a true coefficient is to a reference value 
  !!				before it is thrown away (in %)
  !!		4) del_dist=0.05 is the minimum distance between reference parameters (in %)
  !!		5) verb=1 is the verbage level of the Lloyds code
  !!		6) niter_lloyd=25 the number of iterations of Lloyds code
  !!		7) nloop=4 number of smoothing iterations for merging wavefields 
  !!                      from different references
  !!		8) migapp=1024 migration aperature split evenly about shotpoint
  !!		9) forward=1 Forward/Backscattering switch (0/1)
  !!		10) norm=1 Normalize by Jacobian of matrix
  !!		11/12) zmin/zmax min/max output Cartesian image z-locations
  !!		13/14) xmin/xmax min/max output Cartesian image x-locations
  !!		15/16) dzz/dxx grid spacing in output Cartesian image
  !!
  !! 		Written by Jeff Shragge
  !!				May 15/2006
  !!				contact: jshragge@gmail.com
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Required modules
  use rsf			!! Basic SEPlib I/O 
  use util			!! handles axis I/O and aux file I/O
  use RWE2D_wemig		!! Migration driver module
  use C2R2D			!! Ray-to-Cartesian interpolation 
  use geometry2D !! Calculate geometry coefficients for dispersion relation

  implicit none
  type(myaxis) :: ax,az,aw,ar,as 	!! axes types 

  integer :: ix,iz,ii,iw,ishot(1),iis,mxl,mxr,oldn
  integer :: nref,nsx,nsz,norm,cnx,cnz,migapp
  real    :: forward,xmin,zmin,xmax,zmax,dxx,dzz
  real,    allocatable,dimension(:,:,:):: raysin     !! Input rays		
  real,    allocatable,dimension(:,:)  :: Vel,gnorm  !! Velocity model, norm matrix
  complex, pointer    ,dimension(:,:  ):: rwf, swf   !! receiver and source wavefiels
  real,    allocatable,dimension(:,:)  :: CCimg, Pimg!! Cartesian/ray coordinate images
  logical                              :: verbose,kin
  character(len=128) :: name
  type (file) :: rwfile, rays, image, Rimage, velfile, swfile

  !! Geometry calculations
  integer, allocatable :: mask(:,:),cutmask(:,:)
  real,    allocatable :: refs(:,:,:), fields(:,:,:),cutfields(:,:,:)

  call sf_init()
  rwfile = rsf_input("rwf")
  swfile = rsf_input("swf")
  velfile = rsf_input("vel")
  rays = rsf_input("rays")
  image = rsf_output("image")
  Rimage = rsf_output("Rimage")

  !!		Routines to enable parallel processing with
  !!		SEPlib "Parallel" command
!  call getch_add_string("head=/dev/null")
!  call getch_add_string("noheader=y")
!  call set_no_putch()
!  call sep_begin_prog()

  name="rays"

  !! 		Read in input dimensions from the "rwf" tag
  !!		Determines the x grid parameters
  call from_par(rwfile,"n1",ax%n)
  call from_par(rwfile,"d1",ax%d)
  call from_par(rwfile,"o1",ax%o)
  !! 		Determine the frequencies parameters
  call from_par(rwfile,"n2",aw%n)
  call from_par(rwfile,"o2",aw%o)
  call from_par(rwfile,"d2",aw%d)
  !!		Determine the number of shot parameters
  call from_par(rwfile,"n3",as%n)
  call from_par(rwfile,"d3",as%d)
  call from_par(rwfile,"o3",as%o)

  !! 		Read in input dimensions from the "rays" tag
  !!		Determines the grid for extrapolation direction  
  call from_par(rays,"n1",az%n)
  call from_par(rays,"d1",az%d)
  call from_par(rays,"o1",az%o)

  !! 		Parameter input section.  
  !!	 	Call from command line, Makefile or a "par=pars.P" file
  call from_par("migapp",migapp,ax%n)
  call from_par("forward",forward,0.) 	! Forward scattering option
  call from_par("nref",nref,256)        	! starting number of points for calculating reference velocities
  call from_par("verbose",verbose,.false.)	! level of verbosity
  call from_par("kinematic",  kin,.true.)	! Kinematic approximation 
  call from_par("norm",norm,1)		! Whether (1) or not (0) to normalize by gnorm


  write(0,*) 'DEPTH',az%n,az%o,az%d
  write(0,*) 'X1   ',ax%n,ax%o,ax%d
  write(0,*) 'SHOTS',as%n,as%o,as%d
  write(0,*) 'FREQ ',aw%n,aw%o,aw%d

  !! Define output Cartesian grid
  call from_par("xmin",xmin,ax%o)
  call from_par("zmin",zmin,az%o)
  call from_par("xmax",xmax,(ax%n-1)*ax%d+ax%o)
  call from_par("zmax",zmax,(az%n-1)*az%d+az%o)
  call from_par("dxx",dxx,ax%d)
  call from_par("dzz",dzz,az%d)
  call from_par("nsx",nsx,3)
  call from_par("nsz",nsz,3)
  cnx = floor((xmax - xmin)/dxx)+1
  cnz = floor((zmax - zmin)/dzz)+1

  !! Create output header file associated with "image" tag
  !! . . Output Image in Cartesian coordinates
  call to_par(image,"n1",cnz)
  call to_par(image,"n2",cnx)
  call to_par(image,"n3",as%n)
  call to_par(image,"d1",dzz )
  call to_par(image,"d2",dxx )
  call to_par(image,"d3",as%d)
  call to_par(image,"o1",zmin)
  call to_par(image,"o2",xmin)
  call to_par(image,"o3",as%o)
  call settype(image,sf_float)

  !! Create output header file associated with "image" tag
  !! . . Output Image in RWE coordinates
  call to_par(Rimage,"n1",az%n)
  call to_par(Rimage,"n2",ax%n)
  call to_par(Rimage,"n3",as%n)
  call to_par(Rimage,"d1",az%d)
  call to_par(Rimage,"d2",ax%d)
  call to_par(Rimage,"d3",as%d)
  call to_par(Rimage,"o1",az%o)
  call to_par(Rimage,"o2",ax%o)
  call to_par(Rimage,"o3",as%o)
  call settype(Rimage,sf_float)

  !!	Allocate required arrays
  allocate( rwf(migapp,aw%n), swf(migapp,aw%n) )
  allocate( raysin(az%n,ax%n,2),Vel(az%n,ax%n),gnorm(az%n,ax%n) )
  allocate( mask(ax%n  ,az%n), refs(nref,4,az%n),fields(ax%n,4,az%n))
  allocate( cutmask(migapp,az%n), cutfields(migapp,4,az%n) )
  allocate( CCimg(cnz,cnx) , Pimg(az%n,migapp) )

  !! Zero all arrays
  rwf=0.;  	swf=0.; 	gnorm = 0.; raysin=0
  Vel=0.;	mask=0;	refs=0.;	fields=0.
  cutmask=0.;	cutfields=0.; CCimg=0.; Pimg=0.

  !!	Read in rays file
  call rsf_read(rays,raysin)

  !! Read in Velocity file
  call rsf_read(velfile,Vel)

  write(0,*) 'VELOCITY MIN/MAX', minval(Vel),maxval(Vel)


  !!----------------------------------------------------------------------------------
  !!
  !! GEOMETRY SECTION
  !!
  !! . . Initialize module file

  call geometry_init(az%n,ax%n,az%d,ax%d,kin,nref,verbose)

  !! . . Calculate the appropriate coefficients
  !! . . INPUT:
  !!		1) raysin - the coordinate mesh used to image
  !!		2) Vel - velocity model
  !!
  !!	. . OUTPUT:
  !!		1) fields - the coefficient fields for entire coordinate system (used in SSF)
  !!		2) refs - a sparse representation of 'fields"' used in the PSPI program
  !!		3) mask - an integer matrix listing the closest refs to each field value
  !!		4) gnorm - weight matrix for normalization 

  call geometry_calc(raysin,Vel,fields,refs,mask,gnorm)
  call geometry_close()
  deallocate( Vel )
  write(0,*) 'GNORM',minval(gnorm),maxval(gnorm)

  !!-------------------------------------------------------------------------------
  !!
  !! MIGRATION SECTION
  !!
  !! . . Initialize module and image space
  oldn=ax%n
  ax%n=migapp
  !! . . Initialize the wave-equation migration driver
  call wemig_init(ax,az,aw,ar,forward,kin,nref,verbose)

  !! Initialize the RWE to Cartesian interpolation routine
  call C2R_init(nsx,nsz,norm,dxx,dzz,cnx,cnz,migapp,az%n,xmin,zmin)

  do ii=1,as%n   !! . . Loop over shot points
     swf  = 0.; rwf   = 0.
     call rsf_read(rwfile,rwf)  !! Read in only current receiver wavefield
     call rsf_read(swfile,swf)	 !! Read in only currrent source wavefield
     Pimg = 0.; CCimg = 0.	 !! zero image

     !! Find shot location to determine migration aperature
     !! Not important when aperature is full width

     ishot = maxloc(cabs(swf(:,1)));
     iis=ishot(1)
     mxl = max(   1,iis-migapp/2)
     mxr = min(oldn,iis+migapp/2)
     if (mxl .eq. 1   ) mxr=migapp
     if (mxr .eq. oldn) mxl=oldn-migapp+1
     write(0,*) 'ISHOT',iis,mxl,mxr
     write(0,*) 'RWF MIN/MAX', minval(abs(rwf)),maxval(abs(rwf))
     write(0,*) 'SWF MIN/MAX', minval(abs(swf)),maxval(abs(swf))

     !! Window fields and mask by aperture
     cutfields = fields(mxl:mxr,:,:)
     cutmask   = mask  (mxl:mxr,:  )

     !! Call to migration subroutine

     !! INPUT
     !!		1) rwf - receiver wavefield
     !!		2) swf - source wavefield
     !!		4)	ii - the shot number
     !!		5) cutfields - the values of coefficients (for SSF)
     !!		6) refs - reference values for (PSPI/SSF)
     !! 		7) cutmask - integer for closest reference value
     !!
     !!	INPUT/OUTPUT
     !!		3) Pimg - the RWE-based image

     call wemig(rwf,swf,Pimg,ii,cutfields,refs,cutmask )

     !! Write out RWE-based image field
     call rsf_write(Rimage,Pimg)

     !! RWE-to-Cartesian Interpolation Program
     write(0,*) "Pimg", maxval(Pimg),minval(Pimg)
     write(0,*) "Real Rays", maxval(raysin(:,:,1)),minval(raysin(:,:,1))
     write(0,*) "Imag Rays", maxval(raysin(:,:,2)),minval(raysin(:,:,2))
     write(0,*) "Gnorm", maxval(gnorm),minval(gnorm)

     call C2R_run(Pimg,raysin(:,mxl:mxr,:),CCimg,gnorm(:,mxl:mxr))
     write(0,*) "CCimg",maxval(CCimg),minval(CCimg)

     !! Write out Cartesian-based image
     call rsf_write(image,CCimg)

  end do !! . . End Shot-point loop

!  call sep_end_prog()
  call exit()
end program Mrwe2D
