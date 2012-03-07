! Teleseismic 2D Migration Code for forward- and backscattering modes

!!$  Copyright (C) 2012 The University of Western Australia
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
program Mtelemig2d
  use rsf              !! Basic SEPlib I/O 
  use fft
  use telemig2D_wemig  !! Migration driver module

  implicit none

  integer :: nx,nz,nh,testint,nw,ns,is
  real :: ox,oz,dx,dz,ix,iz,dh,oh,testreal,dw,ow

  real, allocatable,dimension(:,:) :: Svel,Rvel
  real, allocatable,dimension(:,:,:):: img
  complex,allocatable, dimension(:,:) :: rwf,swf

  character (len=128) :: name
 
   type(file) :: infile,outfile,Svelfile,Rvelfile,swffile
  
   call sf_init() 
   Svelfile = rsf_input("Svel") !! Input velocity file SWF velS(nx,nz,ns)
   Rvelfile = rsf_input("Rvel") !! Input velocity file RWF velR(nx,nz,ns)
   infile = rsf_input("in") !! Input receiver wavefield (nx,nw,ns)
   swffile = rsf_input("swf") !! Input source wavefield (nx,nw,ns)
   outfile=rsf_output("out") !! Output image (nx,nh,nz)
   
  !! . . Geometry information
  call from_par(Svelfile,"n1",nx) !! # of Horizontal Samples
  call from_par(Svelfile,"n2",nz) !! # of Depth Samples
  call from_par(Svelfile,"o1",ox) !! Horizontal array offset
  call from_par(Svelfile,"o2",oz) !! Depth array offset
  call from_par(Svelfile,"d1",dx) !! Horizontal sampling interval 
  call from_par(Svelfile,"d2",dz) !! Depth sampling interval
	
  !! . . Wavefield information
   call from_par(infile,"n2",nw) !! # of frequencies
   call from_par(infile,"o2",ow) !! Starting frequencies
   call from_par(infile,"d2",dw) !! Frequency sampling interval
   call from_par(infile,"n3",ns) !! # of sources

  	!! . .Test section .. catch any errors!
  	call from_par(Svelfile,"n1",testint); if (testing .ne. nx) call sf_error("Wrong # of X samples in Svel file")
  	call from_par(Rvelfile,"n1",testint); if (testing .ne. nx) call sf_error("Wrong # of X samples in Rvel file")
  
   	call from_par(swffile,"n1",testint); if (testint .ne. nx) call sf_error("Wrong number of RWF samples")
   	call from_par(swffile,"n1",testint); if (testint .ne. nx) call sf_error("Wrong number of SWF samples")
   	if(ow .eq. 0.) call sf_error("Do not use zero frequency")
    	    call from_par(swffile,"n3",testing) !! # of sources in Velocity
   	if (testint .ne. ns) call sf_error("Wrong number of sources between RWF and SWF"); 
   	     call from_par(velSfile,"n3",testing) !! # of sources in Velocity
   	if (testint .ne. ns) call sf_error("Wrong number of sources between RWF and Velocity model")

  	!! . . Obtain information on offsets
  	call from_param("nh",nh,1) 
  	if (nh .lt. 1 .or. nh .gt. 128) call sf_error("Number of offsets not between 1<nh<128")
  	nh = 2*nh+1; dh = dx;  oh=-nh*dh

  	!! . . Initialize propagation module
  	call wemig_axes_init(nx,ox,dz,nz,oz,dz,nw,ow,dw,nh,oh,dh)
  
  	!! . .  MIGRATION SECTION
    !! . . Allocate space for receiver wavefield and local velocities
  	allocate ( rwf (nx,nw) , swf(nx,nw), Svel(nx,nz), Rvel(nx,nz), img(nx,nh,nw) )

  	!! . . Loop over shot points: inline and streamers
  	SHOTS: do is = 1,ns

		!! . . Zero files each time
		vel = 0.; rwf=0.; swf =0.
	
		!! . . Read in velocity files
		call rsf_read(Svelfile,velS,4*nx*nz)
		call rsf_read(infile,velS,4*nx*nz)
	
		!! . . Read in RWF
		call rsf_read(rwffile,rwf,4*nx*nw)

		!! . . Read in SWF
		call rsf_read(swffile,swf,4*nx*nw)
	
		!! . . Call single-shot migration program
    	call wemig(velS, velR, swf, rwf, img)   

  	end do INLINE_SHOT

  	!! . . write out image file
  	call rsf_write(outfile,img)

  	!! . . Exit file
  	call exit()

end program Mtelemig2d