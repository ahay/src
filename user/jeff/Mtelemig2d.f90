
! 2D Teleseismic migration code allowing: forward- and backscattering imaging conditions and different source- and receiver-side velocity models

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
  use rsf              !! Basic RSF library
  use telemig2d_wemig  !! Migration driver module  
  implicit none

  integer :: nx,nz,nh,testint,nw,ns,is,ix,iz,iw,ntaper
  real  :: ox,dx,oz,dz,oh,dh,ow,dw,ds,os,eps
  logical :: source_norm,verbose,forward
  real, allocatable,dimension(:,:,:) :: Svel,Rvel
  real, allocatable,dimension(:,:) :: img
  real, allocatable, dimension(:,:,:) :: cig
  complex,allocatable, dimension(:,:,:) :: rwf,swf
 
   type(file) :: infile,outfile,Svelfile,rwffile,swffile,cigfile
  
   call sf_init() 
   infile = rsf_input("in") ! Input velocity file RWF velR(nx,nz,ns)
   Svelfile = rsf_input("Svel") ! Input velocity file SWF velS(nx,nz,ns)
   rwffile = rsf_input("rwf") ! Input receiver wavefield (nx,nw,ns)
   swffile = rsf_input("swf") ! Input source wavefield (nx,nw,ns)
   outfile = rsf_output("out") ! Output image (nx,nh,nz)
   cigfile = rsf_output("cig") ! Output CIG file
   
  !! . . Geometry information
  call from_par(Svelfile,"n1",nx) !! # of Horizontal Samples
  call from_par(Svelfile,"o1",ox) !! Horizontal array offset
  call from_par(Svelfile,"d1",dx) !! Horizontal sampling interval 
  call from_par(Svelfile,"n2",nz) !! # of Depth Samples
  call from_par(Svelfile,"o2",oz) !! Depth array offset
  call from_par(Svelfile,"d2",dz) !! Depth sampling interval
	
  !! . . Wavefield information
   call from_par(rwffile,"n2",nw) !! # of frequencies
   call from_par(rwffile,"o2",ow) !! Starting frequencies
   call from_par(rwffile,"d2",dw) !! Frequency sampling interval
   call from_par(rwffile,"n3",ns) !! # of sources
   call from_par(rwffile,"o3",os) !! # of sources
   call from_par(rwffile,"d3",ds) !! # of sources

   !! . .Test section .. catch any errors!
   call from_par(Svelfile,"n1",testint); if (testint .ne. nx) call sf_error("Wrong # of X samples in Svel file")
   call from_par(infile,"n1",testint); if (testint .ne. nx) call sf_error("Wrong number of RWF samples")
   call from_par(swffile,"n1",testint); if (testint .ne. nx) call sf_error("Wrong number of SWF samples")
   call from_par(swffile,"n3",testint) !! # of sources in Velocity
   if (testint .ne. ns) call sf_error("Wrong number of sources between RWF and SWF"); 
   call from_par(Svelfile,"n3",testint) !! # of sources in Velocity
   if (testint .ne. ns) call sf_error("Wrong number of sources between RWF and Velocity model")
   
   !! . . Obtain information on offsets
   call from_par("nh",nh,0)  ! Number of subsurface offsets (between 1 and 128)
   call from_par("source_norm",source_norm,.false.) ! Normalize the image by the power of the SWF (T/F)
   call from_par("verbose", verbose, .false.) ! Verbose (T/F)
   call from_par("forward", forward, .false.) ! Forward scattering (T/F)
   call from_par("ntaper", ntaper, 40) ! Taper on the side boundaries (npts)
   call from_par("eps", eps, 0.01) ! Taper on the side boundaries (npts)
   
   if (nh .lt. 0 .or. nh .gt. 128) call sf_error("Number of offsets not between 1<nh<128")
   nh = 2*nh+1; dh = dx;  oh=-(nh-1)*dh
   
   if (nh .eq. 0)  then
      write(0,*) 'Not computing CIGs'
   else
      write(0,*) 'Using CIGs of size: ',nh
   endif
   
   !! . . Output geometry
   call to_par(outfile,"n1",nx)
   call to_par(outfile,"o1",ox)
   call to_par(outfile,"d1",dx)
   call to_par(outfile,"n2",nz)
   call to_par(outfile,"o2",oz)
   call to_par(outfile,"d2",dz)
   call to_par(outfile,"esize",4)
   
   !! . . Output geometry
   call to_par(cigfile,"n1",nh)
   call to_par(cigfile,"o1",oh)
   call to_par(cigfile,"d1",dh)
   call to_par(cigfile,"n2",nx)
   call to_par(cigfile,"o2",ox)
   call to_par(cigfile,"d2",dx)
   call to_par(cigfile,"n3",nz)
   call to_par(cigfile,"o3",oz)
   call to_par(cigfile,"d3",dz)
   call to_par(cigfile,"esize",4)
   
   !! . . Initialize propagation module
   call wemig_axes_init(nx,ox,dz,nz,oz,dz,nw,ow,dw,nh,oh,dh,source_norm,verbose,forward,ntaper,eps)
   
   write(0,*) 'nx,ox,dx ... ',nx,ox,dx
   write(0,*) 'nz,oz,dz ... ',nz,oz,dz
   write(0,*) 'nw,ow,dw ... ',nw,ow,dw
   write(0,*) 'nh,oh,dh ... ',nh,oh,dh
   write(0,*) 'ns,os,ds ... ',ns,os,ds
   
   !! . .  MIGRATION SECTION
   !! . . Allocate space for receiver wavefield and local velocities
   allocate ( rwf (nx,nw,ns) , swf(nx,nw,ns), Svel(nx,nz,ns), Rvel(nx,nz,ns), img(nx,nz) )
   
   allocate(cig(nh,nx,nz))
   
   !! . . Read in velocity files
   call rsf_read(Svelfile,Svel)
   call rsf_read(infile,Rvel)
   
   !! . . Read in RWF
   call rsf_read(rwffile,rwf)
   
   !! . . Read in SWF
   call rsf_read(swffile,swf)
   
   !! . . Loop over shot points: inline and streamers
   SHOTS: do is = 1,ns
      
      !! . . Call single-shot migration program
      call wemig(Svel(:,:,is) , Rvel(:,:,is) , swf(:,:,is) , rwf(:,:,is) , img, cig)   
      
   end do SHOTS
   
   !! . . write out image file
   call rsf_write(outfile,img)
   
   call rsf_write(cigfile,cig)
   
   !! . . Exit file
   call exit()
   
 end program Mtelemig2d
