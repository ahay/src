!  Inverse velocity spectrum with interpolation by modeling from inversion result
!
!
!
!!$  Copyright (C) 2010 Hansung University
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
program Mvelinvww
  !!
  !! Mvelinvww < input.H ns= ds= os= niter= cmpout= > output.H
  !!
  !!	adj = 0  : from velocity-domain(t,s) to cmp-gather(t,x)
  !!	adj = 1  : from cmp-gather(t,x) to velocity-domain(t,s)
  !!
  !! 
!!!!!!!!!!!!!!!!!!
  !! Required modules

  use rsf
  use velinv

  implicit none
  real, allocatable, dimension (:,:) :: cmp
  real, allocatable, dimension (:,:) :: intcmp
  real, allocatable, dimension (:,:) :: vel
  real, allocatable, dimension (:) :: mask
  real, allocatable, dimension (:) :: oh2
  real, allocatable, dimension (:) :: z2
  real, allocatable, dimension (:) :: s
  integer :: nt,nh,ns
  integer :: it,ih,is
  real    :: dt,dh,ds
  real    :: ot,oh,os, rwt, mwt, srate,eps, h, z
  integer :: niter, savevel, huber,irls,nstep
  type (file) :: infile, outfile, vtr

  call sf_init()

  infile = rsf_input()

  call from_par(infile,"n1",nt);  call from_par(infile,"n2",nh)
  call from_par(infile,"d1",dt);  call from_par(infile,"d2",dh)
  call from_par(infile,"o1",ot);  call from_par(infile,"o2",oh)
 
  call from_par("ns",ns,nh);  call from_par("ds",ds,0.001);  call from_par("os",os,0.00000001)

  call from_par("huber",huber,0)
  call from_par("irls",irls,0)
  call from_par("nstep",nstep,1)
  call from_par("niter",niter,20)
  call from_par("rwt",rwt,0.)
  call from_par("mwt",mwt,0.)
  call from_par("srate",srate,0.01)
  call from_par("epw",eps,0.01)
  call from_par("savevel",savevel,0)

  allocate (cmp(nt,nh))
  allocate (intcmp(nt,nh))
  allocate (vel(nt,ns))
  allocate (mask(nh))
  allocate (oh2(nh))
  allocate (z2(nt))
  allocate (s(ns))

  outfile = rsf_output()
  call to_par(outfile,"n1",nt);  call to_par(outfile,"n2",nh)
  call to_par(outfile,"d1",dt);  call to_par(outfile,"d2",dh)
  call to_par(outfile,"o1",ot);  call to_par(outfile,"o2",oh)

  do ih= 1,nh
     h = oh + dh*(ih-1);        oh2(ih) = h*h
  end do

  do it= 1,nt
     z = ot + dt*(it-1);        z2(it) = z*z
  end do

  do is= 1,ns
     s(is) = os + ds*(is-1);
  end do

  call rsf_read(infile, cmp)

  do ih=1,nh
     mask(ih) = dot(nt,cmp(1,ih),cmp(1,ih))
  end do
  call velinvww( nt,dt,ot, oh2,nh,cmp, s,ns,vel, z2,mask,rwt,mwt,niter&
          &,huber,eps,irls,nstep,srate)

  if ( savevel.eq.1 ) then
     vtr = rsf_output("velout")
     call to_par(vtr,"n1",nt);     call to_par(vtr,"n2",ns)
     call to_par(vtr,"d1",dt);     call to_par(vtr,"d2",ds)
     call to_par(vtr,"o1",ot);     call to_par(vtr,"o2",os)
     call rsf_write(vtr, vel)
  end if

  do ih= 1,nh
     mask(ih) = 1.
  end do

  call velxf( 0,0, nt,dt,ot, oh2,nh,intcmp, mask, s,ns, vel, z2)

  do ih= 1,nh
     mask(ih) = dot(nt,cmp(1,ih),cmp(1,ih))
  end do

  do ih= 1,nh
     if ( mask(ih) .eq. 0. ) then
        do it=1,nt
           cmp(it,ih) = intcmp(it,ih)
        end do
     end if
  end do

  call rsf_write( outfile, intcmp )

end program Mvelinvww
