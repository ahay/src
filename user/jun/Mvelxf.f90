!  Velocity transform for generating velocity spectra and its corresponding hyperbolic modeling
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
program Mvelxf
  !!
  !!	Velxf < input.H adj=[ 1 or 0 ] nx= dx= or ns= ds= > output.H
  !!
  !!	adj = 0  : from velocity-domain(t,s) to cmp-gather(t,x)
  !!	adj = 1  : from cmp-gather(t,x) to velocity-domain(t,s)
  !!
  !! 
!!!!!!!!!!!!!!!!!!
  !! Required modules
  use rsf
  use veltran

  implicit none

  real, allocatable, dimension (:,:) :: data
  real, allocatable, dimension (:,:) :: modl
  real, allocatable, dimension (:) :: mask
  real, allocatable, dimension (:) :: x2
  real, allocatable, dimension (:) :: z2
  real, allocatable, dimension (:) :: s

  integer :: it,ix,is, adj
  integer :: nt,nx,ns
  real    :: dt,dx,ds, x, z
  real	  :: ot,ox,os
  type (file) :: cmp, vtr

  call sf_init()

  call from_par("adj",adj,0) 

  !	adj = 0  : from velocity-domain(t,s) to cmp-gather(t,x)
  !	adj = 1  : from cmp-gather(t,x) to velocity-domain(t,s)

  if (adj.eq.0) then
     vtr = rsf_input()
     call from_par(vtr,"n1",nt)
     call from_par(vtr,"n2",ns)
     call from_par(vtr,"d1",dt)
     call from_par(vtr,"d2",ds)
     call from_par(vtr,"o1",ot)
     call from_par(vtr,"o2",os)
     call from_par("nx",nx,ns) 
     call from_par("dx",dx,10.) 
     call from_par("ox",ox,0.) 
  else
     cmp = rsf_input()   
     call from_par(cmp,"n1",nt)
     call from_par(cmp,"n2",nx)
     call from_par(cmp,"d1",dt)
     call from_par(cmp,"d2",dx)
     call from_par(cmp,"o1",ot)
     call from_par(cmp,"o2",ox)
     call from_par("ns",ns,nx) 
     call from_par("ds",ds,0.001) 
     call from_par("os",os,0.00000001) 
  end if

  allocate (data(nt,nx))
  allocate (modl(nt,ns))
  allocate (mask(nx))
  allocate (x2(nx))
  allocate (z2(nt))
  allocate (s(ns))

  if (adj.eq.0) then
     call rsf_read(vtr, modl)
  else
     call rsf_read(cmp, data)
  end if

  do ix =1,nx  
     mask(ix) = 1.
  end do
  do ix =1,nx  
     x = ox + dx*(ix-1);     x2(ix) = x*x
  end do
  do it =1,nt  
     z = ot + dt*(it-1);     z2(it) = z*z
  end do
  do is =1,ns  
     s(is) = os + ds*(is-1)
  end do

  call velxf( adj,0, nt,dt,ot, x2,nx, data,mask, s,ns, modl,z2)

  if (adj.eq.0) then
     cmp = rsf_output()   
     call to_par(cmp,"n1",nt)
     call to_par(cmp,"n2",nx)
     call to_par(cmp,"d1",dt)
     call to_par(cmp,"d2",dx)
     call to_par(cmp,"o1",ot)
     call to_par(cmp,"o2",ox)
     call rsf_write(cmp, data)
  else
     vtr = rsf_output()
     call to_par(vtr,"n1",nt)
     call to_par(vtr,"n2",ns)
     call to_par(vtr,"d1",dt)
     call to_par(vtr,"d2",ds)
     call to_par(vtr,"o1",ot)
     call to_par(vtr,"o2",os)
     call rsf_write(vtr, modl)
  end if


end program Mvelxf 
