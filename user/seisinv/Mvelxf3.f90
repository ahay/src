!  3-D Velocity transform for generating velocity spectra and its corresponding hyperbolic modeling
!
!
!
!!$  Copyright (C) 2012 China University of Petroleum (East China)
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
  !!	adj = 0  : from velocity-domain(t,v) to cmp-gather(t,x,y)
  !!	adj = 1  : from cmp-gather(t,x,y) to velocity-domain(t,v)
  !!
  !! 
!!!!!!!!!!!!!!!!!!
  !! Required modules
  use rsf
  use hradon

  implicit none

  real, allocatable, dimension (:,:) :: data
  real, allocatable, dimension (:,:) :: modl
  real, allocatable, dimension (:) :: mask
  real, allocatable, dimension (:) :: h2
  real, allocatable, dimension (:) :: z2
  real, allocatable, dimension (:) :: s2

  integer :: it,ix,iy,ih,iv,adj
  integer :: nt,nx,ny,nh,nv
  integer :: ic, nc                               ! number of CMP gathers
  real    :: dt,dx,dy,dv,x,y,z,s
  real    :: ot,ox,oy,ov
  type (file) :: cmp, vtr

  integer(kind=8) sf_leftsize
  external sf_leftsize 

  call sf_init()

  call from_par("adj",adj,0) 

  !	adj = 0  : from velocity-domain(t,s) to cmp-gather(t,x)
  !	adj = 1  : from cmp-gather(t,x) to velocity-domain(t,s)

  if (adj.eq.0) then
     vtr = rsf_input()
     call from_par(vtr,"n1",nt)
     call from_par(vtr,"n2",nv)
     call from_par(vtr,"d1",dt)
     call from_par(vtr,"d2",dv)
     call from_par(vtr,"o1",ot)
     call from_par(vtr,"o2",ov)
     call from_par("nx",nx,nv)
     call from_par("ny",ny,1) 
     call from_par("dx",dx,0.01)
     call from_par("dy",dy,0.01)
     call from_par("ox",ox,0.)
     call from_par("oy",oy,0.) 
  else
     cmp = rsf_input()   
     call from_par(cmp,"n1",nt)
     call from_par(cmp,"n2",nx)
     call from_par(cmp,"n3",ny)
     call from_par(cmp,"d1",dt)
     call from_par(cmp,"d2",dx)
     call from_par(cmp,"d3",dy)
     call from_par(cmp,"o1",ot)
     call from_par(cmp,"o2",ox)
     call from_par(cmp,"o3",oy)
     call from_par("nv",nv,nx) 
     call from_par("dv",dv,0.01) 
     call from_par("ov",ov,1.5) 
  end if

  nh=nx*ny
  allocate (data(nt,nh))
  allocate (modl(nt,nv))
  allocate (mask(nh))
  allocate (h2(nh))
  allocate (z2(nt))
  allocate (s2(nv))

  do ix =1,nh
     mask(ix) = 1.
  end do

  do iy = 1, ny
     y = oy + dy*(iy-1)
     do ix =1,nx
        x = ox + dx*(ix-1)
        ih = ix+(iy-1)*nx
        h2(ih) = x*x+y*y
     end do
  end do

  do it =1,nt
     z = ot + dt*(it-1);     z2(it) = z*z
  end do

  do iv =1,nv
     s = 1./(ov + dv*(iv-1)); s2(iv) = s*s
  end do

  if (adj.eq.0) then
     ! set output
     cmp = rsf_output()

     ! get the number of CMP gathers
     nc = sf_leftsize(vtr,2)

     call to_par(cmp,"n1",nt)
     call to_par(cmp,"n2",nx)
     call to_par(cmp,"n3",ny)
     call to_par(cmp,"n4",nc)
     call to_par(cmp,"d1",dt)
     call to_par(cmp,"d2",dx)
     call to_par(cmp,"d3",dy)
     call to_par(cmp,"o1",ot)
     call to_par(cmp,"o2",ox)
     call to_par(cmp,"o3",oy)

     write(0,*) 'Hyperbolic modeling begin...'

     do ic = 1, nc
        write(0,*) 'ic/nc=',ic,'/',nc
        call rsf_read(vtr, modl)
        call velxf( 0,0, nt,dt,ot, h2,nh, data,mask, s2,nv, modl,z2)
        call rsf_write(cmp, data)
     end do
     write(0,*) 'Hyperbolic modeling end!!!'

  else
     vtr = rsf_output()

     ! get the number of CMP gathers
     nc = sf_leftsize(cmp,3)

     call to_par(vtr,"n1",nt)
     call to_par(vtr,"n2",nv)
     call to_par(vtr,"n3",1)
     call to_par(vtr,"n4",nc)
     call to_par(vtr,"d1",dt)
     call to_par(vtr,"d2",dv)
     call to_par(vtr,"o1",ot)
     call to_par(vtr,"o2",ov)

     write(0,*) 'Hyperbolic tranform begin...'
     
     do ic = 1, nc
        write(0,*) 'ic/nc=',ic,'/',nc
        call rsf_read(cmp, data)
        call velxf( 1,0, nt,dt,ot, h2,nh, data,mask, s2,nv, modl,z2)
        call rsf_write(vtr, modl)
     end do
     write(0,*) 'Hyperbolic transform end!!!'

  end if

end program Mvelxf 
