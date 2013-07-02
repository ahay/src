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
module hradon
use rsf
implicit none
contains

subroutine velxf( adj,add, nt,dt,ot, x2,nx, data,mask, s,ns, modl,z2)
integer           adj,add, nt,it,    ix,nx,           is,ns, iz, endt
real              data(nt,nx), modl(nt,ns)
real              t,dt,ot,    x2(nx), s(ns), mask(nx), z2(nt), x2s
real              ft,gt
call adjnull( adj, add, modl, nt*ns, data, nt*nx )
if ( adj .eq. 0) then
  do is=  1, ns  
    do ix=  1, nx  
      x2s = x2(ix)*s(is)
      if ( x2s > z2(nt) ) then
        exit
      end if
      if ( mask(ix).eq.0. ) then
        cycle
      end if
      endt = int( (sqrt( z2(nt) - x2s) - ot)/dt + 1.)
      do iz=  1, endt  
        t = sqrt( z2(iz) + x2s )
        it = 1.+(t-ot)/dt
        ft = (it*dt-t)/dt 
        gt = 1. - ft
        data(it,  ix) = data(it,  ix) + ft*modl(iz,is)
        data(it+1,ix) = data(it+1,ix) + gt*modl(iz,is)
      end do
    end do
  end do
else
  do is=  1, ns 
    do ix=  1, nx  
      x2s = x2(ix)*s(is)
      if ( x2s > z2(nt) ) then
        exit
      end if
      if ( mask(ix).eq.0. ) then
        cycle
      end if
      endt = int( (sqrt( z2(nt) - x2s) - ot)/dt + 1.)
      do iz=  1, endt  
        t = sqrt( z2(iz) + x2s )
        it = 1.+(t-ot)/dt
        ft = (it*dt-t)/dt 
        gt = 1. - ft
        !write(0,*) iz,z2(iz),endt,t,ot,dt,it
        modl(iz,is) = modl(iz,is) + ft*data(it,ix) + gt*data(it+1,ix)
      end do
    end do
  end do
end if 
return
end subroutine velxf

subroutine adjnull( adj, add, x, nx, y, ny )
integer ix, iy, adj, add, nx, ny
real x( nx), y( ny )
if ( add .eq. 0 ) then
  if ( adj .eq. 0 ) then
    do iy= 1, ny
      y(iy) = 0.
    end do
  else
    do ix= 1, nx
      x(ix) = 0.
    end do
  end if
end if
return
end subroutine adjnull

end module hradon
