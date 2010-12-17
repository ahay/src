module veltran
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
      endt = int( sqrt( z2(nt) - x2s )/dt + 1.)
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
      endt = int( sqrt( z2(nt) - x2s )/dt + 1.)
      do iz=  1, endt  
        t = sqrt( z2(iz) + x2s )
        it = 1.+(t-ot)/dt
        ft = (it*dt-t)/dt 
        gt = 1. - ft
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

end module veltran
