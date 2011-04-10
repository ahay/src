module velinv

use rsf

implicit none

contains

subroutine velinvww( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,rwt,mwt,niter,huber,eps, irls,nstep,srate)
integer nt,it,i, nx,ix, ns,is, huber,niter, iter,irls,nstep
real data(nt,nx),mm(nt,ns)
real dt,t0, dot,rwt,mwt,eps, srate
real x2(nx), s(ns), z2(nt),mask(nx)
real dm(nt,ns), sm(nt,ns)
real rr(nt*nx), dr(nt*nx), sr(nt*nx)
real wrr(nt*nx), trr(nt*nx)
real res(niter+1)

rwt = rwt - 2.
call zero( nt*ns, mm)
call copy( nt*nx, data, rr)
if (irls .eq. 1) then
  rwt = 0.5*rwt
end if
do iter = 0, niter
  if ( irls.eq.0 ) then
    call powwt( iter, rr,nx*nt, rr,wrr, rwt,eps,huber,srate)
  else if ( mod(iter,nstep).eq.0 ) then
    do i=1,nt*nx
      trr(i) = rr(i)
    end do
    call powwt( iter, rr,nx*nt, trr,wrr, rwt,eps,huber,srate)
  else
    call powwt( iter, rr,nx*nt, trr,wrr, rwt,eps,huber,srate)
  end if
  call velxf( 1, 0, nt,dt,t0, x2,nx, wrr,mask, s,ns, dm,z2)
  if ( mwt.ne.0. ) then
    call powwt( iter, dm,ns*nt, mm,dm, mwt,eps,0,srate)
  end if
  call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2)
  if ( irls.eq.1 ) then
    call powwt( iter, dr,nx*nt, trr,dr, rwt,eps,huber,srate)
  end if
  call cgstep(iter, nt*ns, mm,dm,sm, nt*nx, rr,dr,sr)
  res(iter+1) = dot(nt*nx,rr,rr)
end do
! call slice('res.H',4,niter+1,1,res)
return
end subroutine velinvww


subroutine powwt(iter, rr,n, wt, wrr, pow,eps,huber,srate)
integer iter, n, i,huber
real rr(n),wt(n),wrr(n)
real pow,eps,small,minwt,maxwt&
  &,srate
real wwt(n)
do i=1,n
  wwt(i) = abs(wt(i))
end do
if ( pow<=0. ) then
  call minmax(n,wwt,minwt,maxwt)
  small = srate*maxwt
  do i=1,n
    if (wwt(i) < small) then
      wwt(i) = small
    end if
  end do
end if
if ( huber.eq.1 ) then
  call minmax(n,rr,minwt,maxwt)
  small = eps*maxwt
  do i=1,n
    if ( abs(rr(i)) > small ) then
      wrr(i) = rr(i) * wwt(i)**pow
    else
      wrr(i) = rr(i)
    end if
  end do
else
  if (iter > 0 .and. pow .ne. 0) then
    do i=1,n
      wrr(i) = rr(i) * wwt(i)**pow
    end do
  else
    do i=1,n
      wrr(i) = rr(i)
    end do
  end if
end if
return
end subroutine powwt

subroutine cgstep( iter, n, x, g, s, m, rr, gg, ss)
integer i, n, m, iter
real x(n), rr(m), g(n), gg(m), s(n), ss(m)
! solution, residual
! gradient, conjugate gradient
! step, conjugate step
real dot, sds, gdg, gds, determ, gdr, sdr, alfa, beta
if ( iter .eq. 0 ) then
  do i= 1, n
    s(i) = 0.
  end do
  do i= 1, m
    ss(i) = 0.
  end do
  if ( dot(m,gg,gg).eq.0 ) then
    call erexit('cgstep: grad vanishes identically')
  end if
  alfa = dot(m,gg,rr) / dot(m,gg,gg)
  beta = 0.
else
! search plane by solving 2-by-2
  gdg = dot(m,gg,gg) ! G . (R - G*alfa - S*beta) = 0
  sds = dot(m,ss,ss) ! S . (R - G*alfa - S*beta) = 0
  gds = dot(m,gg,ss)
  determ = gdg * sds - gds * gds + 1.e-15
  gdr = dot(m,gg,rr)
  sdr = dot(m,ss,rr)
  alfa = ( sds * gdr - gds * sdr ) / determ
  beta = (-gds * gdr + gdg * sdr ) / determ
end if
do i= 1, n ! s = model step
  s(i) = alfa * g(i) + beta * s(i)
end do
do i= 1, m ! ss = conjugate
  ss(i) = alfa * gg(i) + beta * ss(i)
end do
do i= 1, n ! update solution
  x(i) = x(i) + s(i)
end do
do i= 1, m ! update residual
  rr(i) = rr(i) - ss(i)
end do
return
end subroutine cgstep


real function dot( n, x, y )
integer i, n
real val, x(n), y(n)
val = 0.
do i= 1, n
  val = val + x(i) * y(i)
end do
dot = val
return
end function dot

subroutine minmax( n1, a, min, max)
!
! pick min and max values among the values a(i)
!
integer n1,i1
real a(n1), min, max
max = -1.e10
do i1 = 1, n1
  if ( a(i1) > max ) then
    max = a(i1)
  else
    max = max
  end if
end do
min = 1.e10
do i1 = 1, n1
  if ( a(i1) < min ) then
    min = a(i1)
  else
    min = min
  end if
end do
return
end subroutine minmax


end module velinv
