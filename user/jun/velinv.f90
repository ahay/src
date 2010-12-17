module velinv
use rsf
use veltran
implicit none
contains

subroutine velinvww( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,rwt,mwt,niter,huber,eps, irls,nstep,srate)
  integer nt,it,i, nx,ix, ns,is, huber,niter, iiter,iter,irls,nstep
  integer ndata, nmodl
  real data(nt*nx), mm(nt,ns)
  real dt,t0, dot,rwt,mwt,eps, srate
  real x2(nx), s(ns), z2(nt), mask(nx)
  real dm(nt*ns), sm(nt*ns)
  real wdm(nt*ns), tmm(nt*ns)
  real rr(nt*nx), dr(nt*nx), sr(nt*nx)
  real wrr(nt*nx), trr(nt*nx)
  real res(niter+1)
  
  rwt = 0.5 * (rwt - 2.)   ! residual weight exponent (p-2)/2 for L-p norm

  ndata = nt*nx
  nmodl = nt*ns

  call zero( nmodl, mm)
  call copy( ndata, data,  rr)
  call copy( ndata,   rr, trr)

  if ( irls.eq.0 ) then  ! CGG inversion

     do iter = 0, niter
        call powwt( iter, rr,ndata, rr,wrr, rwt,eps,huber,srate)
        call velxf( 1, 0, nt,dt,t0, x2,nx, wrr,mask, s,ns, dm,z2)
        if ( mwt.ne.0. ) then
           call powwt( iter, dm,nmodl, mm,dm, mwt,eps,0,srate)
        end if
        call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2)
        call cgstep( iter, nmodl, mm,dm,sm, ndata, rr,dr,sr)
        !        call ddot(ndata,rr,rr,res(iter+1))
        !        write(0,*) 'iter ',iter, res(iter+1)
     end do
     
  else                  ! IRLS inversion

     if ( mwt.ne.0 ) then
        mwt = 0.5 * (2. - mwt)   ! model weight exponent (2-p)/2 for L-p norm
     end if

    do iter = 0, niter
       call velxf( 1, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, dm,z2)
       call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2)
       call cgstep( iter, nmodl, mm,dm,sm, ndata, rr,dr,sr)
    end do

    do iiter = 1, nstep

           call copy( nmodl, mm, tmm)
           call copy( ndata, rr, trr)
           call powwt( iiter, mm,nmodl, tmm,wdm, mwt,eps,0,srate)
           call velxf( 0, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, wdm,z2)
           do i = 1, ndata
              rr(i) = data(i) - rr(i)
           end do
           call powwt( iiter, rr,ndata, trr,rr, rwt,eps,huber,srate)

           do iter = 1, niter
              call powwt( iter, rr,ndata, trr,wrr, rwt,eps,huber,srate)        
              call velxf( 1, 0, nt,dt,t0, x2,nx, wrr,mask, s,ns, dm,z2)
              call powwt( iter, dm,nmodl, tmm, dm, mwt,eps,0,srate)
        
              call powwt( iter, dm,nmodl, tmm,wdm, mwt,eps,0,srate)
              call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, wdm,z2)
              call powwt( iter, dr,ndata, trr, dr, rwt,eps,huber,srate)
        
              call cgstep(iter, nmodl, mm,dm,sm, ndata, rr,dr,sr)
              !       call ddot(ndata,rr,rr,res(iter+1))
              !       write(0,*) 'iter ',iter, res(iter+1)
           end do

           call powwt( iter, mm,nmodl, tmm,mm, mwt,eps,0,srate)
     end do
  end if
  return

end subroutine velinvww


subroutine powwt(iter, rr,n, wt, wrr, pow,eps,huber,srate)
  integer iter, n, i,huber
  real rr(n)     !     input vector
  real wt(n)     ! weighting vector
  real wrr(n)    !  weighted vector
  real pow       ! weighting power
  real eps,small,minval,maxval,srate
  real arr(n), awt(n)

  if ( pow.eq.0. .or. iter.eq.0 ) then
     call copy(n,rr,wrr)
     return
  end if

  do i= 1, n
     awt(i) = abs(wt(i))
  end do

  if ( pow < 0. ) then
     call minmax( n,awt,minval,maxval )
     small = srate * maxval
     do i= 1, n
        if (awt(i) < small) then
           awt(i) = small
        end if
     end do
  end if

  if ( huber .eq. 1 ) then
     do i= 1, n
        arr(i) = abs(rr(i))
     end do
     call minmax( n,arr,minval,maxval )
     small = eps * maxval
     do i= 1, n
        if ( arr(i) > small ) then
           wrr(i) = rr(i) * awt(i)**pow
        else
           wrr(i) = rr(i)
        end if
     end do
  else
     do i= 1, n
        wrr(i) = rr(i) * awt(i)**pow
     end do
  end if

  return
end subroutine powwt


subroutine cgstep( iter, n, x, g, s, m, rr, gg, ss)
  integer i, n, m, iter
  real x(n), rr(m), g(n), gg(m), s(n), ss(m)
! solution, residual
! gradient, conjugate gradient
! step, conjugate step
  real dot, den,num, sds, gdg, gds, determ, gdr, sdr, alfa, beta
  if ( iter .eq. 0 ) then
     do i= 1, n
        s(i) = 0.
     end do
     do i= 1, m
        ss(i) = 0.
     end do
     call ddot(m,gg,gg,dot)
     if (dot.eq.0) then
        return
        !    call erexit('cgstep: grad vanishes identically')
     end if
     call ddot(m,gg,rr,den)
     call ddot(m,gg,gg,num)
     alfa = den/num
     beta = 0.
  else
     ! search plane by solving 2-by-2
     call ddot(m,gg,gg,gdg) ! G . (R - G*alfa - S*beta) = 0
     call ddot(m,ss,ss,sds) ! S . (R - G*alfa - S*beta) = 0
     call ddot(m,gg,ss,gds)
     determ = gdg * sds - gds * gds + 1.e-15
     call ddot(m,gg,rr,gdr)
     call ddot(m,ss,rr,sdr)
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

subroutine minmax( n, a, min, max)
!
! pick min and max values among the values a(i)
!
  integer :: n, i
  real    :: a(n), min, max
  max = 0.0
  do i = 1, n
     if ( a(i) > max ) then
        max = a(i)
     end if
  end do
  min = 1.e20
  do i = 1, n
     if ( a(i) < min ) then
        min = a(i)
     end if
  end do
  return
end subroutine minmax

subroutine copy(n,x,y)
integer :: n, i 
real  :: x(n), y(n) 
do i=1,n
y(i) = x(i)
end do
end subroutine copy

subroutine zero(n,x)
integer :: n, i 
real  :: x(n) 
do i=1,n
x(i) = 0.
end do
end subroutine zero 

subroutine ddot(n,x,y,dot)
integer :: i,n
real    :: dot,x(n),y(n)
dot = 0.
do i=1,n
dot = dot + x(i)*y(i)
end do
end subroutine ddot

real function dot( n, x, y ) result (val)
integer :: i, n
real    :: x(n), y(n)
val = 0.
do i= 1, n
  val = val + x(i) * y(i)
end do
end function dot


end module velinv
