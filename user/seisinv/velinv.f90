module velinv
use rsf
use hradon
implicit none
contains

subroutine velinvww( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,rwt,mwt,niter,huber,eps, irls,nstep,srate, res, mm0, mres)
  integer nt,it,i, nx,ix, ns,is, huber,niter, iiter,iter,irls,nstep
  integer ndata, nmodl
  real data(nt*nx), mm(nt*ns), mm0(nt*ns)
  real dt,t0, dot,rwt,mwt,eps, srate, res0, mres0
  real x2(nx), s(ns), z2(nt), mask(nx)
  real dm(nt*ns), sm(nt*ns)
  real wdm(nt*ns), tmm(nt*ns)
  real rr(nt*nx), dr(nt*nx), sr(nt*nx)
  real wrr(nt*nx), trr(nt*nx)
  real res(niter), mres(niter)
  
  rwt = 0.5 * (rwt - 2.)   ! residual weight exponent (p-2)/2 for L-p norm

  ndata = nt*nx
  nmodl = nt*ns

  call zero( nmodl, mm)
  call copy( ndata, data,  rr)
  call copy( ndata,   rr, trr)

  res = 0.
  res0 = sum(rr*rr)
  mres0 = sqrt(sum(mm0*mm0))

  if ( irls.eq.0 ) then  ! CGG inversion

     do iter = 1, niter

        res(iter) = sum(rr*rr)/res0

        mres(iter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

        write(0,*)'iter=',iter,'mres=',mres(iter),'rr=',res(iter)

        call powwt( iter, rr,ndata, rr,wrr, rwt,eps,huber,srate)
        !write(0,*) nt,dt,t0,nx
        call velxf( 1, 0, nt,dt,t0, x2,nx, wrr,mask, s,ns, dm,z2)
        if ( mwt.ne.0. ) then
           call powwt( iter, dm,nmodl, mm,dm, mwt,eps,0,srate)
        end if
        call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2)
        call cgstep( iter, nmodl, mm,dm,sm, ndata, rr,dr,sr)

     end do
     
  else                  ! IRLS inversion

     niter = niter/(nstep+1)

     if ( mwt.ne.0 ) then
        mwt = 0.5 * (2. - mwt)   ! model weight exponent (2-p)/2 for L-p norm
     end if

    do iter = 1, niter
       
       call velxf( 1, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, dm,z2)
       call velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2)
       call cgstep( iter, nmodl, mm,dm,sm, ndata, rr,dr,sr)
       write(0,*) 'iter/niter',iter,niter,sum(rr*rr)
    end do

    do iiter = 1, nstep

           call copy( nmodl, mm, tmm)
           call copy( ndata, rr, trr)
           call powwt( iiter, mm,nmodl, tmm,wdm, mwt,eps,0,srate)
           call velxf( 0, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, wdm,z2)
           do i = 1, ndata
              rr(i) = data(i) - rr(i)
           end do

           res(iiter) = sum(rr*rr)/res0

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
              write(0,*) 'istep/nstep=',iiter,nstep,'iter/niter',iter,niter,sum(rr*rr)

           end do

           call powwt( iter, mm,nmodl, tmm,mm, mwt,eps,0,srate)

           mres(iiter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

           write(0,*)'iter=',iiter,'mres=',mres(iiter),'rr=',res(iiter)

     end do
  end if
  return

end subroutine velinvww

subroutine velinvsh( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,niter,lamda,delta,res, mm0, mres)
  integer nt,it,i, nx,ix, ns,is,niter,iter
  integer ndata, nmodl
  real data(nt*nx), mm(nt*ns),rm(nt*ns), mm0(nt*ns)
  real dt,t0, dot,eps,delta,lamda, res0, mres0
  real x2(nx), s(ns), z2(nt), mask(nx), res(niter), mres(niter)
  real rr(nt*nx), trr(nt*nx), dm(nt*ns)
!  type (file) :: ftest

  ndata = nt*nx
  nmodl = nt*ns

  call zero( nmodl, mm)
  call zero( ndata, rr)
  call zero( ndata, trr)

  res = 0.

  ndata = nt*nx
  nmodl = nt*ns

  mres0 = sqrt(sum(mm0*mm0))

!  lamda=1000.
!  delta=0.0001
  eps = lamda*delta

!  ftest = rsf_output("test")

!  call to_par(ftest,"n1",nt)
!  call to_par(ftest,"n2",ns)
!  call to_par(ftest,"n3",niter)
  
  do iter=1, niter
!     write(0,*) 'iter/niter',iter,'/',niter,'lamda=',lamda

     if(mod(iter,10).eq.0) lamda=lamda/1.5
     ! Ax
     call velxf( 0, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, mm,z2)

     ! rr=Ax-y
     rr = rr-data

     res(iter)=sum(rr*rr)
     if (iter.eq.1) res0 = res(iter)
     res(iter)=res(iter)/res0

     ! dm=At*rr
     call velxf(1,0,nt,dt,t0,x2,nx,rr,mask,s,ns,dm,z2)

     ! 
     do i = 1, nmodl
         rm(i) = mm(i)-delta*dm(i)
     end do

     call shrink(mm,rm,nmodl,eps)

     mres(iter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

     write(0,*)'iter=',iter,'lamda=',lamda,'mres=',mres(iter),'rr=',res(iter)

!     call rsf_write(ftest,mm)

  end do

end subroutine velinvsh

subroutine velinvlbreg( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,niter,step,alpha, res, mm0, mres)
! solve:
!        min ||m||_1+1/(2*alpha)*||x||_2, subject to Ax=b
! at each iteration...
! trr(i+1)=trr(i)-step*[Ax(i)-b]
! x(i+1)=alpha*shrink(A'trr,1.)

  integer nt,it,i, nx,ix, ns,is,niter,iter
  integer ndata, nmodl
  real data(nt*nx), mm(nt*ns),rm(nt*ns),mm0(nt*ns)
  real dt,t0, dot,step,alpha, res0, mres0
  real x2(nx), s(ns), z2(nt), mask(nx)
  real rr(nt*nx), trr(nt*nx), dm(nt*ns), res(niter), mres(niter)

  ndata = nt*nx
  nmodl = nt*ns

  call zero( nmodl, mm)
  call zero( ndata, rr)
  call zero( ndata, trr)

  res0 = sum(data*data)

  res = 0.

  mres0 = sqrt(sum(mm0*mm0))

  ndata = nt*nx
  nmodl = nt*ns

!  step=0.000005
!  alpha=790.635

  do iter=1, niter

     call velxf( 0, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, mm,z2)
     rr = rr-data

     res(iter) = sum(rr*rr)/res0

     do i=1, ndata
        trr(i)=trr(i)-step*rr(i)
     end do

     call velxf(1,0,nt,dt,t0,x2,nx,trr,mask,s,ns,dm,z2)

     call shrink(mm,dm,nmodl,1.)

     do i=1,nmodl
        mm(i)=alpha*mm(i)
     end do

     mres(iter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

     write(0,*)'iter=',iter,'step=',step,'mres=',mres(iter),'rr=',res(iter)

   end do
end subroutine velinvlbreg

subroutine velinvbbls( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,niter,step0,alpha, nres,mm0,mres)
! linearized Bregman iteration with Barzilai-Borwen (BB) stepsize and nonmonotone line search
! solve:
!        min ||m||_1+1/(2*alpha)*||x||_2, subject to Ax=b
! at each iteration...
! trr(i+1)=trr(i)-step*[Ax(i)-b]
! x(i+1)=alpha*shrink(A'trr,1.)

  integer nt,it,i, nx,ix, ns,is,niter,iter
  integer ndata, nmodl
  real data(nt*nx), mm(nt*ns), nres(niter), mm0(nt*ns), mres(niter)
  real dt,t0, dot,step0,step,alpha, isobj, isq, isqp, isrho, isconst, obj,res0,mres0
  real x2(nx), s(ns), z2(nt), mask(nx)
  real res(nt*nx), resp(nt*nx), y(nt*nx), yp(nt*nx), y_diff(nt*nx)
  real dm(nt*ns), rm(nt*ns), atb(nt*ns), dmp(nt*ns)

!  step=0.000005
!  alpha=790.635

  ndata = nt*nx
  nmodl = nt*ns

  y = 0.
  dm = 0.
  res = data
  nres = 0.
  nres(1)=1.
  res0 = sum(res*res)

  mres0 = sqrt(sum(mm0*mm0))

  isobj = 0.
  isq = 1.


  do iter=1, niter

     !-------------------BB step-------------------------------
     if (iter>1) then
        y_diff = y - yp
        step = sum(y_diff*y_diff)/sum(y_diff*(resp-res))
        if(isnan(step)) step=step0
     else
        call velxf(1,0,nt,dt,t0, x2,nx,data,mask, s,ns,atb,z2)
        step = step0+1./maxval(abs(atb))
     endif

     !--------------------y update----------------------------
     yp = y
     dmp = dm                ! dm<-ATy; mm<-x
     y = y + step*res
     call velxf(1,0,nt,dt,t0,x2,nx,y,mask,s,ns,dm,z2)
     rm = dm - dmp           ! rm<-ATy_diff

     !--------------------m update----------------------------
     call shrink(mm,dm,nmodl,1.)
     mm = alpha*mm

     obj = sum(data*y)-sum(mm*mm)/alpha/2.

     !--------------------nonmonotone line search -------------
     isrho = 1.
     isconst = 1e-3*step*(sum(res*res))
     do while(obj<isobj+isrho*isconst)
       isrho = isrho*0.5
       y = yp + (isrho*step)*res
       dm = dmp + isrho*rm
       call shrink(mm,dm,nmodl,1.)
       mm = alpha*mm
       obj = sum(data*y)-sum(mm*mm)/alpha/2.
     enddo

     !-------------------update averaged objective-isobj-for nex line search--
     isqp = isq
     isq = 0.85*isqp+1.
     isobj = (0.85*isqp*isobj+obj)/isq

     !-------------------update other quantities--------------------------
     resp = res
     call velxf( 0, 0, nt,dt,t0, x2,nx, res,mask, s,ns, mm,z2)
     res = data-res
 
     nres(iter) = sum(res*res)/res0

     mres(iter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

     write(0,*)'iter=',iter,'step=',step,'mres=',mres(iter),'rr=',nres(iter),sum(res*res),'obj=',obj

   end do
end subroutine velinvbbls

subroutine velinvaccel( nt,dt,t0, x2,nx,data, s,ns,mm, z2,mask,niter,lip,alpha,reset, nres,mm0,mres)
! linearized Bregman iteration with Nesterov's acceleration and reset (restart or skip)
! solve:
!        min ||m||_1+1/(2*alpha)*||x||_2^2, subject to Ax=b
! at each iteration...
! trr(i+1)=trr(i)-step*[Ax(i)-b]
! x(i+1)=alpha*shrink(A'trr,1.)

  integer nt,it,i, nx,ix, ns,is,niter,iter
  integer ndata, nmodl, reset
  real data(nt*nx), mm(nt*ns), nres(niter), mm0(nt*ns),mres(niter)
  real dt,t0,alpha, theta, beta, lip, res0, mres0
  real x2(nx), s(ns), z2(nt), mask(nx), obj(niter), obj_test
  real res(nt*nx), z(nt*nx), y(nt*nx), z_new(nt*nx)
  real dm(nt*ns)

!  step=0.000005
!  alpha=790.635

  ndata = nt*nx
  nmodl = nt*ns

  y = 0.
  z = 0.
  res = data
  res0 = sum(res*res)
  nres = 0.
  nres(1) =1.

  mres0 = sqrt(sum(mm0*mm0))

  theta = 1.

  do iter=1, niter

    beta = (1-theta)*(sqrt(theta*theta+4.)-theta)/2.

    !  --- restart or skip extrapolation ---
    select case (reset)
      case (1)
        if (iter>2.and.obj(iter-1)<obj(iter-2)) then
           beta=0.
           theta=1.
        endif
      case (2)
        if (iter>2.and.obj(iter-1)<obj(iter-2)) beta=0.
      case (3)
        if (iter>2.and.obj(iter-1)<1.5*obj(iter-2)) then
           beta=0.
           theta=1.
        endif
      case (4)
        if (iter>2.and.obj(iter-1)<1.5*obj(iter-2)) beta=0.
      case (5)
        if (iter>2.and.(sum(res*(y+res/lip-z))<0)) then
           beta=0.
           theta=1.
        endif
      case (6)
        if (iter>2.and.(sum(res*(y+res/lip-z))<0)) beta=0.
    end select

    ! --------------y update -------------------
    z_new = y + res/lip
    y = z_new + beta*(z_new-z)
    z = z_new

    ! -------------x update -------------------
    call velxf(1,0,nt,dt,t0, x2,nx,y,mask, s,ns,dm,z2)
    call shrink(mm,dm,nmodl,1.)
    mm = alpha*mm
    call velxf( 0, 0, nt,dt,t0, x2,nx, res,mask, s,ns, mm,z2)
    res = data-res

    ! -------------theta update--------------------
    theta = theta*(sqrt(theta*theta+4.)-theta)/2.

    ! -------------diagnostics---------------------
    obj(iter) = sum(data*y)-sum(mm*mm)/alpha/2.
    ! obj_test = sum(dprod(data,y))-sum(dprod(mm,mm))/alpha/2.

    nres(iter) = sum(res*res)/res0

    mres(iter) = sqrt(sum((mm-mm0)*(mm-mm0)))/mres0

    write(0,*)'iter=',iter,'mres=',mres(iter),'rr=',nres(iter),'obj=',obj(iter)

   end do
end subroutine velinvaccel


subroutine shrink(mm, rm, nm, delta)
  integer i,nm
  real delta
  real mm(nm),rm(nm)

  mm = 0.
  do i=1, nm
     if (rm(i) > delta) then 
        mm(i) = rm(i)-delta
     elseif (rm(i) < -delta) then
        mm(i) = rm(i)+delta
     else 
        mm(i)=0.
     end if
  end do
end subroutine shrink

subroutine powwt(iter, rr,n, wt, wrr, pow,eps,huber,srate)
  integer iter, n, i,huber
  real rr(n)     !     input vector
  real wt(n)     ! weighting vector
  real wrr(n)    !  weighted vector
  real pow       ! weighting power
  real eps,small,minval,maxval,srate
  real arr(n), awt(n)

  if ( pow.eq.0. .or. iter.eq.1 ) then
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
  if ( iter .eq. 1 ) then
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
