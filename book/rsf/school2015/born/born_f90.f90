! 2-D finite-difference born modeling 

module laplace
  ! Laplacian operator, 4th-order finite-difference 
  implicit none
  real :: c0, c11, c21, c12, c22  
contains
  subroutine laplacian_set(d1,d2)
    real, intent (in) :: d1,d2
    
    c11 = 4.0*d1/3.0
    c12=  -d1/12.0
    c21 = 4.0*d2/3.0
    c22=  -d2/12.0
    c0  = -2.0 * (c11+c12+c21+c22)
  end subroutine laplacian_set

  subroutine laplacian(uin, uout)
    real, dimension (:,:), intent (in)  :: uin
    real, dimension (:,:), intent (out) :: uout
    integer n1, n2
    
    n1 = size(uin,1)
    n2 = size(uin,2)
    
    uout(3:n1-2,3:n2-2) = &
         c11*(uin(2:n1-3,3:n2-2) + uin(4:n1-1,3:n2-2)) + &
         c12*(uin(1:n1-4,3:n2-2) + uin(5:n1,  3:n2-2)) + &
         c21*(uin(3:n1-2,2:n2-3) + uin(3:n1-2,4:n2-1)) + &
         c22*(uin(3:n1-2,1:n2-4) + uin(3:n1-2,5:n2  )) + &
         c0*uin(3:n1-2,3:n2-2)
  end subroutine laplacian
end module laplace

subroutine deriv2(ww, nt)
! second time derivative
integer it, nt
real ww(nt), temp(nt)

temp(1)=ww(1)
do it=2, nt
    temp(it) = ww(it)-ww(it-1)
end do
do it=1, nt-1
    ww(it) = temp(it+1)-temp(it)
end do
ww(nt)=temp(nt)-temp(nt-1)
return
end subroutine deriv2

program Born_modeling 
  use rsf ! use module
  use laplace ! use module
  implicit none
  
  logical :: born
  integer :: it, nt, ft, jt, n1, n2
  real    :: dt, d1, d2, dt2
  integer :: nr, r0, rz ! receiver setup

  real, dimension (:),   allocatable   :: ww
  real, dimension (:,:), allocatable   :: vv,rr,ref,dd
  real, dimension (:,:), allocatable   :: u0,u1,u2,ud,vv2
  real, dimension (:,:,:), allocatable  :: wave

  type(file) :: Fw,Fv,Fr,Fo,Fd,Fdetv  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! born modeling or forward modeling
  call from_par("born", born, .true.)

  ! setup I/O files 
  Fr = rsf_input ("in")   ! source position 
  Fd = rsf_output("out")  ! output shot record

  Fo = rsf_output("snapshot") ! output wavefield
  Fw = rsf_input ("wav")  ! source wavelet 
  Fv = rsf_input ("v")    ! velocity 
  
  ! velocity perturbation/reflectivity
  if(born) then
      Fdetv = rsf_input("detv")
  end if
  
  ! Read/Write axes
  call from_par(Fr,"n1",n1)
  call from_par(Fr,"n2",n2)
  call from_par(Fr,"d1",d1)
  call from_par(Fr,"d2",d2)

  call from_par(Fw,"n1",nt)
  call from_par(Fw,"d1",dt)
  
  call from_par("nr",nr,n2) ! trace number
  call from_par("r0",r0,0) ! starting position
  call from_par("rz",rz,0) ! depth of shot record
  call from_par("ft",ft,0) ! first recorded time
  call from_par("jt",jt,1) ! time interval

  ! set dimension of output data file
  call to_par(Fd,"n1",nt)
  call to_par(Fd,"d1",dt)
  call to_par(Fd,"o1",0.)
  call to_par(Fd,"label1","Time")
  call to_par(Fd,"unit1","s")

  call to_par(Fd,"n2",nr)
  call to_par(Fd,"d2",d2)
  call to_par(Fd,"o2",r0*d2)

  ! set dimension of output wavefield file
  call to_par(Fo,"n3",(nt-ft)/jt)
  call to_par(Fo,"d3",jt*dt)
  call to_par(Fo,"o3",ft*dt)

  dt2 = dt*dt
  ft  = ft+1

  ! set Laplacian coefficients 
  call laplacian_set(1.0/(d1*d1),1.0/(d2*d2))
  
  ! read wavelet, velocity, source position & reflectivity
  allocate (ww(nt), vv(n1,n2), rr(n1,n2))
  call rsf_read(Fw,ww)
  call rsf_read(Fv,vv)
  call rsf_read(Fr,rr)
  if(born) then
      allocate (ref(n1,n2))
      call rsf_read(Fdetv,ref)
  end if

  ! allocate wavefield and data arrays
  allocate (wave(n1,n2,nt), dd(nt,nr), vv2(n1,n2))

  ! allocate temporary arrays */
  allocate (u0(n1,n2), u1(n1,n2), u2(n1,n2), ud(n1,n2))
    
  ! initialize temporary arrays
  u0=0.0
  u1=0.0
  u2=0.0
  ud=0.0
  vv2 = vv*vv*dt2

  ! second time derivative
  if(born) then
      call deriv2(ww, nt)
  end if

  ! Time loop 1 (U_0)
  do it=1,nt
      write(0,'(a,i4)',advance='no') achar(13), it
      
      ! wavefield storage
      wave(:,:,it) = u1
      
      ! call function
      call laplacian(u1,ud)
      
      ! scale by velocity 
      ud = ud*vv2

      ! inject wavelet 
      ud = ud + ww(it) * rr
      
      ! time step 
      u2 = 2*u1 - u0 + ud 
      u0 = u1
      u1 = u2
  end do ! end of it

  ! if forward modeling, solve wave equation once
  if(.not. born) then
      do it=1,nt
          ! write shot record to output
          dd(it,1:nr) = wave(rz,r0:r0+nr,it)
          
          ! write wavefield to output 
          if (it >= ft .and. 0 == mod(it-ft,jt)) then
              call rsf_write(Fo,wave(:,:,it))
          end if
      end do

      call rsf_write(Fd, dd)

      ! if forward modeling, only solve wave equation once
      call exit(0)
  end if

  ! second initialization
  u0=0.0
  u1=0.0
  u2=0.0
  ud=0.0

  ! Time loop 2 (det U)
  do it=1,nt
      write(0,'(a,i4)',advance='no') achar(13), it
      
      ! write shot record to output
      dd(it,1:nr) = u1(rz,r0:r0+nr)
      
      ! write wavefield to output 
      if (it >= ft .and. 0 == mod(it-ft,jt)) then
          call rsf_write(Fo, u1)
      end if
      
      ! call function
      call laplacian(u1,ud)
      
      ! scale by velocity 
      ud = ud*vv2

      ! inject source term
      ud = ud + wave(:,:,it)*ref*2./vv
      
      ! time step 
      u2 = 2*u1 - u0 + ud 
      u0 = u1
      u1 = u2
  
  end do ! end of it

  call rsf_write(Fd, dd)
    
  call exit(0)
end program Born_modeling 
