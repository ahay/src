! 2-D finite-difference acoustic wave propagation 
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

program Wave
  use rsf ! use module
  use laplace ! use module
  implicit none
  
  integer :: it, nt, ft, jt, n1, n2
  real    :: dt, d1, d2, dt2

  real, dimension (:),   allocatable :: ww
  real, dimension (:,:), allocatable :: vv, rr
  real, dimension (:,:), allocatable :: u0,u1,u2,ud
  
  type(file) :: Fw,Fv,Fr,Fo  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fr = rsf_input ("in")   ! source position 
  Fo = rsf_output("out")  ! output wavefield 

  Fw = rsf_input ("wav")  ! source wavelet 
  Fv = rsf_input ("v")    ! velocity 

  ! Read/Write axes
  call from_par(Fr,"n1",n1)
  call from_par(Fr,"n2",n2)
  call from_par(Fr,"d1",d1)
  call from_par(Fr,"d2",d2)

  call from_par(Fw,"n1",nt)
  call from_par(Fw,"d1",dt)
  
  call from_par("ft",ft,0) ! first recorded time
  call from_par("jt",jt,0) ! time interval

  ! set dimension of output file
  call to_par(Fo,"n3",(nt-ft)/jt)
  call to_par(Fo,"d3",jt*dt)
  call to_par(Fo,"o3",ft*dt)

  dt2 = dt*dt
  ft  = ft+1

  ! set Laplacian coefficients 
  call laplacian_set(1.0/(d1*d1),1.0/(d2*d2))
  
  ! read wavelet, velocity & source position 
  allocate (ww(nt), vv(n1,n2), rr(n1,n2))
  call rsf_read(Fw,ww)
  call rsf_read(Fv,vv)
  call rsf_read(Fr,rr)

  ! allocate temporary arrays */
  allocate (u0(n1,n2), u1(n1,n2), u2(n1,n2), ud(n1,n2))
    
  ! initialize temporary arrays
  u0=0.0
  u1=0.0
  u2=0.0
  ud=0.0
  vv = vv*vv*dt2

  ! Time loop 
  do it=1,nt

     ! call function
     call laplacian(u1,ud)

     ! scale by velocity 
     ud = ud*vv

     ! inject wavelet 
     ud = ud + ww(it) * rr
     
     ! time step 
     u2 = 2*u1 - u0 + ud 
     u0 = u1
     u1 = u2
	
     ! write wavefield to output 
     if (it >= ft .and. 0 == mod(it-ft,jt)) then
        write(0,'(a,i4)',advance='no') achar(13), it 
        call rsf_write(Fo,u1)
     end if
  end do
  write (0,*) 
    
  call exit(0)
end program Wave

