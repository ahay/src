module fft
  use rsf

  implicit none

contains

  subroutine ftu (signi, cx, scale)
    !   complex fourier transform with unitary scaling
    !
    !               1         nx          signi*2*pi*i*(j-1)*(k-1)/nx
    !   cx(k)  =  -------- * sum cx(j) * e
    !             sqrt(nx)   j=1             for k=1,2,...,nx=2**integer
    !
    complex, dimension (:)         :: cx
    logical, intent (in)           :: signi
    real,    intent (in), optional :: scale

    integer :: nx, i, j, k, m, istep
    real    :: arg
    complex :: cw, cdel, ct

    nx = size (cx) 
    if (nx /= pad2(nx)) call sf_error('ftu: nx not a power of 2')
    cx = cx / sqrt( 1.*nx)

    j = 1; k = 1
    do i= 1, nx 
       if (i<=j) then
          ct = cx(j); cx(j) = cx(i); cx(i) = ct 
       end if
       m = nx/2
       do while (j>m .and. m>1) 
          j = j-m; m = m/2 
       end do
       j = j+m
    end do
    do
       istep = 2*k;   cw = 1.;   arg = 3.14159265/k
       if (signi) arg = - arg
       cdel = cmplx(cos(arg), sin(arg))
       do m= 1, k 
          do i= m, nx, istep
             ct=cw*cx(i+k);  cx(i+k)=cx(i)-ct;  cx(i)=cx(i)+ct 
          end do
          cw = cw * cdel
       end do
       k = istep
       if (k>=nx) exit
    end do
    if (present (scale)) then
       if (signi) then 
          cx = cx/scale
       else
          cx = cx*scale
       end if
    end if
  end subroutine ftu

  function pad2 (n) result (n2)
    integer :: n2
    integer, intent (in) :: n
    n2 = 1
    do while( n2 < n )
       n2 = n2 * 2 
    end do
  end function pad2

  ! FT a vector in a matrix, with first omega = - pi
  subroutine  fth (adj, signi, cx)
    logical, intent (in)   :: adj, signi
    complex, dimension (:) :: cx

    if (adj) then
       call ftu (.not. signi, cx)
       cx (2::2) = - cx (2::2)
    else
       cx (2::2) = - cx (2::2)
       call ftu (signi, cx)
    end if
  end subroutine fth

  ! 1D Fourier transform on a 2D data set along the 1-axis
  !
  subroutine  ft1axis( adj, sign1, cx)
    logical, intent (in)     :: adj, sign1
    complex, dimension (:,:) :: cx
 
    integer                  :: i2
    do i2= 1, size (cx,2)
       call fth( adj, sign1, cx(:,i2))
    end do
  end subroutine ft1axis

  ! 1D Fourier transform on a 2D data set along the 2-axis
  !
  subroutine  ft2axis( adj, sign2, cx)
    logical, intent (in)     :: adj, sign2
    complex, dimension (:,:) :: cx
 
    integer                  :: i1
    do i1= 1, size (cx,1)
       call fth( adj, sign2, cx(i1,:))
    end do
  end subroutine ft2axis

end module fft

