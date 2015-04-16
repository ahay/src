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

