module helixcartmod
! vector to matrix index transform and its inverse

  implicit none
  contains
  subroutine helix2cart( nn, i, ii)
    integer, dimension( :), intent( in) :: nn
! cartesian axes (n1,n2,n3,...)

    integer, dimension( :), intent(out) :: ii
! cartesn coords (i1,i2,i3,...)

    integer               , intent( in) :: i
! equivalent 1-D linear index

    integer                             :: axis, n123
    n123 = 1
    do axis = 1, size(  nn)
      ii( axis) = mod( ( i-1)/n123, nn( axis)) + 1
      n123 = n123 * nn( axis)
    end do
  end subroutine
  subroutine cart2helix( nn, ii, i)
    integer, dimension( :), intent( in) :: nn, ii
    integer                             :: i, axis, n123
    n123 = 1
    i = 1
    do axis = 1, size(  nn)
      i = i + ( ii( axis)-1)*n123
      n123 = n123 * nn( axis)
    end do
  end subroutine
end module
