module cub
  use rsf

  implicit none

  type, public :: axis
     integer   :: n
     real      :: o,d
  end type axis

contains

!------------------------------------------------------------

  subroutine iaxis(FF,AA,i)
    type(file), intent(in) :: FF
    integer   , intent(in) :: i
    type(axis), intent(out):: AA
    character(len=128)     :: BB

    write(BB,"(a,i1)" ) 'n',i
    call from_par(FF,BB,AA%n,1)
    write(BB,"(a,i1)" ) 'o',i
    call from_par(FF,BB,AA%o,0.)
    write(BB,"(a,i1)" ) 'd',i
    call from_par(FF,BB,AA%d,1.)

  end subroutine iaxis

!------------------------------------------------------------

  subroutine oaxis(FF,AA,i)
    type(file), intent(in) :: FF
    integer   , intent(in) :: i
    type(axis), intent(in) :: AA
    character(len=128)     :: BB

    write(BB,"(a,i1)" ) 'n',i
    call to_par(FF,BB,AA%n)
    write(BB,"(a,i1)" ) 'o',i
    call to_par(FF,BB,AA%o)
    write(BB,"(a,i1)" ) 'd',i
    call to_par(FF,BB,AA%d)

  end subroutine oaxis

!------------------------------------------------------------

end module cub
