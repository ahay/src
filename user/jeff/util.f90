module util
  use rsf
  implicit none

  type, public :: myaxis
     integer   :: n
     real      :: o,d
     character(len=128) :: l
  end type myaxis

  interface grab_param
     module procedure grab_param_real
     module procedure grab_param_logical
     module procedure grab_param_integer
  end interface

contains

  !-----------------------------------------
  integer function sig(m) result(s)
    real            :: m
    if (m>0.) then
       s=+1
    else if (m .eq. 0.) then
       s=0
    else 

       s=-1
    end if
  end function sig
  !----------------------------------------------------------------
!  subroutine hinaxis(axis,i)
!   integer    , intent(in)  :: i
!   type(myaxis), intent(out) :: axis
!   character (len=128)      :: buffer

!   write (buffer, "(a,i1)" ) 'n',i
!   call from_history(buffer,axis%n,1)
!   write (buffer, "(a,i1)" ) 'o',i
!   call from_history(buffer,axis%o,0.)
!   write (buffer, "(a,i1)" ) 'd',i
!   call from_history(buffer,axis%d,1.)

!   write (buffer, "(a,i1)") 'label',i
!   call from_history(buffer,axis%l," ")
!write(0,*) "hin",axis%n,axis%o,axis%d,axis%l
! end subroutine  hinaxis
! !----------------------------------------------------------------
! subroutine houaxis(axis,i)
!   integer    , intent(in)  :: i
!   type(myaxis),intent(in)  :: axis
!   character (len=128)      :: buffer
!   write (buffer, "(a,i1)" ) 'n',i
!   call to_history(buffer,axis%n)
!   write (buffer, "(a,i1)" ) 'o',i
!   call to_history(buffer,axis%o)
!   write (buffer, "(a,i1)" ) 'd',i
!   call to_history(buffer,axis%d)

!   write (buffer, "(a,i1)") 'label',i
!   call to_history(buffer,axis%l)
! end subroutine houaxis
! !----------------------------------------------------------------
! subroutine pinaxis(axis)
!   type(myaxis), intent(out) :: axis
!   character (len=128)       :: buffer

!   write (buffer, "(a,a)" ) trim(axis%l),'n'
!   call from_param(trim(buffer),axis%n,1)

!   write (buffer, "(a,a)" ) trim(axis%l),'o'
!   call from_param(buffer,axis%o,0.)

!   write (buffer, "(a,a)" ) trim(axis%l),'d'
!   call from_param(buffer,axis%d,1.)
! end subroutine  pinaxis
! !----------------------------------------------------------------
! subroutine ainaxis(axis,i,file)
!   character(len=128)       :: file
!   integer    , intent(in)  :: i
!   type(myaxis),intent(out) :: axis
!   character (len=128)      :: buffer

!   write (buffer, "(a,i1)" ) 'n',i
!   call from_aux(file,buffer,axis%n)
!   write (buffer, "(a,i1)" ) 'o',i
!   call from_aux(file,buffer,axis%o)
!   write (buffer, "(a,i1)" ) 'd',i
!   call from_aux(file,buffer,axis%d)
!write(0,*) "ain",axis%n,axis%o,axis%d
! end subroutine ainaxis
  !----------------------------------------------------------------
  subroutine printaxis(axis)
    type(myaxis), intent(in) :: axis
    write(0,*) axis%n,axis%o,axis%d,trim(axis%l)
  end subroutine printaxis
  !----------------------------------------------------------------
  subroutine grab_param_real(pname,pval,pdefault)
    real,               intent(inout) :: pval
    character(len=128), intent(in)    :: pname
    real,               intent(in)    :: pdefault

    call from_par(pname,pval,pdefault)
    write(0,*) trim(pname),pval
  end subroutine grab_param_real


  subroutine grab_param_logical(pname,pval,pdefault)
    logical,            intent(inout) :: pval
    character(len=128), intent(in)    :: pname
    logical,            intent(in)    :: pdefault

    call from_par(pname,pval,pdefault)
    write(0,*) trim(pname),pval
  end subroutine grab_param_logical


  subroutine grab_param_integer(pname,pval,pdefault)
    integer,            intent(inout) :: pval
    character(len=128), intent(in)    :: pname
    integer,            intent(in)    :: pdefault

    call from_par(pname,pval,pdefault)
    write(0,*) trim(pname),pval
  end subroutine grab_param_integer

end module util













