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
  !----------------------------------------------------------------
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
  subroutine hinaxis(inp,axis,i)
    type(file), intent(in)   :: inp
    integer    , intent(in)  :: i
    type(myaxis), intent(out) :: axis
    character (len=128)      :: buffer

    write (buffer, "(a,i1)" ) 'n',i
    call from_par(inp,buffer,axis%n,1)
    write (buffer, "(a,i1)" ) 'o',i
    call from_par(inp,buffer,axis%o,0.)
    write (buffer, "(a,i1)" ) 'd',i
    call from_par(inp,buffer,axis%d,1.)

    write (buffer, "(a,i1)") 'label',i
    call from_par(inp,buffer,axis%l," ")
  end subroutine  hinaxis
  !----------------------------------------------------------------
  subroutine houaxis(outp,axis,i)
    type(file), intent(in)   :: outp
    integer    , intent(in)  :: i
    type(myaxis), intent(in)  :: axis
    character (len=128)      :: buffer
    write (buffer, "(a,i1)" ) 'n',i
    call to_par(outp,buffer,axis%n)
    write (buffer, "(a,i1)" ) 'o',i
    call to_par(outp,buffer,axis%o)
    write (buffer, "(a,i1)" ) 'd',i
    call to_par(outp,buffer,axis%d)

    write (buffer, "(a,i1)") 'label',i
    call to_par(outp,buffer,axis%l)
  end subroutine houaxis
  !----------------------------------------------------------------
  subroutine pinaxis(axis)
    type(myaxis), intent(out) :: axis
    character (len=128)      :: buffer

    write (buffer, "(a,a)" ) trim(axis%l),'n'
    call from_par(trim(buffer),axis%n,1)

    write (buffer, "(a,a)" ) trim(axis%l),'o'
    call from_par(buffer,axis%o,0.)

    write (buffer, "(a,a)" ) trim(axis%l),'d'
    call from_par(buffer,axis%d,1.)
  end subroutine  pinaxis
  !----------------------------------------------------------------
  subroutine ainaxis(axis,i,fil)
    type(file)                :: fil
    integer    , intent(in)   :: i
    type(myaxis), intent(out) :: axis
    character (len=128)       :: buffer

    write (buffer, "(a,i1)" ) 'n',i
    call from_par(fil,buffer,axis%n)
    write (buffer, "(a,i1)" ) 'o',i
    call from_par(fil,buffer,axis%o)
    write (buffer, "(a,i1)" ) 'd',i
    call from_par(fil,buffer,axis%d)
  end subroutine ainaxis
  !----------------------------------------------------------------
  !! write axis info to an aux file
  subroutine aouaxis(axis,i,fil)
    type(file),       intent(in)  :: fil
    integer,          intent(in)  :: i
    type(myaxis),     intent(in)  :: axis
    character(len=128)            :: buffer

    write (buffer, "(a,i1)") 'n',i
    call to_par(fil,trim(buffer),axis%n)
    write (buffer, "(a,i1)") 'o',i
    call to_par(fil,trim(buffer),axis%o)
    write (buffer, "(a,i1)") 'd',i
    call to_par(fil,trim(buffer),axis%d)
  end subroutine aouaxis
  !----------------------------------------------------------------

  !  subroutine printaxis(axis)
  !    type(myaxis), intent(out) :: axis
  !    write(0,*) axis%n,axis%o,axis%d,trim(axis%l)
  !  end subroutine printaxis

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
