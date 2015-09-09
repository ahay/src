! Fortran-90 interface

module RSF
  implicit none

  ! File below to be created by config. Defines integer PTRKIND (kind of
  ! integer used for pointers. Necessary to do it this way because F90
  ! cannot figure out if platform is 32-bit or 62-bit
  include "ptr_sz.f90"

  ! Kind of integer for representing positions (offsets) inside a file
  ! Equivalent in C defined in fortran.c
  integer, parameter :: OFFKIND=8

  ! External string functions cannot return LEN=* strings
  ! Choose something large enough for most functions
  integer, parameter :: FSTRLEN=256

  ! Types of data in RSF binary. Correspond to sf_datatype in filt/lib/file.c
  integer, parameter :: sf_uchar=0, sf_char   =1, sf_int=2
  integer, parameter :: sf_float=3, sf_complex=4, sf_short=5

  ! For calling sf_seek. For portability, these should actually be extracted
  ! from stdio.h during the configure step of the installation
  integer, parameter :: sf_seek_set=0, sf_seek_cur=1, sf_seek_end=2

  type file
     private
     integer(kind=PTRKIND) :: tag
  end type file

  type, public :: axa
     integer   :: n
     real      :: o,d
  end type axa

  interface from_par
     module procedure from_history_int
     module procedure from_history_int_array
     module procedure from_history_real
     module procedure from_history_string
!!$     module procedure from_history_real_array
!!$     module procedure from_history_char
     module procedure from_history_dim
     module procedure from_param_int
     module procedure from_param_int_array
     module procedure from_param_real
     module procedure from_param_real_array
     module procedure from_param_bool
     module procedure from_param_bool_array
     module procedure from_param_string
  end interface

  interface from_either
     module procedure from_either_int
     module procedure from_either_real
     module procedure from_either_dim
!!$     module procedure from_either_bool
!!$     module procedure from_either_char
  end interface

  interface to_par
     module procedure to_history_int
     module procedure to_history_int_array
     module procedure to_history_real
     module procedure to_history_string
!!$     module procedure to_history_real_array
!!$     module procedure to_history_bool
!!$     module procedure to_history_char
  end interface

  interface rsf_read
     module procedure rsf_read_1d_n
     module procedure rsf_read_1d
     module procedure rsf_read_2d
     module procedure rsf_read_3d
     module procedure rsf_read_4d
     module procedure rsf_read_5d
     module procedure rsf_read_complex_1d
     module procedure rsf_read_complex_2d 
     module procedure rsf_read_complex_3d 
     module procedure rsf_read_complex_4d 
     module procedure rsf_read_complex_5d
  end interface

  interface rsf_write
     module procedure rsf_write_1d_n
     module procedure rsf_write_1d
     module procedure rsf_write_2d
     module procedure rsf_write_3d
     module procedure rsf_write_4d
     module procedure rsf_write_5d
     module procedure rsf_write_complex_1d
     module procedure rsf_write_complex_2d 
     module procedure rsf_write_complex_3d 
     module procedure rsf_write_complex_4d 
     module procedure rsf_write_complex_5d 
  end interface

contains
  function gettype (f) result (t)
    integer                  :: t
    type (file), intent (in) :: f

    integer sf_gettype
    external sf_gettype

    t = sf_gettype(f%tag)
  end function gettype

  subroutine settype (f,t)
    integer,     intent (in) :: t
    type (file), intent (in) :: f

    call sf_settype(f%tag,t)
  end subroutine settype

  function filesize (f, dim) result (s)
    integer                  :: s 
    type (file), intent (in) :: f
    integer,     intent (in) :: dim
    optional                 :: dim

    integer(kind=OFFKIND) sf_filesize, sf_leftsize
    external sf_filesize, sf_leftsize

    if (present (dim)) then
       s = sf_leftsize(f%tag,dim)
    else
       s = sf_filesize(f%tag)
    end if
  end function filesize

  function rsf_input (tag) result (f)
    type (file)                    :: f
    character (len=*), intent (in) :: tag
    optional                       :: tag

    integer(kind=PTRKIND) sf_input
    external sf_input

    if (present (tag)) then
       f%tag = sf_input(tag)
    else 
       f%tag = sf_input("in")
    end if
  end function rsf_input

  function rsf_output (tag) result (f)
    type (file)                    :: f
    character (len=*), intent (in) :: tag
    optional                       :: tag

    integer(kind=PTRKIND) sf_output
    external sf_output

    if (present (tag)) then
       f%tag = sf_output(tag)
    else 
       f%tag = sf_output("out")
    end if
  end function rsf_output

  function dimension (hist, n) result (nd)
    integer                              :: nd
    integer, dimension (:), intent (out) :: n
    type (file),            intent (in)  :: hist
    integer sf_filedims
    external sf_filedims

    nd = sf_filedims(hist%tag,n)
  end function dimension
  
  function axisname (i,var) result (name)
    character (len=3) :: name
    integer           :: i
    character (len=*) :: var
    optional          :: var
    if (present (var)) then
       if (i < 10) then
          write (name,"(a,i1)") var,i
       else
          write (name,"(a,i2)") var,i
       end if
    else
       if (i < 10) then
          write (name,"(a,i1)") "n",i
       else
          write (name,"(a,i2)") "n",i
       end if
    end if
  end function axisname

  subroutine from_history_dim (hist, n1, n2, n3, n4, type)
    integer,             intent (out) :: n1, n2, n3, n4
    type (file),         intent (in)  :: hist
    character (len = *), intent (in)  :: type    
    optional                          :: n2, n3, n4, type

    integer     :: ftype

    integer sf_gettype
    external sf_gettype

    call from_history_int (hist,"n1",n1)
    if (present (n2)) call from_history_int (hist,"n2",n2,1)
    if (present (n3)) call from_history_int (hist,"n3",n3,1)
    if (present (n4)) call from_history_int (hist,"n4",n4,1)
    
    ftype = sf_gettype(hist%tag)
    if (present (type)) then
       select case (type(1:1))
       case ("f")
          if (ftype /= 3) call sf_error("Need float input")
       case ("i")
          if (ftype /= 2) call sf_error("Need int input")
       case ("c")
          if (ftype /= 4) call sf_error("Need complex input")
       end select
    else if (ftype /= 3) then
       call sf_error("Need float input")
    end if
  end subroutine from_history_dim

  subroutine from_history_int (hist, name, value, default)
    type (file),         intent (in) :: hist
    character (len = *), intent (in) :: name
    integer, intent (out)            :: value
    integer, intent (in), optional   :: default

    logical sf_histint
    external sf_histint

    if(.not. sf_histint(hist%tag,name,value)) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing history value: " // name)
       end if
    end if
  end subroutine from_history_int

!!$  subroutine from_history_char (hist, name, value, default)
!!$    type (file),         intent (in)  :: hist
!!$    character (len = *), intent (in)  :: name, default
!!$    character (len = *), intent (out) :: value
!!$    optional                          :: default
!!$
!!$    character (len = *) sf_histtsring
!!$    external sf_histstring
!!$
!!$    if(hetch(name,'s',value) == 0) then
!!$       if (present (default)) then
!!$          value = default
!!$       else
!!$          call sf_error("missing history value: " // name)
!!$       end if
!!$    end if
!!$  end  subroutine from_history_char

  subroutine from_history_int_array (hist, name, value, default)
    type (file),         intent (in)               :: hist
    character (len = *), intent (in)               :: name
    integer, intent (out), dimension (:)           :: value
    integer, intent (in),  dimension (:), optional :: default

    logical sf_histints
    external sf_histints

    if(.not. sf_histints(hist%tag,name,value,size(value))) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing history value: " // name)
       end if
    end if
  end subroutine from_history_int_array

  subroutine from_history_real (hist, name, value, default)
    type (file),         intent (in)  :: hist
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default

    logical sf_histfloat
    external sf_histfloat

    if(.not. sf_histfloat(hist%tag,name,value)) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing history value: " // name)
       end if
    end if
  end subroutine from_history_real

  subroutine from_history_string (hist, name, value, default)
    type (file),         intent (in)  :: hist
    character (len = *), intent (in)  :: name
    character (len = *), intent (out) :: value
    character (len = *), intent (in), optional :: default

    character(len=FSTRLEN) sf_histstring
    external sf_histstring
    
    if (present (default)) value = default
    value = sf_histstring(hist%tag,name)
  end subroutine from_history_string

!!$  subroutine from_history_real_array (hist,name, value, default)
!!$    type (file),         intent (in)            :: hist
!!$    character (len = *), intent (in)            :: name
!!$    real, intent (out), dimension (:)           :: value
!!$    real, intent (in),  dimension (:), optional :: default
!!$
!!$    logical sf_histfloats
!!$    external sf_histfloats
!!$
!!$    if(.not. sf_histfloats(hist%tag,name,value,size(value))) then
!!$       if (present (default)) then
!!$          value = default
!!$       else
!!$          call sf_error("missing history value: " // name)
!!$       end if
!!$    end if
!!$  end subroutine from_history_real_array

  subroutine from_either_dim (hist, n1, n2, n3, type)
    type (file),         intent (in)  :: hist
    integer,             intent (out) :: n1, n2, n3
    character (len = *), intent (in)  :: type
    optional                          :: n2, n3, type

    integer                 :: ftype

    integer sf_gettype
    external sf_gettype

    call from_either (hist,"n1",n1)
    if (present (n2)) call from_either (hist,"n2",n2)
    if (present (n3)) call from_either (hist,"n3",n3)
    
    ftype = sf_gettype(hist%tag)
    if (present (type)) then
       select case (type(1:1))
       case ("f")
          if (ftype /= 2) call sf_error("Need float input")
       case ("i")
          if (ftype /= 1) call sf_error("Need int input")
       case ("c")
          if (ftype /= 3) call sf_error("Need complex input")
       end select
    else if (ftype /= 2) then
       call sf_error("Need float input")
    end if
  end subroutine from_either_dim

  subroutine from_either_int (hist, name, value, default)
    type (file),         intent (in)  :: hist
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value
    integer, intent (in), optional    :: default

    logical sf_histint, sf_getint
    external sf_histint, sf_getint

    if(.not. sf_histint(hist%tag,name,value)) then
       if (.not. sf_getint(name,value)) then
          if (present (default)) then
             value = default
          else
             call sf_error("missing either value: " // name)
          end if
       end if
    end if
  end  subroutine from_either_int

!!$  subroutine from_either_char (name, value, default)
!!$    character (len = *), intent (in)  :: name, default
!!$    character (len = *), intent (out) :: value
!!$    optional                          :: default
!!$
!!$    integer fetch
!!$    external fetch
!!$
!!$    if(fetch(name,'s',value) == 0) then
!!$       if (present (default)) then
!!$          value = default
!!$       else
!!$          call sf_error("missing either value: " // name)
!!$       end if
!!$    end if
!!$    if(putch_no .ne. "no putch") &
!!$    call putch("From either: " // name,'s',value)
!!$  end  subroutine from_either_char

  subroutine from_either_real (hist, name, value, default)
    type (file),         intent (in)  :: hist
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default

    logical sf_histfloat, sf_getfloat
    external sf_histfloat, sf_getfloat

    if(.not. sf_histfloat(hist%tag,name,value)) then
       if (.not. sf_getfloat(name,value)) then
          if (present (default)) then
             value = default
          else
             call sf_error("missing either value: " // name)
          end if
       end if
    end if
  end subroutine from_either_real

  subroutine from_param_int (name, value, default)
    character (len = *), intent (in)  :: name
    integer, intent (out)             :: value
    integer, intent (in), optional    :: default

    logical sf_getint
    external sf_getint

    if(.not. sf_getint(name,value)) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_int

  subroutine from_param_bool (name, value, default)
    character (len = *), intent (in)  :: name
    logical, intent (out)             :: value
    logical, intent (in), optional    :: default
    logical                           :: val

    logical sf_getbool
    external sf_getbool

    if(.not. sf_getbool(name,value)) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_bool

  subroutine from_param_int_array (name, value, default)
    character (len = *),    intent (in)           :: name
    integer, dimension (:), intent (out)          :: value
    integer, dimension (:), intent (in), optional :: default

    logical sf_getints
    external sf_getints

    if(.not. sf_getints(name,value,size(value))) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_int_array

  subroutine from_param_real (name, value, default)
    character (len = *), intent (in)  :: name
    real, intent (out)                :: value
    real, intent (in), optional       :: default

    logical sf_getfloat
    external sf_getfloat

    if(.not. sf_getfloat(name,value)) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_real

  subroutine from_param_real_array (name, value, default)
    character (len = *), intent (in)           :: name
    real, dimension (:), intent (out)          :: value
    real, dimension (:), intent (in), optional :: default

    logical sf_getfloats
    external sf_getfloats

    if(.not. sf_getfloats(name,value,size(value))) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_real_array

  subroutine from_param_bool_array (name, value, default)
    character (len = *),    intent (in)           :: name
    logical, dimension (:), intent (out)          :: value
    logical, dimension (:), intent (in), optional :: default

    logical sf_getbools
    external sf_getbools

    if(.not. sf_getbools(name,value,size(value))) then
       if (present (default)) then
          value = default
       else
          call sf_error("missing parameter value: " // name)
       end if
    end if
  end subroutine from_param_bool_array

  subroutine from_param_string (name, value, default)
    character (len = *), intent (in)  :: name, default
    character (len = *), intent (out) :: value
    optional                          :: default

    character(len=FSTRLEN) sf_getstring
    external sf_getstring

    value = sf_getstring(name)
    if (value == "" .and. present (default)) value = default
  end subroutine from_param_string

  function exist_par (name) result (test)
    logical                           :: test
    character (len = *), intent (in)  :: name
    character (len = 80)              :: par

    call from_par(name,par)
    test = (par /= "")
  end function exist_par

  subroutine to_history_int (hist, name, value)
    type (file),         intent (in) :: hist
    character (len = *), intent (in) :: name
    integer,             intent (in) :: value

    call sf_putint (hist%tag,name,value)
  end subroutine to_history_int

  subroutine to_history_int_array (hist, name, value)
    type (file),            intent (in) :: hist
    character (len = *),    intent (in) :: name
    integer, dimension (:), intent (in) :: value

    call sf_putints (hist%tag,name,value,size(value))
  end subroutine to_history_int_array

  subroutine to_history_real (hist, name, value)
    type (file),         intent (in) :: hist
    character (len = *), intent (in) :: name
    real,                intent (in) :: value
   
    call sf_putfloat (hist%tag,name,value)
  end subroutine to_history_real

  subroutine to_history_string (hist, name, value)
    type (file),         intent (in) :: hist
    character (len = *), intent (in) :: name
    character (len = *), intent (in) :: value
   
    call sf_putstring (hist%tag,name,value)
  end subroutine to_history_string

!!$ subroutine to_history_real_array (name, value,tag)
!!$    character (len = *), intent (in) :: name
!!$    real, dimension (:), intent (in) :: value
!!$    character (len=*), optional                :: tag
!!$    character (len = 1024)           :: buffer,tag_out
!!$    character (len = 4)              :: fmt
!!$
!!$    if(present(tag)) then
!!$       tag_out=tag
!!$    else
!!$      tag_out="out"
!!$    end if
!!$    write (fmt, "(i4)") size (value)
!!$    write (buffer, "(a," // fmt // '(e14.6:",")' // ")") name // '=', value
!!$    call strip_blanks (buffer)
!!$    call auxputlin (tag_out,buffer)
!!$  end  subroutine to_history_real_array
!!$
!!$  subroutine to_history_bool (name, value, file)
!!$    character (len = *), intent (in)  :: file, name
!!$    logical, intent (in)              :: value
!!$    optional                          :: file
!!$
!!$    if (present (file)) then
!!$       call auxputch(name,'l',value,file)
!!$    else
!!$    if(putch_no .ne. "no putch") &
!!$       call putch (name,'l',value)
!!$    end if
!!$  end  subroutine to_history_bool
!!$
!!$  subroutine to_history_char (name, value, file)
!!$    character (len = *), intent (in)  :: file, name, value
!!$    optional                          :: file
!!$
!!$    if (present (file)) then
!!$       call auxputch(name,'s',value,file)
!!$    else
!!$    if(putch_no .ne. "no putch") &
!!$       call putch (name,'s',value)
!!$    end if
!!$  end  subroutine to_history_char

  subroutine rsf_read_1d_n (hist, array, n)
    type (file),         intent (in)  :: hist
    integer,             intent (in)  :: n
    real, dimension (n), intent (out) :: array

    call sf_floatread(array, n, hist%tag)
  end subroutine rsf_read_1d_n

  subroutine rsf_read_1d (hist, array)    
    type (file),         intent (in)  :: hist
    real, dimension (:), intent (out) :: array

    call sf_floatread(array, size(array), hist%tag)
  end subroutine rsf_read_1d

  subroutine rsf_read_2d (hist, array)    
    type (file),           intent (in)  :: hist
    real, dimension (:,:), intent (out) :: array

    call rsf_read_1d_n(hist, array, size(array))
  end subroutine rsf_read_2d

  subroutine rsf_read_3d (hist, array)    
    type (file),             intent (in)  :: hist
    real, dimension (:,:,:), intent (out) :: array
    
    call rsf_read_1d_n(hist, array, size(array))
  end subroutine rsf_read_3d

  subroutine rsf_read_4d (hist, array)    
    type (file),               intent (in)  :: hist
    real, dimension (:,:,:,:), intent (out) :: array
    
    call rsf_read_1d_n(hist, array, size(array))
  end subroutine rsf_read_4d

  subroutine rsf_read_5d (hist, array)    
    type (file),                 intent (in)  :: hist
    real, dimension (:,:,:,:,:), intent (out) :: array
    
    call rsf_read_1d_n(hist, array, size(array))
  end subroutine rsf_read_5d
  
  subroutine rsf_read_complex_1d_n (hist, array, n)
    type (file),         intent (in)  :: hist
    integer,             intent (in)  :: n
    complex, dimension (n), intent (out) :: array

    call sf_complexread(array, n, hist%tag)
  end subroutine rsf_read_complex_1d_n

  subroutine rsf_read_complex_1d (hist, array)    
    type (file),         intent (in)  :: hist
    complex, dimension (:), intent (out) :: array

    call rsf_read_complex_1d_n(hist,array,size(array))
  end subroutine rsf_read_complex_1d

  subroutine rsf_read_complex_2d (hist, array)    
    type (file),           intent (in)  :: hist
    complex, dimension (:,:), intent (out) :: array

    call rsf_read_complex_1d_n(hist,array,size(array))
  end subroutine rsf_read_complex_2d

  subroutine rsf_read_complex_3d (hist, array)    
    type (file),             intent (in)  :: hist
    complex, dimension (:,:,:), intent (out) :: array
    
    call rsf_read_complex_1d_n(hist,array,size(array))
  end subroutine rsf_read_complex_3d

  subroutine rsf_read_complex_4d (hist, array)    
    type (file),               intent (in)  :: hist
    complex, dimension (:,:,:,:), intent (out) :: array
    
    call rsf_read_complex_1d_n(hist,array,size(array))
  end subroutine rsf_read_complex_4d

  subroutine rsf_read_complex_5d (hist, array)    
    type (file),                 intent (in)  :: hist
    complex, dimension (:,:,:,:,:), intent (out) :: array
    
    call rsf_read_complex_1d_n(hist,array,size(array))
  end subroutine rsf_read_complex_5d
  
  subroutine rsf_write_1d_n (hist, array, n)
    type (file),         intent (in)  :: hist
    integer,             intent (in)  :: n
    real, dimension (n), intent (in) :: array

    call sf_floatwrite(array, n, hist%tag)
  end subroutine rsf_write_1d_n

  subroutine rsf_write_1d (hist, array)    
    type (file),         intent (in)  :: hist
    real, dimension (:), intent (in) :: array

    call rsf_write_1d_n(hist,array,size(array))
  end subroutine rsf_write_1d

  subroutine rsf_write_2d (hist, array)    
    type (file),           intent (in)  :: hist
    real, dimension (:,:), intent (in) :: array

    call rsf_write_1d_n(hist,array,size(array))
  end subroutine rsf_write_2d

  subroutine rsf_write_3d (hist, array)    
    type (file),             intent (in)  :: hist
    real, dimension (:,:,:), intent (in) :: array
    
    call rsf_write_1d_n(hist,array,size(array))
  end subroutine rsf_write_3d

  subroutine rsf_write_4d (hist, array)    
    type (file),               intent (in)  :: hist
    real, dimension (:,:,:,:), intent (in) :: array

    call rsf_write_1d_n(hist,array,size(array))
  end subroutine rsf_write_4d

  subroutine rsf_write_5d (hist, array)    
    type (file),                 intent (in)  :: hist
    real, dimension (:,:,:,:,:), intent (in) :: array
    
    call rsf_write_1d_n(hist,array,size(array))
  end subroutine rsf_write_5d
  
  subroutine rsf_write_complex_1d_n (hist, array, n)
    type (file),            intent (in)  :: hist
    integer,                intent (in)  :: n
    complex, dimension (n), intent (in) :: array

    call sf_complexwrite(array, n, hist%tag)
  end subroutine rsf_write_complex_1d_n

  subroutine rsf_write_complex_1d (hist, array)    
    type (file),            intent (in)  :: hist
    complex, dimension (:), intent (in) :: array

    call sf_complexwrite(array, size(array), hist%tag)
  end subroutine rsf_write_complex_1d

  subroutine rsf_write_complex_2d (hist, array)    
    type (file),              intent (in)  :: hist
    complex, dimension (:,:), intent (in) :: array

    call rsf_write_complex_1d_n(hist,array,size(array))
  end subroutine rsf_write_complex_2d

  subroutine rsf_write_complex_3d (hist, array)    
    type (file),                intent (in)  :: hist
    complex, dimension (:,:,:), intent (in) :: array
    
    call rsf_write_complex_1d_n(hist,array,size(array))
  end subroutine rsf_write_complex_3d

  subroutine rsf_write_complex_4d (hist, array)    
    type (file),                  intent (in)  :: hist
    complex, dimension (:,:,:,:), intent (in) :: array
    
    call rsf_write_complex_1d_n(hist,array,size(array))
  end subroutine rsf_write_complex_4d

  subroutine rsf_write_complex_5d (hist, array)    
    type (file),                    intent (in)  :: hist
    complex, dimension (:,:,:,:,:), intent (in) :: array
    
    call rsf_write_complex_1d_n(hist,array,size(array))
  end subroutine rsf_write_complex_5d

  !------------------------------------------------------------

  subroutine iaxa(FF,AA,i)
    type(file), intent(in) :: FF
    integer   , intent(in) :: i
    type(axa),  intent(out):: AA
    character(len=128)     :: BB

    write(BB,"(a,i1)" ) 'n',i
    call from_par(FF,BB,AA%n,1)
    write(BB,"(a,i1)" ) 'o',i
    call from_par(FF,BB,AA%o,0.)
    write(BB,"(a,i1)" ) 'd',i
    call from_par(FF,BB,AA%d,1.)

  end subroutine iaxa

  !------------------------------------------------------------

  subroutine oaxa(FF,AA,i)
    type(file), intent(in) :: FF
    integer   , intent(in) :: i
    type(axa),  intent(in) :: AA
    character(len=128)     :: BB

    write(BB,"(a,i1)" ) 'n',i
    call to_par(FF,BB,AA%n)
    write(BB,"(a,i1)" ) 'o',i
    call to_par(FF,BB,AA%o)
    write(BB,"(a,i1)" ) 'd',i
    call to_par(FF,BB,AA%d)

  end subroutine oaxa
  !------------------------------------------------------------

  subroutine raxa(AA)
    type(axa),  intent(in) :: AA

    write(0,*) AA%n,AA%o,AA%d
  end subroutine raxa

  !------------------------------------------------------------

  subroutine rsf_fileflush(hist, hist_src)
    type(file) :: hist
    type(file), optional :: hist_src
    integer(kind=PTRKIND) :: M_NULL = 0
    if (present(hist_src)) then
      call sf_fileflush(hist%tag, hist_src%tag)
    else
      call sf_fileflush(hist%tag, M_NULL)
    end if
  end subroutine rsf_fileflush

  !------------------------------------------------------------

  subroutine rsf_seek(hist, offset, whence)
    type(file) :: hist
    integer(kind=OFFKIND) :: offset
    integer :: whence
    call sf_seek(hist%tag, offset, whence)
  end subroutine rsf_seek

  
end module RSF

