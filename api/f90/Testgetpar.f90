!      a=100 Xa=5
!      float=5.625 cc=fgsg
!      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"

program Test
  use rsf
  
  implicit none
  integer                :: i
  real                   :: f
  real,    dimension (4) :: d
  character (len=100)    :: str
  logical                :: yes
  logical, dimension (2) :: no

  call sf_init()
  call from_par("a",i)
  call assert(i==100)
  call from_par("c",i,0)
  call assert(i==0)
  call from_par("float",f)
  call assert(f==5.625)
  call from_par("dd",d)
  call assert(all(d==(/ 1.,4.,4.,2.25/)))
  call from_par("true",yes)
  call assert(yes)
  call from_par("false",no)
  call assert(.not. any (no))
!  str = sf_getstring("label")      

contains
  subroutine assert(cond)
    logical, intent (in) :: cond
    if (.not. cond) call sf_error("failed")
  end subroutine assert
end program Test

!	$Id: Testgetpar.f90 982 2005-01-30 23:38:22Z shan $	
