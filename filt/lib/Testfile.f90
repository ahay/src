program Test
  use rsf

  implicit none
  type (file)                      :: in, out
  integer                          :: n1, i
!  character (len=100)              :: label1
  real, dimension (:), allocatable :: trace

  call sf_init()
  in = rsf_input()
  out = rsf_output("out")

  call from_par(in,"n1",n1)
  allocate (trace (n1))

  call to_par(out,"n2",5)

  call rsf_read(in,trace)
  do i=1,5
     call rsf_write(out,trace)
  end do
end program Test

!	$Id: Testfile.f90,v 1.3 2003/09/29 14:34:55 fomels Exp $	
