program Clipit
  use rsf

  implicit none
  type (file)                      :: in, out
  integer                          :: n1, n2, i1, i2
  real                             :: clip
  real, dimension (:), allocatable :: trace

  call sf_init()
  in = rsf_input()
  out = rsf_output()

  call from_par(in,"n1",n1)
  n2 = filesize(in,1)

  call from_par("clip",clip)

  allocate (trace (n1))

  do i2=1, n2
     call rsf_read(in,trace)
     
     where (trace >  clip) trace =  clip
     where (trace < -clip) trace = -clip

     call rsf_write(out,trace)
  end do
end program Clipit
