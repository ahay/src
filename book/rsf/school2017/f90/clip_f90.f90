! Clip the data.
program Clip
  use rsf ! use module


  implicit none
  type (file)                      :: in, out
  integer                          :: n1, n2, i1, i2
  real                             :: upper, lower
  real, dimension (:), allocatable :: trace

  ! initialize RSF
  call sf_init() 
  ! Standard input file (void)
  in = rsf_input()
  ! Standard output file (void)
  out = rsf_output()

  ! check that the input is float
  if (sf_float /= gettype(in)) call sf_error("Need floats")

  ! n1 is the fastest dimension
  call from_par(in,"n1",n1)

  ! leftsize gets n2*n3*n4*...
  n2 = filesize(in,1)

  ! parameter form the command line
  call from_par("upper",upper,10000.) 
  call from_par("lower",lower,-10000.) 

  ! allocate floating point array
  allocate (trace (n1))


  ! loop over traces
  do i2=1, n2 
     
      ! read a trace
     call rsf_read(in,trace)
    
     ! loop over samples
     where (trace > upper) trace = upper
     where (trace < lower) trace = lower

     ! write a trace
     call rsf_write(out,trace)
  end do
end program Clip
