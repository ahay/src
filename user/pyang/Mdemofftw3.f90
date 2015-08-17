! The demofftw3.f90 is an example to demonstrate the use of
! FFTW3 in legacy Fortran (Fortran 90).
!    1. DO NOT use implicit none when calling fftw3 
!    2. changes have to be made:
!	the starting charaters: 'fftw' to be 'dfftw'
!	fftw_execute(plan)-->dfftw_execute_dft(plan,in,out)
!
! More details of the use of FFTW3 can be found in FFW3 mannual:
!	chapter 8: Calling FFTW from Legacy Fortran		

program Mdemofftw3
	!include 'fftw3.h'
	implicit none

	integer, parameter::N=32
  	integer(kind=4), parameter :: fftw_estimate = 64 ! do NOT use fftw_measure=0;
	integer(kind=4), parameter:: fftw_forward=-1
	
	double complex, dimension(N):: in,out
	integer*8 plan
	integer :: i


	call generate_array(in,N)

	call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_estimate)
	call dfftw_execute_dft(plan, in, out)
	call dfftw_destroy_plan(plan)

	do i = 1,N
        	write(*,*) i, in(i), out(i)
    	end do
end program Mdemofftw3

subroutine generate_array(input,N)
	implicit none
	integer::N
	double complex::input(N)
	 
	integer::i
	!real::mydata(2,2)
	!mydata=[3,4,5] !equivalent to: mydata=(/3,4,5/)
	!data mydata /4*1.0/

	input=(/(i*(1.,0.),i=1,N)/)
	!data in /N*(1.,0.)/
end subroutine generate_array
