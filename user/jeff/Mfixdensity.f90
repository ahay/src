program fixdensity
  use rsf
  implicit none
  
  integer :: nx,nz,iz,ix,fix
  real, dimension(:,:), allocatable :: in

  type(file) :: infile,outfile

  call sf_init()
  infile  = rsf_input("in")
  outfile = rsf_output("out")

  call from_par(infile,"n1",nz)
  call from_par(infile,"n2",nx)
  allocate( in(nz,nx) )

  call rsf_read(infile,in)

  do ix=1,nx
     fix=0
     do iz=1,nz
        if( in(iz,ix) .gt. 1.1) then
           in(iz:iz+2,ix) = in(1:3,ix)
           fix=1
        end if
        if (fix .eq. 1) exit
     end do
  end do

  call rsf_write(outfile,in)

end program fixdensity
