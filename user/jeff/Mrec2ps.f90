! Program to Rotate from field [Ux,Uy,Uz] coordinates to [P,Sv,Sh] (FSRM - Kennett, 198x)
program rec2ps 
  use rsf
  use rotate !! Rotation module

  implicit none 

  integer                                 :: n1,n2,n3		
  real, dimension(:,:,:), allocatable     :: in,out
  real	                                  :: alpha,beta,p1,p2,irot

  type(file) :: infile,outfile
  
  call sf_init()
  infile = rsf_input ("in")
  outfile= rsf_output("out")

  !! Parameters from header file
  call from_par(infile,"n1",n1) ! Number of time samples
  call from_par(infile,"n2",n2) ! Number of components
  call from_par(infile,"n3",n3) ! Number of traces
  call from_par(infile,"p1",p1) ! Ray parameter of planewave in 1 direction
  call from_par(infile,"p2",p2) ! Ray parameter of planewave in 2 direction

  !! Parameters required from command line (with default values)
  call from_par("alpha", alpha,6.2) ! P-wave Velocity at surface
  call from_par("beta",  beta,3.5) ! S-wave velocity at surface
  call from_par("irot",  irot, 0.) ! Rotation of array w.r.t. 1st axis

  allocate( in(n1,n2,n3),out(n1,n2,n3) )

!! . . Read in File
  call rsf_read(infile,in)
  call rotatedata_init(n1,n2,n3,alpha,beta,p1,p2,irot)
  call rotatedata(.true.,.false.,in,out)

  !! . . Output file
  call rsf_write(outfile,out)
  deallocate(in,out)
  call exit()
end program rec2ps 
