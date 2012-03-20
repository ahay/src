! Program to Rotate from field [Ux,Uy,Uz] coordinates to [P,Sv,Sh] (FSRM - Kennett, 198x)
program rec2ps 
  use rsf
  use rotate !! Rotation module

  implicit none 

  integer                            :: n1,n2,n3		
  real, dimension(:,:,:), allocatable:: data,modl
  real	                             :: alpha,beta,p1,p2,irot
  logical                            :: adj
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
  call from_par("adj", adj,.true.)
  allocate( data(n1,n2,n3),modl(n1,n2,n3) )

  call rotatedata_init(n1,n2,n3,alpha,beta,p1,p2,irot)

!! . . Read in File
  if (adj) then
     call rsf_read(infile,data)
  else
     call rsf_read(infile,modl)
  end if

  call rotatedata(adj,.false.,data,modl)

  !! . . Output file
  if (adj) then
     call rsf_write(outfile,modl)
  else
     call rsf_write(outfile,data)
  end if

  deallocate(data,modl)
  call exit()
end program rec2ps 
