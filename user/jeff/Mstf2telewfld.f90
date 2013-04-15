program stf2telewfld
  use rsf
  implicit none
  
  ! Projects a source time function into a 2-D space, along a planar wavefront (as in a teleseismic arrival) for a given ray parameter and background velocity.
  ! This is adapted from Jeff's original version, which used a single Gaussian pulse (and which is now appended -singleGaussian.f90).
  ! This version uses an arbitrary time series as input.
  
  integer :: nx,nz,iz,ix,fix,nt
  real    :: dx,dz,ox,oz,p1,z0,delz0,v0,x,z,dt,dts,L,T,m,b
  real, dimension(:),     allocatable :: stf
  real, dimension(:,:),   allocatable :: in
  real, dimension(:,:,:), allocatable :: uo,um

  type(file) :: infile,stffile,outfile,oldfile

  call sf_init()
  infile = rsf_input("in") !! Input file with correct z,x dimensions
  stffile = rsf_input("stf")  !! Source Time Function
  outfile = rsf_output("out") !! Current wavefield
  oldfile = rsf_output("old") !! Previous wavefield

  !! . . Input dimensions
  call from_par(infile,"n1",nz)
  call from_par(infile,"n2",nx)
  call from_par(infile,"o1",oz)
  call from_par(infile,"o2",ox)
  call from_par(infile,"d1",dz)
  call from_par(infile,"d2",dx)
  call from_par(stffile,"n1",nt)  ! n of samples in the source
  call from_par(stffile,"d1",dts) ! dt of the source
 
  call from_par("p1",p1,0.) !! Inline planewave dip (s/km)
  call from_par("z0",z0,40.) !! Depth
  call from_par("v0",v0,8.0)
  call from_par("dt",dt,0.001)

  !! . . Output files
  !! . . current wfld
  call to_par(outfile,"n1",nz)
  call to_par(outfile,"n2",nx)
  call to_par(outfile,"n3",2)
  call to_par(outfile,"o1",oz)
  call to_par(outfile,"o2",ox)
  call to_par(outfile,"o3",0.)
  call to_par(outfile,"d1",dz)
  call to_par(outfile,"d2",dx)
  call to_par(outfile,"d3",1.)
  !! . . Previous wfld
  call to_par(oldfile,"n1",nz)
  call to_par(oldfile,"n2",nx)
  call to_par(oldfile,"n3",2)
  call to_par(oldfile,"o1",oz)
  call to_par(oldfile,"o2",ox)
  call to_par(oldfile,"o3",0.)
  call to_par(oldfile,"d1",dz)
  call to_par(oldfile,"d2",dx)
  call to_par(oldfile,"d3",1.)

  !! . . Allocation
  allocate( in(nz,nx) )
  allocate( stf(nt) )
  allocate( uo(nz,nx,2) )
  allocate( um(nz,nx,2) )

  !! . . Read in Source Time Function
  call rsf_read(stffile,stf)
  
  !! . . Set up current
  
  !! . . Equation
  !! . . z = z0+(x-ox)*p1*v0
  !! . . Aexp(-((x-ox)*(x-ox)+((x-ox)*p1*v0)*((x-ox)*p1*v0))/std)
  !! . . amplitude
  !! . . . . UZ = cos(p1*v0)
  !! . . . . UX = sin(p1*v0)
  
  write(0,*) 'Creating planewave with dip: ',asin(p1*v0)*180./3.1415926
  
  ! . . Calculate some stuff ahead of time
  delz0  = v0*dt*cos(asin(p1*v0)) ! magnitude of shift in Z-dimension by one time step for wavefield at previous timestep
  m = tan(asin(p1*v0)) ! slope of wavefront
  
  if (p1 .gt. 0.) then   ! Ray parameter is positive
    b = -z0 - (m*nx*dx)  ! Y-intercept, for distance calculation
  else                   ! Ray parameter is negative
    b = -z0
  end if

  do iz=1,nz
    z = - (oz+(iz-1)*dz) ! Make negative so we can use a sane coordinate system
      do ix=1,nx
        x = ox+(ix-1)*dx
        if ( abs(z) .ge. abs(m*x+b) ) then ! we've reached the wavefront
          L=abs(z-m*x-b)/sqrt(m*m+1.) ! distance b/wn current point and wavefront
          T=L/v0 ! Time along the wavetrain to current point
!          if (nint(T/dts).gt.nt) then
!            write(0,*) 'Warning: Not enough samples in the source time series!!!'
!            write(0,*) 'L, T, nt:', L, T, nt
!          end if
          
 	      uo(iz,ix,1) = cos(p1*v0) * stf(nint(T/dts))
          uo(iz,ix,2) = sin(p1*v0) * stf(nint(T/dts))
        else
          uo(iz,ix,1) = 0.
          uo(iz,ix,2) = 0.
	      end if
	      if ( abs(z) .ge. abs(m*x+b)+delz0 ) then ! we've reached the previous step's wavefront
	        L=abs( z-m*x-(b-delz0) )/sqrt(m*m+1.) ! distance b/wn current point and wavefront at previous timestep
	        T=L/v0 ! Time along the wavetrain to current point
	        um(iz,ix,1) = cos(p1*v0) * stf(nint(T/dts))
	        um(iz,ix,2) = sin(p1*v0) * stf(nint(T/dts))
        else
          uo(iz,ix,1) = 0.
          uo(iz,ix,2) = 0.
	      end if
	   end do
  end do

  !! . . Current wavefield
  call rsf_write(outfile,uo)

  !! . . Previous wavefield
  call rsf_write(oldfile,um)

  !! . . Deallocation
  deallocate( in, stf, uo, um )

end program stf2telewfld



