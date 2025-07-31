program rfspecdiv
  use rsf
  implicit none

  ! Calculates receiver functions via spectral division a la Langston 1979.
  ! Warren and Sjoerd, 3/23/09
  ! Modified by Warren for Madagascar, May 2012
  ! Langston uses vertical and radial; we use P and SV
  ! P component is standard input. Radial (or tangential) is 'file1'
  ! must also specify 'a=' and 'c=' (Gaussian width and water level)
  ! Jesse's example in /Volumes/fuzz/SRC/math.subs/decon.f


  integer             :: n1, n2
  real                :: d1, d2, o1, o2, w, pi
  real                :: a, c, phi, pd ! See Langston 1979: 'a' = Gaussian width, 'c' = water level.
  real,allocatable    :: maxphi(:)
  complex             :: gaus, wdel
  integer             :: i , j , k
  complex,allocatable :: dataV(:,:), dataF(:,:), dataOut(:,:)

  type(file) :: infile,file1,outfile

  call sf_init()
  infile = rsf_input("in") !! Input file with correct z,x dimensions
  file1 = rsf_input("file1") !! Tangential components
  outfile = rsf_output("out")


  ! Read and write History
  call from_par(infile,'n1',n1 ) ! 1-axis is freq
  call from_par(infile,'n2',n2 ) ! 2-axis is station
  call from_par(infile,'d1',d1)
  call from_par(infile,'d2',d2)
  call from_par(infile,'o1',o1)
  call from_par(infile,'o2',o2)

  call from_par("a",a,1.25)
  call from_par("c",c,0.01)
  call from_par("pd",pd,10.) ! phase delay
  
  call to_par  ( outfile,"n1", n1 )
  call to_par  ( outfile,"n2", n2 )
  call to_par  ( outfile,"d1", d1 )
  call to_par  ( outfile,"d2", d2 )
  call to_par  ( outfile,"o1", o1 )
  call to_par  ( outfile,"o2", o2 )

  write(0,*) '------- RFspecdiv -------'
  write(0,*) 'n1, n2:', n1,n2
  write(0,*) 'o1, o2:', o1,o2
  write(0,*) 'a:', a

  ! Allocate memory
  allocate ( dataOut(n1,n2), dataF(n1,n2), dataV(n1,n2) )
  allocate ( maxphi(n2) )

  ! Read the data
!  write(0,*) '-------------------------'
!  write(0,*) 'Read Data'
  call rsf_read(infile,dataV)
  call rsf_read(file1,dataF)

  ! Computations
!  write(0,*) '-------------------------'
!  write(0,*) 'Computations'
  
    ! determine water level:
  do j=1,n2 ! loop over traces
     maxphi(j) = 0.
     do i=1,n1 ! loop over samples
        phi = real(dataV(i,j)*conjg(dataV(i,j)))
!        write(0,*) "i,j,phi: ",i,j,phi
        if ( phi > maxphi(j) ) maxphi(j) = phi
        phi = 0.
     end do
  end do

  pi = atan(1.)*4.

  do j=1,n2 ! loop over traces
     do i=1,n1 ! loop over samples (in frequency domain)

        phi = real(dataV(i,j)*conjg(dataV(i,j)))
        !write(0,*) phi
        if ( phi < c*maxphi(j) ) phi = c*maxphi(j)
		if ( phi < c*maxphi(j) ) write(0,*) "water level correction:", i,j
        w  = 2. * pi * (float(i-1)*d1)
        gaus = cmplx( exp( -(w)**2 / (4*a**2) ), 0. )
        wdel = exp( cmplx( 0., -1*pd*w) ) ! -10iw makes +10sec phase delay
        dataOut(i,j) = ( dataF(i,j)*conjg(dataV(i,j)) / phi ) * gaus * wdel

     end do
  end do

  ! Write the data
!  write(0,*) '-------------------------'
!  write(0,*) 'Write Data'
  call rsf_write(outfile,dataOut)

  ! Deallocate memory
  deallocate(dataOut,dataF,dataV)
!  write(0,*) '-------------------------'

end program rfspecdiv
