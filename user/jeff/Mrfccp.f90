program rfccp
  use rsf
  implicit none
  
  ! Back-projection of receiver functions in 1D velocity model for CCP stacking.
  ! Adapted from I.W. Bailey's Matlab implementation of ray path calcutation.
  ! Warren Caldwell, August 2012
  
  integer             :: n1,n2,nz
  real                :: d1,d2,dz,dt,o1,o2,oz,toff,p,dx,x0,z0,dSx2,tcum
  integer             :: i,j,k,jcurr,roff,tcoff,tcnt1,tcnt2,cnt,verb
  real,allocatable    :: vp(:),vs(:),v(:),dSx(:),q(:)
  real,allocatable    :: RFs(:,:),CCP(:,:)
  integer, parameter  :: maxcnt = 2500 ! Infinite loop killer

  type(file) :: infile,vpfile,vsfile,outfile

  call sf_init()
  infile  = rsf_input( "in" )        !! Input file with receiver functions
  vpfile  = rsf_input( "vp" )        !! Vp (km/s) at each layer
  vsfile  = rsf_input( "vs" )        !! Vs (km/s) at each layer
  outfile = rsf_output("out")

  ! Read and write History
  call from_par(infile,'n1',n1) ! 1-axis is time
  call from_par(infile,'n2',n2) ! 2-axis is traces
!  call from_par(infile,'n3',n3) ! 3-axis is events
  call from_par(infile,'d1',d1)
  call from_par(infile,'d2',d2)
  call from_par(infile,'o1',o1)
  call from_par(infile,'o2',o2)
  call from_par(infile,'p1',p)  ! ray parameter, horizontal slowness (s/km)
  call from_par(vpfile,'n1',nz)
  call from_par(vpfile,'d1',dz)
  call from_par(vpfile,'o1',oz)
  call from_par('t_off',toff,10.) ! time offset of receiver functions
  call from_par('verb',verb,0)
  ! 'verb' is verbosity level:
  ! 0 = print nothing
  ! 1 = print information pertaining to each run of the program
  ! 2 = print information after each run of the subroutine
  ! 3 = print information from within each run of the subroutine

  call to_par(outfile,"n1",nz)
  call to_par(outfile,"n2",n2)
  call to_par(outfile,"d1",dz)
  call to_par(outfile,"d2",d2)
  call to_par(outfile,"o1",oz)
  call to_par(outfile,"o2",o2)
  
  ! Allocate memory
  allocate ( vp(nz),vs(nz),v(nz),dSx(nz),q(nz) )
  allocate ( RFs(n1,n2),CCP(nz,n2) )
  CCP=0.
  write(0,*) 'CCP is MxN:', size(CCP,1), size(CCP,2)

  ! Read the data
  call rsf_read(infile,RFs)
  call rsf_read(vpfile,vp)
  call rsf_read(vsfile,vs)

  dx=d2                         ! For clarity
  dt=d1                         ! For clarity
  tcoff=toff/dt                 ! Offset of t=0 of RF traces (in samples)
  v=(vp*vs)/(vp-vs)             ! V_p-s, the P-S offset velocity for each layer (km/s)
  if ( verb.ge.1 ) write(0,*) 'v: ', v
  q = sqrt((1/v)**2 - p**2)     ! Vertical slowness in each layer in S/km
  if ( verb.ge.1 ) write(0,*) 'q: ', q
  dSx = p * dz / q              ! Horizontal travel steps in each layer in km
  if ( verb.ge.1 ) write(0,*) 'dSx: ', dSx
  
  do j=1,n2                     ! Loop over traces
    jcurr=j                     ! Start in current column
    cnt=0                       ! Initialize infinite loop killer
 !   write(0,*) '******** TRACE NUMBER: ', j
    dSx2=dSx(1)                 ! Remaining horiz dist in row; seed with first value
    x0=dx/2.                    ! Starting coordinates of trace
    z0=0.                       ! Starting coordinates of trace
    tcnt1=1                     ! First time sample is first unused time sample
    tcum=0.                     ! Cumulative time count for trace
    i=1                         ! Start in topmost layer
    
 !   if ( verb.ge.0 ) write(0,*) 'pre-sub: ',x0,z0,dSx2
    
    ! Loop over layers and/or time samples:
    do while ( i.le.nz .and. tcnt1+tcoff.le.n1 .and. jcurr.le.n2 .and. jcurr.ge.1) !JCS was jcurr.ge.1
      cnt=cnt+1
      if ( cnt > maxcnt ) STOP  ! Kill infinite loop
 !     if ( verb.ge.0 ) write(0,*) cnt,j,i
      
      call trace_block(dx,dz,dt,1/v(i),p,q(i),x0,z0,dSx2,tcum,tcnt1,tcnt2,roff,verb)
      
!      if ( verb.ge.0 ) write(0,*) 'post-sub:',x0,z0,dSx2,tcum,tcnt1,tcnt2,roff
      
      if ( tcnt1 .ne. tcnt2 ) then
        do k=tcnt1,tcnt2-1        ! Loop through relevant time samples
!          if ( verb.ge.0 ) write(0,*) 'samples',i,jcurr,j
          CCP(i,jcurr) = CCP(i,jcurr) + RFs(k+tcoff,j) ! Add to CCP stack
        end do
        tcnt1=tcnt2             ! Update highest unused time sample
      end if
  	  
  	  if ( roff .eq. 0. ) then  ! We reached the bottom of the row
  	    i=i+1                   ! Move to next layer of the model
  	    if ( i.le.nz ) dSx2=dSx(i) ! Set the new X distance for the next layer
  	  else                      ! Ray moves into an adjacent column, specified by 'off'
  	    jcurr=jcurr+roff        ! Shift current column for next iteration
  	  end if

    end do
  end do
  
  call rsf_write(outfile,CCP)
!  deallocate(RFs,CCP)
!  deallocate(vp,vs,v,dSx,q)
  
end program rfccp
  
   ! Traces a ray through a block of the model
  subroutine trace_block(dx,dz,dt,p,px,pz,x0,z0,dSx,tcum,tcnt1,tcnt2,off,verb)
    real,    intent(in)    :: dx,dz,dt ! Dimensions of block and time sampling
    real,    intent(in)    :: p,px,pz  ! Slownesses
    integer, intent(in)    :: verb     ! Verbosity level
    real,    intent(inout) :: x0,z0    ! Where trace enters the block
    real,    intent(inout) :: dSx      ! (Remaining) Distance (to be) traveled
    real,    intent(inout) :: tcum     ! Cumulative time for current trace
    integer, intent(in)    :: tcnt1    ! Index of highest previously unused time sample
    integer, intent(out)   :: tcnt2    ! Index of highest sample falling beyond current cell
    integer, intent(out)   :: off      ! Row offset: -1,1,0 for left, right, down
    real                   :: t        ! Time in current block
    integer                :: i
    
    if ( x0+dSx .lt. dx .and. x0+dSx .gt. 0. ) then ! Case 1: ray hits bottom of cell
      off=0                       ! No row offset; move to next row on next iteration
      t=(dz-z0)*pz                ! Time in current cell
!      dS=sqrt(dSx(1)**2+dz**2)   ! Distance traveled in current cell
!      tc=floor(t/dt)              ! Count of time samples in current block
      x0=x0+dSx                   ! New x-offset
      z0=0.                       ! Reset z-offset (ray will start at top of next layer)
      if ( verb.ge.3 ) write(0,'(a,F5.4,"=(",F4.3,"-",F4.3,")*",F4.3)') 'Ray hit bottom. t=(dz-z0)*pz: ',t,dz,z0,pz

    elseif (x0+dSx .gt. dx) then  ! Case 2: ray hits right side of cell
      off=1                       ! Next iteration will move one column to the right
      t=(dx-x0)*px                ! Time in current cell
      if ( verb.ge.3 ) write(0,'(a,F5.4,"=(",F4.3,"-",F4.3,")*",F4.3)') 'Ray hit R side. t=(dx-x0)*px: ',t,dx,x0,px
!      dS=sqrt((dx-x0)**2+dSz**2) ! Distance traveled in current cell
!      tc=floor(t/dt)              ! Count of time samples in current block
      dSx=dSx-(dx-x0)             ! Adjustment to x-distance left to travel in layer
      x0=0.                       ! New x-offset (ray will start at L. of adjacent block)
      z0=z0+(t/pz)                ! New z-offset
      if (z0.gt.dz) write(0,*) "Something f'd up!" ! Sanity check

    elseif (x0+dSx .lt. 0.) then  ! Case 3: ray hits left side of cell
      off=-1                      ! Next iteration will move one column to the left
      t=x0*px                     ! Time in current cell
     if ( verb.ge.3 ) write(0,'(a,F5.4,"=",F4.3,"*",F4.3)') 'Ray hit L side. t=x0*px: ',t,x0,px
!      dS=sqrt(x0**2+dSz**2)      ! Distance traveled in current cell
!      tc=floor(t/dt)              ! Count of time samples in current block
      dSx=dSx+x0                  ! Adjustment to x-distance left to travel in layer
      x0=dx                       ! New x-offset (ray will start at R. of adjacent block)
      z0=z0+(t/pz)                ! New z-offset
      if (z0.gt.dz) write(0,*) "Something f'd up!" ! Sanity check
    end if

    i=tcnt1                       ! Start index for testing time samples
    tcnt2=tcnt1                   ! Default to zero time samples in current cell
    do while ( i*dt .le. tcum+t )
      if ( i*dt .gt. tcum ) then
        tcnt2=tcnt1+1
      end if
      i=i+1
    end do
    tcum=tcum+t                   ! Update cumulative time
    
  end subroutine trace_block
  
