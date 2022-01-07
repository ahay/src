! time-domain acoustic FD modeling
program AFDMf90
  use rsf

  implicit none

  ! Laplacian coefficients
  real :: c0=-30./12.,c1=+16./12.,c2=- 1./12.

  logical    :: verb         ! verbose flag
  type(file) :: Fw,Fv,Fr,Fo  ! I/O files
  type(axa)  :: at,az,ax     ! cube axes
  integer    :: it,iz,ix     ! index variables
  integer    :: nt,nz,nx
  real       :: dt,dz,dx
  real       :: idx,idz,dt2

  real, allocatable :: vv(:,:),rr(:,:),ww(:)           ! I/O arrays
  real, allocatable :: um(:,:),uo(:,:),up(:,:),ud(:,:) ! tmp arrays

  call sf_init() ! init RSF
  call from_par("verb",verb,.false.)

  ! setup I/O files
  Fw=rsf_input ("in")
  Fv=rsf_input ("vel")
  Fr=rsf_input ("ref")
  Fo=rsf_output("out")

  ! Read/Write axes
  call iaxa(Fw,at,1); nt = at%n; dt = at%d
  call iaxa(Fv,az,1); nz = az%n; dz = az%d
  call iaxa(Fv,ax,2); nx = ax%n; dx = ax%d
  
  call oaxa(Fo,az,1)
  call oaxa(Fo,ax,2)
  call oaxa(Fo,at,3)

  dt2 =    dt*dt
  idz = 1/(dz*dz)
  idx = 1/(dx*dx) 

  ! read wavelet, velocity & reflectivity
  allocate(ww(nt));    call rsf_read(Fw,ww)
  allocate(vv(nz,nx)); call rsf_read(Fv,vv)
  allocate(rr(nz,nx)); call rsf_read(Fr,rr)

  ! allocate temporary arrays
  allocate(um(nz,nx)); um=0.
  allocate(uo(nz,nx)); uo=0.
  allocate(up(nz,nx)); up=0.
  allocate(ud(nz,nx)); ud=0.

  ! MAIN LOOP
  do it=1,nt
     if(verb) write (0,*) it

     ud(3:nz-2,3:nx-2) = &
          c0* uo(3:nz-2,3:nx-2) * (idx + idz)            + &
          c1*(uo(3:nz-2,2:nx-3) + uo(3:nz-2,4:nx-1))*idx + &
          c2*(uo(3:nz-2,1:nx-4) + uo(3:nz-2,5:nx  ))*idx + &
          c1*(uo(2:nz-3,3:nx-2) + uo(4:nz-1,3:nx-2))*idz + &
          c2*(uo(1:nz-4,3:nx-2) + uo(5:nz  ,3:nx-2))*idz

     ! inject wavelet
     ud = ud - ww(it) * rr

     ! scale by velocity
     ud= ud *vv*vv

     ! time step
     up = 2*uo - um + ud * dt2
     um =   uo
     uo =   up

     ! write wavefield to output
     call rsf_write(Fo,uo)
  end do

  call exit(0)
end program AFDMf90
