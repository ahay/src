! 
! time-domain acoustic FD modeling
! 

program ERFDMf90
  use rsf
  use cub

  implicit none

  logical    :: verb ! verbose flag

  ! I/O files
  type(file) :: Fw,Fv,Fr
  type(file) :: Fo

  ! cube axes
  type(axis) :: at,az,ax
  integer    :: it,iz,ix
  real       :: idx,idz,dt2

  ! arrays
  real, allocatable :: ww(:), vv(:,:), rr(:,:)
  real, allocatable :: um(:,:), uo(:,:), up(:,:), ud(:,:)

  ! Laplacian coefficients
  real :: c0=-30./12.
  real :: c1=+16./12.
  real :: c2=- 1./12.

  ! init RSF
  call sf_init()
  call from_par("verb",verb,.false.)

  ! setup I/O files
  Fw=rsf_input ("in")
  Fv=rsf_input ("vel")
  Fr=rsf_input ("ref")
  Fo=rsf_output("out")

  ! read axes
  call iaxis(Fw,at,1)
  call iaxis(Fv,az,1)
  call iaxis(Fv,ax,2)

  ! setup output header
  call oaxis(Fo,az,1)
  call oaxis(Fo,ax,2)
  call oaxis(Fo,at,3)

  dt2 =    at%d*at%d
  idz = 1/(az%d*az%d)
  idx = 1/(ax%d*ax%d) 

  ! read wavelet, velocity & reflectivity
  allocate(ww(at%n     )); ww=0.
  allocate(vv(az%n,ax%n)); vv=0.
  allocate(rr(az%n,ax%n)); rr=0.

  call rsf_read(Fw,ww)
  call rsf_read(Fv,vv)
  call rsf_read(Fr,rr)

  ! allocate temporary arrays
  allocate(um(az%n,ax%n)); um=0.
  allocate(uo(az%n,ax%n)); uo=0.
  allocate(up(az%n,ax%n)); up=0.
  allocate(ud(az%n,ax%n)); ud=0.

  ! 
  ! MAIN LOOP
  ! 
  do it=1,at%n
     if(verb) write (0,*) it

     ! 4th order laplacian
     do iz=2,az%n-2
        do ix=2,ax%n-2
           ud(iz,ix) = &
                c0* uo(iz,  ix  ) * (idx + idz)        + &
                c1*(uo(iz  ,ix-1) + uo(iz  ,ix+1))*idx + &
                c2*(uo(iz  ,ix-2) + uo(iz  ,ix+2))*idx + &
                c1*(uo(iz-1,ix  ) + uo(iz+1,ix  ))*idz + &
                c2*(uo(iz-2,ix  ) + uo(iz+2,ix  ))*idz
        end do
     end do

     ud= ud *vv*vv

     ! inject wavelet
     ud = ud - ww(it) * rr

     ! time step
     up = 2*uo - um + ud * dt2
     um =   uo
     uo =   up

     ! write wavefield to output
     call rsf_write(Fo,uo)
  end do

  call exit(0)
end program ERFDMf90
