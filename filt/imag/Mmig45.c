program FinDif
  use ctridiagonal
  use sep

  implicit none
  logical                               :: inv
  integer                               :: nw,nz,nx, iw,ix,iz
  real                                  :: dw,dz,dx, vel0,pi, eps, den, beta
  complex                               :: w, a
  real,    dimension (:,:), allocatable :: depth, vel, voff
  real,    dimension (:),   allocatable :: time
  complex, dimension (:),   allocatable :: ctime, tt, diag1,diag2, offd1,offd2
  type (ctris)                          :: slv

  call sep_init ()
  call from_par("inv",inv,.false.)
  call from_par("eps",eps,0.01)
  call from_par("beta",beta,1./12.)
  if (inv) then ! modeling
     call from_history (nz,nx)
     call from_history ("d1",dz)
     call from_history ("d2",dx) 
     call from_par ("nt",nw)
     call from_par ("dt",dw)
     dw = 1./(dw*(nw-1))
     call to_history ("n2",nw); call to_history ("d2",dw)
     call to_history ("n1",nx); call to_history ("d1",dx)
  else          ! migration
     call from_history (nx,nw)
     call from_history ("d2",dw) 
     call from_history ("d1",dx) 
     if (exist_file ("velocity")) then
        call from_aux ("velocity", "n1", nz)
        call from_aux ("velocity", "d1", dz)
     else
        call from_par ("nz",nz) 
        call from_par ("dz",dz)
     end if
     call to_history ("n1",nz); call to_history ("d1",dz)
     call to_history ("n2",nx); call to_history ("d2",dx)
  end if

  allocate (vel (nz,nx), voff(nz,nx))
  if (exist_file ("velocity")) then
     call sep_read (vel, "velocity")
  else ! constant velocity
     call from_par ("vel", vel0)
     vel = vel0
  end if
  call sep_close ()
  pi = acos (-1.)
  dw = dw*pi
  vel = 0.5*vel
  voff (:,:nx-1) = sqrt(vel(:,:nx-1)*vel(:,2:)) 
  dx = 0.25/(dx*dx)

  allocate (depth (nz,nx))
  allocate (time (nx), ctime (nx), tt (nx))
  allocate (diag1 (nx), offd1 (nx), diag2 (nx), offd2 (nx))
  if (inv) then
     call sep_read (depth)
  else
     depth = 0.
  end if
  call ctridiagonal_init (nx, slv)

!           (1. + k2*0.25*(1.-w*dz)/(w*w))/ &
!           (1. + k2*0.25*(1.+w*dz)/(w*w))

  do iw = 1, nw
     write (0,*) iw, nw
     if (inv) then ! modeling
        w = cmplx(eps*dw,(iw-1)*dw)
        ctime (:) = depth (nz,:)
        do iz = nz-1, 1, -1
           diag1 =   -2.*(beta - (vel(iz,:)/w-dz)*vel(iz,:)*dx/w)
           diag2 = 1.-2.*(beta - (vel(iz,:)/w+dz)*vel(iz,:)*dx/w)
           
           a = cexp(-0.5*w/(vel(iz,1)*sqrt(dx)))

           diag1(1) =    (a-2.)*(beta - (vel(iz,1)/w-dz)*vel(iz,1)*dx/w)
           diag2(1) = 1.+(a-2.)*(beta - (vel(iz,1)/w+dz)*vel(iz,1)*dx/w)

           a = cexp(-0.5*w/(vel(iz,nx)*sqrt(dx)))

           diag1(nx) =    (a-2.)*(beta - (vel(iz,nx)/w-dz)*vel(iz,nx)*dx/w)
           diag2(nx) = 1.+(a-2.)*(beta - (vel(iz,nx)/w+dz)*vel(iz,nx)*dx/w)

           offd1 = beta - (voff(iz,:)/w-dz)*voff(iz,:)*dx/w
           offd2 = beta - (voff(iz,:)/w+dz)*voff(iz,:)*dx/w

           tt(1) = diag1(1)*ctime(1) + offd1(1)*ctime(2)
           tt(2:nx-1) = &
                offd1(1:nx-2)*ctime(1:nx-2) + &
                diag1(2:nx-1)*ctime(2:nx-1) + &
                offd1(2:nx-1)*ctime(3:nx)
           tt(nx) = offd1(nx-1)*ctime(nx-1) + diag1(nx)*ctime(nx)
 
           ctime = ctime + tt

           call ctridiagonal_define (slv, diag2, offd2)
           call ctridiagonal_solve (slv, ctime)

           ctime = ctime * cexp(-w*dz/vel(iz,:)) + depth (iz,:)     
        end do
        time = real (ctime)
        call sep_write (time)
     else ! migration
        w = cmplx(eps*dw,-(iw-1)*dw)
        call sep_read (time)
        ctime = time
        do iz = 1, nz
           depth (iz,:) = depth(iz,:) + real (ctime)

           diag1 =   -2.*(beta - (vel(iz,:)/w-dz)*vel(iz,:)*dx/w)
           diag2 = 1.-2.*(beta - (vel(iz,:)/w+dz)*vel(iz,:)*dx/w)

           a = cexp(-0.5*w/(vel(iz,1)*sqrt(dx)))

           diag1(1) =    (a-2.)*(beta - (vel(iz,1)/w-dz)*vel(iz,1)*dx/w)
           diag2(1) = 1.+(a-2.)*(beta - (vel(iz,1)/w+dz)*vel(iz,1)*dx/w)

           a = cexp(-0.5*w/(vel(iz,nx)*sqrt(dx)))

           diag1(nx) =    (a-2.)*(beta - (vel(iz,nx)/w-dz)*vel(iz,nx)*dx/w)
           diag2(nx) = 1.+(a-2.)*(beta - (vel(iz,nx)/w+dz)*vel(iz,nx)*dx/w)

           offd1 = beta - (voff(iz,:)/w-dz)*voff(iz,:)*dx/w
           offd2 = beta - (voff(iz,:)/w+dz)*voff(iz,:)*dx/w

           tt(1) = diag1(1)*ctime(1) + offd1(1)*ctime(2)
           tt(2:nx-1) = &
                offd1(1:nx-2)*ctime(1:nx-2) + &
                diag1(2:nx-1)*ctime(2:nx-1) + &
                offd1(2:nx-1)*ctime(3:nx)
           tt(nx) = offd1(nx-1)*ctime(nx-1) + diag1(nx)*ctime(nx)
  
           ctime = ctime + tt

           call ctridiagonal_define (slv, diag2, offd2)
           call ctridiagonal_solve (slv, ctime)

           ctime = ctime * cexp(-w*dz/vel(iz,:)) 
        end do
     end if
  end do
  call ctridiagonal_close (slv)
  if (.not. inv) call sep_write (depth)

  deallocate (depth, time, ctime, vel, voff, tt)
  deallocate (diag1, offd1, diag2, offd2)
  call exit (0)
end program FinDif


