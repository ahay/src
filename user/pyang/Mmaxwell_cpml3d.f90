!  Backward reconstruction based upon Maxwell attenuation model and CPML ABC
!  It is shown that even with attenuation, reverse reconstruction is still a 
!  feasible way to build the incident wavefield, using boundary saving scheme.
!  Allowing for the heavy burden of the boundary saving in 3D, we prefer 4-th 
!  order FD in space.
!
!!$  Copyright (C) 2015 University Joseph Fourier, Grenoble (Pengliang Yang)
!!$  
!!$  This program is free software; you can redistribute it and/or modify
!!$  it under the terms of the GNU General Public License as published by
!!$  the Free Software Foundation; either version 2 of the License, or
!!$  (at your option) any later version.
!!$  
!!$  This program is distributed in the hope that it will be useful,
!!$  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!$  GNU General Public License for more details.
!!$  
!!$  You should have received a copy of the GNU General Public License
!!$  along with this program; if not, write to the Free Software
!!$  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
program mexwell_cpml2_backward
  use rsf
  implicit none

  logical :: order1,attenuating
  integer :: ib, it, kt, nt, nz, nx, ny, nb, sz, sx, sy, nzpad, nxpad, nypad
  real  :: dt, dz, dx, dy, fm, tmp, idz, idx, idy
  real*8, parameter::PI=4.*atan(1.)
  real, dimension(:), allocatable :: wlt,bndr
  real, dimension(:,:,:), allocatable :: v0, vv, rho, eta
  real, dimension(:,:,:), allocatable :: p, vz, vx, vy
  real, dimension(:,:,:,:),allocatable :: conv_pz,conv_px,conv_py
  real, dimension(:,:,:,:),allocatable :: conv_vz,conv_vx,conv_vy
  real, dimension(:,:,:,:,:),allocatable:: bvz, bvx, bvy
  type(file) :: Fv, Fw1, Fw2, Frho, Feta  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fv = rsf_input("in")    ! source position 
  Frho=rsf_input("rho")   ! density
  Feta=rsf_input("eta")   ! Pascal
  Fw1 =rsf_output("out")  ! output forward wavefield 
  Fw2 =rsf_output("back") ! output backward reconstructed wavefield

  ! Read/Write axes
  call from_par(Fv,"n1",nz) ! velocity model: nz
  call from_par(Fv,"n2",nx) ! velocity model: nx
  call from_par(Fv,"n3",ny) ! velocity model: ny
  call from_par(Fv,"d1",dz) ! velocity model: dz
  call from_par(Fv,"d2",dx) ! velocity model: dx
  call from_par(Fv,"d3",dy) ! velocity model: dy

  call from_par("nb", nb, 30) ! thinkness of sponge ABC
  call from_par("nt", nt, 1000) !number of time steps
  call from_par("dt", dt, 0.001) ! time sampling interval
  call from_par("fm", fm, 20.) ! domainant frequency for ricker wavelet
  call from_par("order1",order1,.true.) ! 1st order or 2nd order accuracy
  call from_par("attenuating",attenuating,.true.) ! add attenuation or not
  call from_par("kt", kt, 500) ! recording time of snapshot

  idz=1./dz
  idx=1./dx
  idy=1./dy
  nzpad=nz+2*nb
  nxpad=nx+2*nb
  nypad=ny+2*nb
  sz=nzpad/2
  sx=nxpad/2
  sy=nypad/2

  allocate(wlt(nt))
  allocate(bndr(nb))
  allocate(v0(nz,nx,ny))
  allocate(vv(nzpad,nxpad,nypad))
  allocate(rho(nzpad,nxpad,nypad))
  allocate(eta(nzpad,nxpad,nypad))
  allocate(p(nzpad,nxpad,nypad))
  allocate(vz(nzpad,nxpad,nypad))
  allocate(vx(nzpad,nxpad,nypad))
  allocate(vy(nzpad,nxpad,nypad))
  allocate(conv_pz(nb,nxpad,nypad,2))
  allocate(conv_px(nzpad,nb,nypad,2))
  allocate(conv_py(nzpad,nxpad,nb,2))
  allocate(conv_vz(nb,nxpad,nypad,2))
  allocate(conv_vx(nzpad,nb,nypad,2))
  allocate(conv_vy(nzpad,nxpad,nb,2))
  allocate(bvz(3,nx,ny,2,nt))
  allocate(bvx(nz,3,ny,2,nt))
  allocate(bvy(nz,nx,3,2,nt))

  !generate ricker wavelet with a delay
  do it=1,nt  
     tmp=PI*fm*(it*dt-1.0/fm)
     tmp=tmp*tmp
     wlt(it)=(1.0-2.0*tmp)*exp(-tmp)
  enddo
  !generate coefficients for the absorbing boundary
  call cpmlcoeff_init(bndr,dx,nb)
  call rsf_read(Fv,v0)
  call check_sanity(maxval(v0),dt,dz,dx,dy)
  call expand3d(vv, v0, nz, nx, ny, nb)
  call rsf_read(Frho,v0)
  call expand3d(rho, v0, nz, nx, ny, nb)
  call rsf_read(Feta, v0)
  call expand3d(eta, v0, nz, nx, ny, nb)
  p=0.
  vx=0.
  vz=0.
  vy=0.
  conv_pz=0.
  conv_px=0.
  conv_py=0.
  conv_vz=0.
  conv_vx=0.
  conv_vy=0.


  !forward modeling
  do it=1,nt
     write(0,*) it
     if(it==kt) then ! record snapshot at it=kt
        call window3d(v0, p, nz, nx, ny, nb)
        call rsf_write(Fw1,v0)
     endif

     write(0,*) "sum(conv_pz)",sum(abs(conv_pz))
     write(0,*) "sum(conv_px)",sum(abs(conv_px))
     write(0,*) "sum(conv_py)",sum(abs(conv_py))

     if (order1) then ! scheme 1, 1st order accuracy, default
        call step_forward_v(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_vzvxvy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
        call step_forward_p(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_pzpxpy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
        if(attenuating) call apply_attenuation(p, eta, rho, vv, dt, nzpad, nxpad, nypad)
     else !2nd order scheme
        if(attenuating) call apply_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad, nypad)
        call step_forward_v(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_vzvxvy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
        call step_forward_p(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_pzpxpy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
        if(attenuating) call apply_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad, nypad)
     endif
     call add_sources(p, dt, wlt(it), sz, sx, sy, nzpad, nxpad, nypad)
     call boundary_rw(.true.,bvz(:,:,:,:,it),bvx(:,:,:,:,it),bvy(:,:,:,:,it),vz,vx,vy,nz,nx,ny,nb)
  enddo

  !backward reconstruction
  do it=nt,1,-1     
     write(0,*) it
     call boundary_rw(.false.,bvz(:,:,:,:,it),bvx(:,:,:,:,it),bvy(:,:,:,:,it),vz,vx,vy,nz,nx,ny,nb)
     if(it==kt) then ! record snapshot at it=kt
        call window3d(v0, p, nz, nx, ny, nb)
        call rsf_write(Fw2,v0)
     endif
     call add_sources(p, -dt, wlt(it), sz, sx, sy, nzpad, nxpad, nypad)
     if (order1) then ! scheme 1, 1st order accuracy, default
        if(attenuating) call apply_attenuation(p, eta, rho, vv, -dt, nzpad, nxpad, nypad)
        call update_cpml_pzpxpy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, -dt, nz, nx, ny, nb)
        call step_forward_p(p, vz, vx, vy, vv, rho, -dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_vzvxvy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, -dt, nz, nx, ny, nb)
        call step_forward_v(p, vz, vx, vy, vv, rho, -dt, idz, idx, idy, nzpad, nxpad, nypad)
     else !2nd order scheme
        if(attenuating) call apply_attenuation(p, eta, rho, vv, -0.5*dt, nzpad, nxpad,nypad)
        call update_cpml_pzpxpy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, -dt, nz, nx, ny, nb)
        call step_forward_p(p, vz, vx, vy, vv, rho, -dt, idz, idx, idy, nzpad, nxpad, nypad)
        call update_cpml_vzvxvy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, -dt, nz, nx, ny, nb)
        call step_forward_v(p, vz, vx, vy, vv, rho, -dt, idz, idx, idy, nzpad, nxpad, nypad)
        if(attenuating) call apply_attenuation(p, eta, rho, vv,-0.5*dt, nzpad, nxpad, nypad)
     endif
  enddo

  deallocate(wlt)
  deallocate(bndr)
  deallocate(v0)
  deallocate(vv)
  deallocate(rho)
  deallocate(eta)
  deallocate(p)
  deallocate(vz)
  deallocate(vx)  
  deallocate(vy)
  deallocate(conv_pz)
  deallocate(conv_px)
  deallocate(conv_py)
  deallocate(conv_vz)
  deallocate(conv_vx)
  deallocate(conv_vy)
  deallocate(bvz)
  deallocate(bvx)
  deallocate(bvy)

  call exit(0)
end program mexwell_cpml2_backward

!------------------------------------------------------------------------------
! check the CFL/stability condition is satisfied or not
subroutine check_sanity(vpmax,dt,dz,dx,dy)
  implicit none
  
  real::vpmax,dt,dx,dz,dy
  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  real CFL,tmp

  tmp=c1-c2
  CFL=vpmax*dt/(max(dz,dx,dy)/sqrt(3.*tmp))

  if (CFL>=1) then 
     write(0,*)'do NOT satisfy CFL condition'
     call exit(0)
  else
     write(0,*)'CFL=',CFL
  endif  
end subroutine check_sanity

!------------------------------------------------------------------------------
! expand the model with artificial boundaries
subroutine expand3d(tgt, src, nz, nx, ny, nb)
  implicit none

  integer::i1,i2,i3
  integer::nz,nx,ny,nb,nzpad,nxpad,nypad
  real::tgt(nz+2*nb,nx+2*nb,ny+2*nb),src(nz,nx,ny)

  nzpad=nz+2*nb
  nxpad=nx+2*nb
  nypad=ny+2*nb

  !first copy from source to inner part of target
  do i3=1,ny
     do i2=1,nx
        do i1=1,nz
           tgt(i1+nb,i2+nb,i3+nb)=src(i1,i2,i3)
        enddo
     enddo
  enddo
  !then pad the boundaries
  do i3=1,nypad
     do i2=1,nxpad
        do i1=1,nb
           tgt(i1,i2,i3)=tgt(nb+1,i2,i3)
           tgt(i1+nz+nb,i2,i3)=tgt(nz+nb,i2,i3)
        enddo
     enddo
  enddo
  do i3=1,nypad
     do i2=1,nb
        do i1=1,nzpad
           tgt(i1,i2,i3)=tgt(i1,nb+1,i3)
           tgt(i1,i2+nx+nb,i3)=tgt(i1,nx+nb,i3)
        enddo
     enddo
  enddo
  do i3=1,nb
     do i2=1,nxpad
        do i1=1,nzpad
           tgt(i1,i2,i3)=tgt(i1,i2,nb+1)
           tgt(i1,i2,i3+ny+nb)=tgt(i1,i2,ny+nb)
        enddo
     enddo
  enddo
end subroutine expand3d

!------------------------------------------------------------------------------
!window the inner part from the expanded model
!the source is assumed to be larger in size than the target
subroutine window3d(tgt, src, nz, nx, ny, nb)
  implicit none

  integer::i1, i2, i3, nz, nx, ny, nb
  real::src(nz+2*nb,nx+2*nb,ny+2*nb),tgt(nz,nx,ny)
  do i3=1,ny
     do i2=1,nx
        do i1=1,nz
           tgt(i1,i2,i3)=src(i1+nb,i2+nb,i3+nb)
        enddo
     enddo
  enddo
end subroutine window3d

!-------------------------------------------------------------------------------
subroutine step_forward_v(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
  implicit none

  integer::i1, i2, i3
  real::tmp,diff1,diff2,diff3

  integer:: nzpad, nxpad, nypad
  real::idz,idx,idy, dt
  real,dimension(nzpad,nxpad,nypad)::p, vz, vx, vy, vv, rho

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.
  
  do i3=2,nypad-2
     do i2=2,nxpad-2
        do i1=2,nzpad-2
           diff1=c1*(p(i1+1,i2,i3)-p(i1,i2,i3))&
                +c2*(p(i1+2,i2,i3)-p(i1-1,i2,i3))
           diff2=c1*(p(i1,i2+1,i3)-p(i1,i2,i3))&
                +c2*(p(i1,i2+2,i3)-p(i1,i2-1,i3))
           diff3=c1*(p(i1,i2,i3+1)-p(i1,i2,i3))&
                +c2*(p(i1,i2,i3+2)-p(i1,i2,i3-1))
           vz(i1,i2,i3)=vz(i1,i2,i3)-dt*idz*diff1/rho(i1,i2,i3)
           vx(i1,i2,i3)=vx(i1,i2,i3)-dt*idx*diff2/rho(i1,i2,i3)
           vy(i1,i2,i3)=vy(i1,i2,i3)-dt*idy*diff3/rho(i1,i2,i3)
        enddo
     enddo
  enddo
end subroutine step_forward_v

!------------------------------------------------------------------------------
subroutine step_forward_p(p, vz, vx, vy, vv, rho, dt, idz, idx, idy, nzpad, nxpad, nypad)
  implicit none

  integer::i1, i2,i3
  real::tmp,diff1,diff2,diff3

  integer:: nzpad, nxpad, nypad
  real::idz,idx,idy,dt
  real,dimension(nzpad,nxpad,nypad)::p, vz, vx, vy, vv, rho

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.
  
  do i3=3,nzpad-1
     do i2=3,nxpad-1
        do i1=3,nzpad-1
           tmp=vv(i1,i2,i3)
           tmp=rho(i1,i2,i3)*tmp*tmp
           diff1=c1*(vz(i1,i2,i3)-vz(i1-1,i2,i3))&
                +c2*(vz(i1+1,i2,i3)-vz(i1-2,i2,i3))
           diff2=c1*(vx(i1,i2,i3)-vx(i1,i2-1,i3))&
                +c2*(vx(i1,i2+1,i3)-vx(i1,i2-2,i3))
           diff3=c1*(vy(i1,i2,i3)-vy(i1,i2,i3-1))&
                +c2*(vy(i1,i2,i3+1)-vy(i1,i2,i3-2))
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*(idz*diff1+idx*diff2+idy*diff3)
        enddo
     enddo
  enddo
end subroutine step_forward_p

!------------------------------------------------------------------------------
! initialize PML damping profile
subroutine cpmlcoeff_init(bndr,dx,nb)
  implicit none
  
  integer::nb
  real::dx
  real,dimension(nb)::bndr

  integer::ib
  real::x,L,d0
  real,parameter::Rc=1.e-4
  
  L=nb*dx
  d0=-3.*log(Rc)/(2.*L*L*L)

  do ib=1,nb
     x=(ib-nb)*dx     !x=1.-cos(0.5*(nb-ib)*PI/nb) !1-cos(x)~=0.5*x^2 when x is small
     bndr(ib)=d0*x*x 
  enddo
end subroutine cpmlcoeff_init

!------------------------------------------------------------------------------
subroutine update_cpml_vzvxvy(p, vz, vx, vy, conv_pz, conv_px, conv_py, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
  implicit none

  integer::nz,nx,ny,nb
  real::idz,idx,idy,dt
  real::bndr(nb)
  real,dimension(nz+2*nb,nx+2*nb,ny+2*nb)::p,vz,vx,vy,rho,vv
  real::conv_pz(nb,nx+2*nb,ny+2*nb,2),conv_px(nz+2*nb,nb,ny+2*nb,2),conv_py(nz+2*nb,nx+2*nb,nb,2)

  integer::nzpad,nxpad,nypad,i1,i2,i3,ib
  real*8::b,diff1,diff2,diff3

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  nzpad=nz+2*nb
  nxpad=nx+2*nb
  nypad=ny+2*nb

  !update conv_pz
  do i3=1,nypad
     do i2=1,nxpad
        do i1=2,nb !top
           b=exp(-bndr(i1)*vv(i1,i2,i3)*dt)
           diff1=c1*(p(i1+1,i2,i3)-p(i1,i2,i3)) &
                +c2*(p(i1+2,i2,i3)-p(i1-1,i2,i3))
           conv_pz(i1,i2,i3,1)=b*conv_pz(i1,i2,i3,1)+(b-1.)*diff1*idz
        enddo
        do i1=nz+nb+1,nzpad-2 !bottom
           ib=nzpad-i1+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff1=c1*(p(i1+1,i2,i3)-p(i1,i2,i3)) &
                +c2*(p(i1+2,i2,i3)-p(i1-1,i2,i3))
           conv_pz(ib,i2,i3,2)=b*conv_pz(ib,i2,i3,2)+(b-1.)*diff1*idz
        enddo
     enddo
  enddo
  !update conv_px
  do i3=1,nypad
     do i1=1,nzpad
        do i2=2,nb !left
           b=exp(-bndr(i2)*vv(i1,i2,i3)*dt)
           diff2=c1*(p(i1,i2+1,i3)-p(i1,i2,i3))&
                +c2*(p(i1,i2+2,i3)-p(i1,i2-1,i3)) 
           conv_px(i1,i2,i3,1)=b*conv_px(i1,i2,i3,1)+(b-1.)*diff2*idx
        enddo
        do i2=nx+nb+1,nxpad-2 !right
           ib=nxpad-i2+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff2=c1*(p(i1,i2+1,i3)-p(i1,i2,i3)) &
                +c2*(p(i1,i2+2,i3)-p(i1,i2-1,i3))
           conv_px(i1,ib,i3,2)=b*conv_px(i1,ib,i3,2)+(b-1.)*diff2*idx
        enddo
     enddo
  enddo
  !update conv_py
  do i2=1,nxpad
     do i1=1,nzpad
        do i3=2,nb
           b=exp(-bndr(i3)*vv(i1,i2,i3)*dt)
           diff3=c1*(p(i1,i2,i3+1)-p(i1,i2,i3))&
                +c2*(p(i1,i2,i3+2)-p(i1,i2,i3-1))
           conv_py(i1,i2,i3,1)=b*conv_py(i1,i2,i3,1)+(b-1.)*diff3*idy
        enddo
        do i3=ny+nb+1,nypad-2
           ib=nypad-i3+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff3=c1*(p(i1,i2,i3+1)-p(i1,i2,i3))&
                +c2*(p(i1,i2,i3+2)-p(i1,i2,i3-1))
           conv_py(i1,i2,ib,2)=b*conv_py(i1,i2,ib,2)+(b-1.)*diff3*idy
        enddo
     enddo
  enddo

  !update vz
  do i3=1,nypad
     do i2=1,nxpad
        do i1=1,nb !top
           vz(i1,i2,i3)=vz(i1,i2,i3)-dt*conv_pz(i1,i2,i3,1)/rho(i1,i2,i3)
        enddo
        do i1=nz+nb+1,nzpad !bottom
           ib=nzpad-i1+1
           vz(i1,i2,i3)=vz(i1,i2,i3)-dt*conv_pz(ib,i2,i3,2)/rho(i1,i2,i3)
        enddo
     enddo
  enddo
  !update vx
  do i3=1,nypad
     do i1=1,nzpad
        do i2=1,nb !left
           vx(i1,i2,i3)=vx(i1,i2,i3)-dt*conv_px(i1,i2,i3,1)/rho(i1,i2,i3)
        enddo
        do i2=nx+nb+1,nxpad !right
           ib=nxpad-i2+1
           vx(i1,i2,i3)=vx(i1,i2,i3)-dt*conv_px(i1,ib,i3,2)/rho(i1,i2,i3)
        enddo
     enddo
  enddo
  !update vy
  do i2=1,nxpad
     do i1=1,nzpad
        do i3=1,nb
           vy(i1,i2,i3)=vy(i1,i2,i3)-dt*conv_py(i1,i2,i3,1)/rho(i1,i2,i3)
        enddo
        do i3=ny+nb+1,nypad
           ib=nypad-i3+1
           vy(i1,i2,i3)=vy(i1,i2,i3)-dt*conv_py(i1,i2,ib,2)/rho(i1,i2,i3)
        enddo
     enddo
  enddo
end subroutine update_cpml_vzvxvy

!------------------------------------------------------------------------------
subroutine update_cpml_pzpxpy(p, vz, vx, vy, conv_vz, conv_vx, conv_vy, rho, vv, bndr, idz, idx, idy, dt, nz, nx, ny, nb)
  implicit none
  
  integer::nz,nx,ny,nb
  real::idz,idx,idy,dt
  real,dimension(nb)::bndr
  real,dimension(nz+2*nb,nx+2*nb,ny+2*nb)::p,vz,vx,vy,rho,vv
  real::conv_vz(nb,nx+2*nb,ny+2*nb,2),conv_vx(nz+2*nb,nb,ny+2*nb,2),conv_vy(nz+2*nb,nx+2*nb,nb,2)

  integer::i1,i2,i3,ib,nzpad,nxpad,nypad
  real*8::diff1,diff2,diff3,b,tmp

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  nzpad=nz+2*nb
  nxpad=nx+2*nb
  nypad=ny+2*nb
  
  !update conv_vz
  do i3=1,nypad
     do i2=1,nxpad
        do i1=3,nb !top
           b=exp(-bndr(i1)*vv(i1,i2,i3)*dt)
           diff1=c1*(vz(i1,i2,i3)-vz(i1-1,i2,i3))&
                +c2*(vz(i1+1,i2,i3)-vz(i1-2,i2,i3))
           conv_vz(i1,i2,i3,1)=b*conv_vz(i1,i2,i3,1)+(b-1.)*diff1*idz
        enddo
        do i1=nz+nb+1,nzpad-1 !bottom
           ib=nzpad-i1+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff1=c1*(vz(i1,i2,i3)-vz(i1-1,i2,i3))&
                +c2*(vz(i1+1,i2,i3)-vz(i1-2,i2,i3))
           conv_vz(ib,i2,i3,2)=b*conv_vz(ib,i2,i3,2)+(b-1.)*diff1*idz
        enddo
     enddo
  enddo
  !update conv_vx
  do i3=1,nypad
     do i1=1,nzpad
        do i2=3,nb !left
           b=exp(-bndr(i2)*vv(i1,i2,i3)*dt)
           diff2=c1*(vx(i1,i2,i3)-vx(i1,i2-1,i3))&
                +c2*(vx(i1,i2+1,i3)-vx(i1,i2-2,i3))
           conv_vx(i1,i2,i3,1)=b*conv_vx(i1,i2,i3,1)+(b-1.)*diff2*idx
        enddo
        do i2=nx+nb+1,nxpad-1 !right
           ib=nxpad-i2+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff2=c1*(vx(i1,i2,i3)-vx(i1,i2-1,i3))&
                +c2*(vx(i1,i2+1,i3)-vx(i1,i2-2,i3))
           conv_vx(i1,ib,i3,2)=b*conv_vx(i1,ib,i3,2)+(b-1.)*diff2*idx
        enddo
     enddo
  enddo
  !update conv_vy
  do i2=1,nxpad
     do i1=1,nzpad
        do i3=3,nb
           b=exp(-bndr(i3)*vv(i1,i2,i3)*dt)
           diff3=c1*(vy(i1,i2,i3)-vy(i1,i2,i3-1))&
                +c2*(vy(i1,i2,i3+1)-vy(i1,i2,i3-2))
           conv_vy(i1,i2,i3,1)=b*conv_vy(i1,i2,i3,1)+(b-1.)*diff3*idy
        enddo
        do i3=ny+nb+1,nypad-1
           ib=nypad-i3+1
           b=exp(-bndr(ib)*vv(i1,i2,i3)*dt)
           diff3=c1*(vy(i1,i2,i3)-vy(i1,i2,i3-1))&
                +c2*(vy(i1,i2,i3+1)-vy(i1,i2,i3-2))
           conv_vy(i1,i2,ib,2)=b*conv_vy(i1,i2,ib,2)+(b-1.)*diff3*idy
        enddo
     enddo
  enddo

  !update pz
  do i3=1,nypad
     do i2=1,nxpad
        do i1=1,nb !top
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vz(i1,i2,i3,1)
        enddo
        do i1=nz+nb+1,nzpad !bottom
           ib=nzpad-i1+1
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vz(ib,i2,i3,2)
        enddo
     enddo
  enddo
  !update px
  do i3=1,nypad
     do i1=1,nzpad
        do i2=1,nb !left
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vx(i1,i2,i3,1)
        enddo
        do i2=nx+nb+1,nxpad !right
           ib=nxpad-i2+1
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vx(i1,ib,i3,2)
        enddo
     enddo
  enddo
  !update py
  do i2=1,nxpad
     do i1=1,nzpad
        do i3=1,nb
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vy(i1,i2,i3,1)
        enddo
        do i3=ny+nb+1,nypad
           ib=nypad-i3+1
           tmp=vv(i1,i2,i3);        tmp=rho(i1,i2,i3)*tmp*tmp
           p(i1,i2,i3)=p(i1,i2,i3)-dt*tmp*conv_vy(i1,i2,ib,2)
        enddo
     enddo
  enddo
end subroutine update_cpml_pzpxpy

!------------------------------------------------------------------------------
subroutine apply_attenuation(p, eta, rho, vv, dt, nzpad, nxpad,nypad)
  implicit none

  integer::i1,i2,i3
  real*8::a,tau

  integer::nzpad,nxpad,nypad
  real::dt
  real, dimension(nzpad,nxpad,nypad)::p, eta, rho, vv
  
  do i3=1,nypad
     do i2=1,nxpad
        do i1=1,nzpad
           a=rho(i1,i2,i3)*vv(i1,i2,i3)*vv(i1,i2,i3)
           tau=eta(i1,i2,i3)/a
           a=exp(-dt/tau)
           p(i1,i2,i3)=a*p(i1,i2,i3)
        enddo
     enddo
  enddo
end subroutine apply_attenuation

!-------------------------------------------------------------------------------
subroutine add_sources(p, dt, wlt, sz, sx, sy,nzpad, nxpad, nypad)
  implicit none

  integer::sz,sx,sy,nzpad, nxpad,nypad
  real::dt,wlt
  real,dimension(nzpad, nxpad,nypad)::p

  p(sz,sx,sy)=p(sz,sx,sy)+dt*wlt
end subroutine add_sources

!-------------------------------------------------------------------------------
subroutine boundary_rw(v2b, bvz, bvx, bvy, vz, vx, vy, nz, nx, ny, nb)
  implicit none
  
  logical::v2b !v to bourndary or reverse
  integer::nz,nx,ny,nb
  real,dimension(nz+2*nb,nx+2*nb,ny+2*nb)::vz,vx,vy
  real::bvz(3,nx,ny,2),bvx(nz,3,ny,2),bvy(nz,nx,3,2)

  integer::i1,i2,i3
  
  if(v2b) then !v to bourndary
     do i3=1,ny
        do i2=1,nx
           do i1=1,3
              bvz(i1,i2,i3,1)=vz(i1+nb-2,i2+nb,i3+nb)
              bvz(i1,i2,i3,2)=vz(i1+nz+nb-2,i2+nb,i3+nb)
           enddo
        enddo
     enddo
     do i3=1,ny
        do i1=1,nz
           do i2=1,3
              bvx(i1,i2,i3,1)=vx(i1+nb,i2+nb-2,i3+nb)
              bvx(i1,i2,i3,2)=vx(i1+nb,i2+nz+nb-2,i3+nb)
           enddo
        enddo
     enddo
     do i2=1,nx
        do i1=1,nz
           do i3=1,3
              bvy(i1,i2,i3,1)=vy(i1+nb,i2+nb,i3+nb-2)
              bvy(i1,i2,i3,2)=vy(i1+nb,i2+nb,i3+nz+nb-2)
           enddo
        enddo
     enddo
  else !boundary to v
     do i3=1,ny
        do i2=1,nx
           do i1=1,3
              vz(i1+nb-2,i2+nb,i3+nb)=bvz(i1,i2,i3,1)
              vz(i1+nz+nb-2,i2+nb,i3+nb)=bvz(i1,i2,i3,2)
           enddo
        enddo
     enddo
     do i3=1,ny
        do i1=1,nz
           do i2=1,3
              vx(i1+nb,i2+nb-2,i3+nb)=bvx(i1,i2,i3,1)
              vx(i1+nb,i2+nz+nb-2,i3+nb)=bvx(i1,i2,i3,2)
           enddo
        enddo
     enddo
     do i2=1,nx
        do i1=1,nz
           do i3=1,3
              vy(i1+nb,i2+nb,i3+nb-2)=bvy(i1,i2,i3,1)
              vy(i1+nb,i2+nb,i3+ny+nb-2)=bvy(i1,i2,i3,2)
           enddo
        enddo
     enddo
  endif
end subroutine boundary_rw

