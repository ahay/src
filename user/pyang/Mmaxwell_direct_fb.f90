!  Forward modeling based upon Maxwell model and sponge ABC
!  This code is only valid with attenuation. No attenuation fails.
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
program mexwell_direct_backward
  use rsf
  implicit none

  integer :: ib, it, nt, nz, nx, nb, sx, sz, nxpad, nzpad
  real    :: dt, dz, dx, fm, tmp, idx, idz
  real*8, parameter::PI=4.*atan(1.)
  real, dimension (:),   allocatable :: wlt,bndr,ef,eb, wltf, wltb
  real, dimension (:,:), allocatable :: v0, vv, rho, eta
  real, dimension (:,:), allocatable :: p, vz, vx
  real, dimension(:,:,:,:),allocatable :: bvz,bvx
  
  type(file) :: Fv, Frho, Feta, Fw1, Fw2, Fef,Feb, Fwltf, Fwltb  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fv = rsf_input ("in")   ! velocity model
  Fw1 = rsf_output("out") ! output forward wavefield 
  Fw2 = rsf_output("back")! output forward wavefield 
  Frho=rsf_input("rho")   ! density
  Feta=rsf_input("eta")   ! Pascal
  Fef=rsf_output("ef")    ! energy (kinematic+deformation) when forward modeling
  Feb=rsf_output("eb")    ! energy (kinematic+deformation) when back-reconstruction
  Fwltf=rsf_output("wltf") ! wavelet in forward modeling
  Fwltb=rsf_output("wltb") ! wavelet in backward reconstruction

  ! Read/Write axes
  call from_par(Fv,"n1",nz) ! velocity model: nz
  call from_par(Fv,"n2",nx) ! velocity model: nx
  call from_par(Fv,"d1",dz) ! velocity model: dz
  call from_par(Fv,"d2",dx) ! velocity model: dx

  call from_par("nb", nb, 30) ! thinkness of sponge ABC
  call from_par("nt", nt, 1000) !number of time steps
  call from_par("dt", dt, 0.001) ! time sampling interval
  call from_par("fm", fm, 20.) ! domainant frequency for ricker wavelet

  call to_par(Fw1,"n1",nz)
  call to_par(Fw1,"n2",nx)
  call to_par(Fw1,"d1",dz)
  call to_par(Fw1,"d2",dx)
  call to_par(Fw1,"n3",nt)
  call to_par(Fw1,"o3",0)
  call to_par(Fw1,"d3",dt)

  call to_par(Fw2,"n1",nz)
  call to_par(Fw2,"n2",nx)
  call to_par(Fw2,"d1",dz)
  call to_par(Fw2,"d2",dx)
  call to_par(Fw2,"n3",nt)
  call to_par(Fw2,"o3",(nt-1)*dt)
  call to_par(Fw2,"d3",-dt)

  call to_par(Fef,"n1",nt)
  call to_par(Fef,"o1",0)
  call to_par(Fef,"d1",dt)
  call to_par(Fef,"n2",1)

  call to_par(Feb,"n1",nt)
  call to_par(Feb,"o1",0)
  call to_par(Feb,"d1",dt)
  call to_par(Feb,"n2",1)

  call to_par(Fwltf,"n1",nt)
  call to_par(Fwltf,"o1",0)
  call to_par(Fwltf,"d1",dt)
  call to_par(Fwltf,"n2",1)

  call to_par(Fwltb,"n1",nt)
  call to_par(Fwltb,"o1",0)
  call to_par(Fwltb,"d1",dt)
  call to_par(Fwltb,"n2",1)

  idx=1./dx
  idz=1./dz
  nzpad=nz+2*nb
  nxpad=nx+2*nb
  sx=nxpad/2
  sz=nzpad/2

  allocate(wlt(nt))
  allocate(bndr(nb))
  allocate(v0(nz,nx))
  allocate(vv(nzpad,nxpad))
  allocate(rho(nzpad,nxpad))
  allocate(eta(nzpad,nxpad))
  allocate(p(nzpad,nxpad))
  allocate(vz(nzpad,nxpad))
  allocate(vx(nzpad,nxpad))
  allocate(bvz(7,nx,2,nt))
  allocate(bvx(nz,7,2,nt))
  allocate(ef(nt))
  allocate(eb(nt))
  allocate(wltf(nt))
  allocate(wltb(nt))

  !generate ricker wavelet with a delay
  do it=1,nt  
     tmp=PI*fm*(it*dt-1.0/fm)
     tmp=tmp*tmp
     wlt(it)=(1.0-2.0*tmp)*exp(-tmp)
  enddo
  !generate coefficients for the absorbing boundary
  do ib=1,nb 
     tmp=0.015*(nb-ib)
     bndr(ib)=exp(-tmp*tmp)
  enddo
  call rsf_read(Fv,v0)
  call check_sanity(maxval(v0),dt,dx,dz)
  call expand2d(vv, v0, nz, nx, nb)
  call rsf_read(Frho,v0)
  call expand2d(rho, v0, nz, nx, nb)
  call rsf_read(Feta, v0)
  call expand2d(eta, v0, nz, nx, nb)
  p=0.
  vx=0.
  vz=0.

  !forward modeling
  do it=1,nt
     call step_forward(p, vz, vx, vv, rho, eta, dt, idz, idx, nzpad, nxpad)
     call add_sources(p, dt, wlt(it), sz, sx, nzpad, nxpad,wltf(it))

     ! apply sponge ABC
     call apply_sponge(p, bndr,nz,nx,nb)
     call apply_sponge(vz,bndr,nz,nx,nb)
     call apply_sponge(vx,bndr,nz,nx,nb)

     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw1,v0)

     call compute_energy(ef(it),vz(nb+1:nb+nz,nb+1:nb+nx),&
          vx(nb+1:nb+nz,nb+1:nb+nx),p(nb+1:nb+nz,nb+1:nb+nx),&
          rho(nb+1:nb+nz,nb+1:nb+nx),vv(nb+1:nb+nz,nb+1:nb+nx),nz,nx)
     call boundary_rw(.true.,bvz(:,:,:,it),bvx(:,:,:,it),vz,vx,nz,nx,nb)
  enddo
  call rsf_write(Fef,ef) ! store forward energy history
  call rsf_write(Fwltf,wltf) ! store wavelet series during forward modeling

  !backward reconstruction
  do it=nt,1,-1
     call boundary_rw(.false.,bvz(:,:,:,it),bvx(:,:,:,it),vz,vx,nz,nx,nb)
     call compute_energy(eb(it),vz(nb+1:nb+nz,nb+1:nb+nx),&
          vx(nb+1:nb+nz,nb+1:nb+nx),p(nb+1:nb+nz,nb+1:nb+nx),&
          rho(nb+1:nb+nz,nb+1:nb+nx),vv(nb+1:nb+nz,nb+1:nb+nx),nz,nx)

     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw2,v0)
     
     call add_sources(p, -dt, wlt(it), sz, sx, nzpad, nxpad,wltb(it))
     call step_forward(p, vz, vx, vv, rho, eta, -dt, idz, idx, nzpad, nxpad)
  enddo
  call rsf_write(Feb,eb) !store backward energy history
  call rsf_write(Fwltb,wltb) ! store wavelet series during backward reconstruction  

  deallocate(wlt)
  deallocate(bndr)
  deallocate(v0)
  deallocate(vv)
  deallocate(rho)
  deallocate(eta)
  deallocate(p)
  deallocate(vz)
  deallocate(vx)  
  deallocate(bvz)
  deallocate(bvx)
  deallocate(ef)
  deallocate(eb)
  deallocate(wltf)
  deallocate(wltb)

  call exit(0)
end program mexwell_direct_backward

!--------------------------------------------------------------------------------
! check the CFL/stability condition is satisfied or not
subroutine check_sanity(vpmax,dt,dx,dz)
  implicit none
  
  real::vpmax,dt,dx,dz
  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

  real CFL,tmp

  tmp=c1-c2+c3-c4
  CFL=vpmax*dt/(max(dx,dz)/sqrt(2.*tmp))

  if (CFL>=1) then 
     write(0,*)'do NOT satisfy CFL condition'
     call exit(0)
  else
     write(0,*)'CFL=',CFL
  endif  
end subroutine check_sanity

!------------------------------------------------------------------------------
! expand the model with artificial boundaries
subroutine expand2d(tgt, src, nz, nx, nb)
  implicit none

  integer::i1,i2
  integer::nz,nx,nb,nzpad,nxpad
  real::tgt(nz+2*nb,nx+2*nb),src(nz,nx)

  nzpad=nz+2*nb
  nxpad=nx+2*nb

  !first copy from source to inner part of target
  do i2=1,nx
     do i1=1,nz
        tgt(i1+nb,i2+nb)=src(i1,i2)
     enddo
  enddo
  !then pad the boundaries
  do i2=1,nxpad
     do i1=1,nb
        tgt(i1,i2)=tgt(nb+1,i2)
        tgt(i1+nz+nb,i2)=tgt(nz+nb,i2)
     enddo
  enddo
  do i2=1,nb
     do i1=1,nzpad
        tgt(i1,i2)=tgt(i1,nb+1)
        tgt(i1,i2+nx+nb)=tgt(i1,nx+nb)
     enddo
  enddo
end subroutine expand2d

!------------------------------------------------------------------------------
!window the inner part from the expanded model
!the source is assumed to be larger in size than the target
subroutine window2d(tgt, src, nz, nx, nb)
  implicit none

  integer::i1, i2, nz, nx, nb
  real::src(nz+2*nb,nx+2*nb),tgt(nz,nx)

  do i2=1,nx
     do i1=1,nz
        tgt(i1,i2)=src(i1+nb,i2+nb)
     enddo
  enddo

  return
end subroutine window2d

!-------------------------------------------------------------------------------
subroutine step_forward(p, vz, vx, vv, rho, eta, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2,tau

  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho, eta

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

  if(dt>0) then
     do i2=4,nxpad-4
        do i1=4,nzpad-4
           diff1=c1*(p(i1+1,i2)-p(i1,i2))&
                +c2*(p(i1+2,i2)-p(i1-1,i2))&
                +c3*(p(i1+3,i2)-p(i1-2,i2))&
                +c4*(p(i1+4,i2)-p(i1-3,i2))
           diff2=c1*(p(i1,i2+1)-p(i1,i2))&
                +c2*(p(i1,i2+2)-p(i1,i2-1))&
                +c3*(p(i1,i2+3)-p(i1,i2-2))&
                +c4*(p(i1,i2+4)-p(i1,i2-3))
           vz(i1,i2)=vz(i1,i2)-dt*idz*diff1/rho(i1,i2)
           vx(i1,i2)=vx(i1,i2)-dt*idx*diff2/rho(i1,i2)
        enddo
     enddo

     do i2=5,nxpad-3
        do i1=5,nzpad-3
           tmp=vv(i1,i2)
           tmp=rho(i1,i2)*tmp*tmp
           tau=eta(i1,i2)/tmp
           diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
                +c2*(vz(i1+1,i2)-vz(i1-2,i2))&
                +c3*(vz(i1+2,i2)-vz(i1-3,i2))&
                +c4*(vz(i1+3,i2)-vz(i1-4,i2))
           diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
                +c2*(vx(i1,i2+1)-vx(i1,i2-2))&
                +c3*(vx(i1,i2+2)-vx(i1,i2-3))&
                +c4*(vx(i1,i2+3)-vx(i1,i2-4))
           p(i1,i2)=((1.-0.5*dt/tau)*p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2))/(1.+0.5*dt/tau)
        enddo
     enddo
  else
     do i2=5,nxpad-3
        do i1=5,nzpad-3
           tmp=vv(i1,i2)
           tmp=rho(i1,i2)*tmp*tmp
           tau=eta(i1,i2)/tmp
           diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
                +c2*(vz(i1+1,i2)-vz(i1-2,i2))&
                +c3*(vz(i1+2,i2)-vz(i1-3,i2))&
                +c4*(vz(i1+3,i2)-vz(i1-4,i2))
           diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
                +c2*(vx(i1,i2+1)-vx(i1,i2-2))&
                +c3*(vx(i1,i2+2)-vx(i1,i2-3))&
                +c4*(vx(i1,i2+3)-vx(i1,i2-4))
           p(i1,i2)=((1.-0.5*dt/tau)*p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2))/(1.+0.5*dt/tau)
        enddo
     enddo
     do i2=4,nxpad-4
        do i1=4,nzpad-4
           diff1=c1*(p(i1+1,i2)-p(i1,i2))&
                +c2*(p(i1+2,i2)-p(i1-1,i2))&
                +c3*(p(i1+3,i2)-p(i1-2,i2))&
                +c4*(p(i1+4,i2)-p(i1-3,i2))
           diff2=c1*(p(i1,i2+1)-p(i1,i2))&
                +c2*(p(i1,i2+2)-p(i1,i2-1))&
                +c3*(p(i1,i2+3)-p(i1,i2-2))&
                +c4*(p(i1,i2+4)-p(i1,i2-3))
           vz(i1,i2)=vz(i1,i2)-dt*idz*diff1/rho(i1,i2)
           vx(i1,i2)=vx(i1,i2)-dt*idx*diff2/rho(i1,i2)
        enddo
     enddo
  endif
end subroutine step_forward

!-------------------------------------------------------------------------------
subroutine add_sources(p, dt, wlt, sz, sx, nzpad, nxpad, wlt_actual)
  implicit none

  integer::sz,sx,nzpad, nxpad
  real::dt,wlt,wlt_actual
  real,dimension(nzpad, nxpad)::p

  wlt_actual=dt*wlt
  p(sz,sx)=p(sz,sx)+wlt_actual
end subroutine add_sources

!-------------------------------------------------------------------------------
subroutine apply_sponge(p, bndr, nz, nx, nb)
  implicit none

  integer::nb,nz,nx
  real,dimension(nb)::bndr
  real,dimension(nz+2*nb,nx+2*nb)::p

  integer::i1,i2,nzpad,nxpad
  nzpad=nz+2*nb
  nxpad=nx+2*nb

  do i2=1,nxpad
     do i1=1,nb !top
        p(i1,i2)=bndr(i1)*p(i1,i2)
     enddo
     do i1=nz+nb+1,nzpad !bottom
        p(i1,i2)=bndr(nzpad-i1+1)*p(i1,i2)
     enddo
  enddo
  do i1=1,nzpad
     do i2=1,nb !left
        p(i1,i2)=bndr(i2)*p(i1,i2)
     enddo
     do i2=nb+nx+1,nxpad !right
        p(i1,i2)=bndr(nxpad-i2+1)*p(i1,i2)
     enddo
  enddo

  return
end subroutine apply_sponge


!-------------------------------------------------------------------------------
subroutine boundary_rw(v2b,bvz,bvx,vz,vx,nz,nx,nb)
  implicit none
  
  logical::v2b !v to bourndary or reverse
  integer::nz,nx,nb
  real,dimension(nz+2*nb,nx+2*nb)::vz,vx
  real::bvz(7,nx,2),bvx(nz,7,2)

  integer::i1,i2
  
  if(v2b) then !v to bourndary
     do i2=1,nx
        do i1=1,7
           bvz(i1,i2,1)=vz(i1+nb-4,i2+nb)
           bvz(i1,i2,2)=vz(i1+nz+nb-4,i2+nb)
        enddo
     enddo
     do i1=1,nz
        do i2=1,7
           bvx(i1,i2,1)=vx(i1+nb,i2+nb-4)
           bvx(i1,i2,2)=vx(i1+nb,i2+nz+nb-4)
        enddo
     enddo
  else !boundary to v
     do i2=1,nx
        do i1=1,7
           vz(i1+nb-4,i2+nb)=bvz(i1,i2,1)
           vz(i1+nz+nb-4,i2+nb)=bvz(i1,i2,2)
        enddo
     enddo
     do i1=1,nz
        do i2=1,7
           vx(i1+nb,i2+nb-4)=bvx(i1,i2,1)
           vx(i1+nb,i2+nz+nb-4)=bvx(i1,i2,2)
        enddo
     enddo
  endif
end subroutine boundary_rw

!-------------------------------------------------------------------------------
subroutine compute_energy(e,vz,vx,p,rho,vv,nz,nx)
  implicit none

  integer::nz,nx
  real::e
  real,dimension(nz,nx)::rho,vv,p,vz,vx

  e=0.5*sum((vz*vz+vx*vx)*rho+p*p/(rho*vv*vv))
end subroutine compute_energy
