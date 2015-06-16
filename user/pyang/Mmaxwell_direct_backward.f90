!  Forward modeling based upon Maxwell model and sponge ABC
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
  real, parameter::PI=3.14159265
  real, dimension (:),   allocatable :: wlt,bndr
  real, dimension (:,:), allocatable :: v0, vv, rho, eta
  real, dimension (:,:), allocatable :: p, vz, vx
  real, dimension(:,:,:,:),allocatable :: bndrtb,bndrlr
  
  type(file) :: Fv, Frho, Feta, Fw1, Fw2  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fv = rsf_input ("in")   ! source position 
  Fw1 = rsf_output("out")  ! output forward wavefield 
  Fw2 = rsf_output("back")  ! output forward wavefield 
  Frho=rsf_input("rho")   ! density
  Feta=rsf_input("eta")   ! Pascal

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
  allocate(bndrtb(7,nx,2,nt))
  allocate(bndrlr(nz,7,2,nt))

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
     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw1,v0)

     call step_forward(.true.,p, vz, vx, vv, rho, eta, dt, idz, idx, nzpad, nxpad)
     call add_sources(p, eta, rho, vv, dt, wlt(it), sz, sx, nzpad, nxpad)

     ! apply sponge ABC
     call apply_sponge(p,bndr,nz,nx,nb)
     call apply_sponge(vz,bndr,nz,nx,nb)
     call apply_sponge(vx,bndr,nz,nx,nb)

     ! save the boundaries
     !call boundary_rw(.true.,p,bndrtb(:,:,:,it),bndrlr(:,:,:,it),nz,nx,nb)
  enddo

  !backward reconstruction
  do it=nt,1,-1

     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw2,v0)

     !read the saved boundaries into p
     !call boundary_rw(.false.,p,bndrtb(:,:,:,it),bndrlr(:,:,:,it),nz,nx,nb)
     
     call add_sources(p, eta, rho, vv, -dt, wlt(it), sz, sx, nzpad, nxpad)
     call step_forward(.false.,p, vz, vx, vv, rho, eta, -dt, idz, idx, nzpad, nxpad)
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
  deallocate(bndrtb)
  deallocate(bndrlr)

  call exit(0)
end program mexwell_direct_backward

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
subroutine step_forward(forw, p, vz, vx, vv, rho, eta, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2,tau

  logical::forw
  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho, eta

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857
  if(forw) then
     do i2=4,nxpad-4
        do i1=4,nzpad-4
           diff1=c1*(p(i1+1,i2)-p(i1,i2))+c2*(p(i1+2,i2)-p(i1-1,i2)) &
                +c3*(p(i1+3,i2)-p(i1-2,i2))+c4*(p(i1+4,i2)-p(i1-3,i2))
           diff2=c1*(p(i1,i2+1)-p(i1,i2))+c2*(p(i1,i2+2)-p(i1,i2-1)) &
                +c3*(p(i1,i2+3)-p(i1,i2-2))+c4*(p(i1,i2+4)-p(i1,i2-3))
           vz(i1,i2)=vz(i1,i2)-dt*idz*diff1/rho(i1,i2)
           vx(i1,i2)=vx(i1,i2)-dt*idx*diff2/rho(i1,i2)
        enddo
     enddo

     do i2=5,nxpad-3
        do i1=5,nzpad-3
           tmp=vv(i1,i2)
           tmp=rho(i1,i2)*tmp*tmp
           tau=eta(i1,i2)/tmp
           diff1=c1*(vz(i1,i2)-vz(i1-1,i2))+c2*(vz(i1+1,i2)-vz(i1-2,i2)) &
                +c3*(vz(i1+2,i2)-vz(i1-1,i2))+c4*(vz(i1+3,i2)-vz(i1-2,i2))
           diff2=c1*(vx(i1,i2)-vx(i1,i2-1))+c2*(vx(i1,i2+1)-vx(i1,i2-2)) &
                +c3*(vx(i1,i2+2)-vx(i1,i2-3))+c4*(vx(i1,i2+3)-vx(i1,i2-4))
           p(i1,i2)=((1-0.5*dt/tau)*p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2))/(1.+0.5*dt/tau)
        enddo
     enddo
  else
     do i2=5,nxpad-3
        do i1=5,nzpad-3
           tmp=vv(i1,i2)
           tmp=rho(i1,i2)*tmp*tmp
           tau=eta(i1,i2)/tmp
           diff1=c1*(vz(i1,i2)-vz(i1-1,i2))+c2*(vz(i1+1,i2)-vz(i1-2,i2)) &
                +c3*(vz(i1+2,i2)-vz(i1-1,i2))+c4*(vz(i1+3,i2)-vz(i1-2,i2))
           diff2=c1*(vx(i1,i2)-vx(i1,i2-1))+c2*(vx(i1,i2+1)-vx(i1,i2-2)) &
                +c3*(vx(i1,i2+2)-vx(i1,i2-3))+c4*(vx(i1,i2+3)-vx(i1,i2-4))
           p(i1,i2)=((1-0.5*dt/tau)*p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2))/(1.+0.5*dt/tau)
        enddo
     enddo
     do i2=4,nxpad-4
        do i1=4,nzpad-4
           diff1=c1*(p(i1+1,i2)-p(i1,i2))+c2*(p(i1+2,i2)-p(i1-1,i2)) &
                +c3*(p(i1+3,i2)-p(i1-2,i2))+c4*(p(i1+4,i2)-p(i1-3,i2))
           diff2=c1*(p(i1,i2+1)-p(i1,i2))+c2*(p(i1,i2+2)-p(i1,i2-1)) &
                +c3*(p(i1,i2+3)-p(i1,i2-2))+c4*(p(i1,i2+4)-p(i1,i2-3))
           vz(i1,i2)=vz(i1,i2)-dt*idz*diff1/rho(i1,i2)
           vx(i1,i2)=vx(i1,i2)-dt*idx*diff2/rho(i1,i2)
        enddo
     enddo
  endif
  return
end subroutine step_forward

!-------------------------------------------------------------------------------
subroutine add_sources(p, eta, rho, vv, dt, wlt, sz, sx, nzpad, nxpad)
  implicit none

  integer::sz,sx,nzpad, nxpad
  real::dt,wlt
  real,dimension(nzpad, nxpad)::p, eta, rho, vv

  real::a, tau

  a=rho(sz,sx)*vv(sz,sx)*vv(sz,sx)
  tau=eta(sz,sx)/a
  a=dt/(1.+0.5*dt/tau)
  p(sz,sx)=p(sz,sx)+a*wlt

  return
end subroutine add_sources

!-------------------------------------------------------------------------------
subroutine apply_sponge(p,bndr,nz,nx,nb)
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
subroutine boundary_rw(p2b, p,bndrtb,bndrlr,nz,nx,nb)
  implicit none
  
  logical::p2b !p to bourndary or reverse
  integer::nz,nx,nb
  real::p(nz+2*nb,nx+2*nb),bndrtb(7,nx,2),bndrlr(nz,7,2)

  integer::i1,i2
  
  if(p2b) then !p to bourndary
     do i2=1,nx
        do i1=1,7
           bndrtb(i1,i2,1)=p(nb-i1+1,i2+nb) !top
           bndrtb(i1,i2,2)=p(nb+nz+i1,i2+nb) !bottom
        enddo
     enddo

     do i2=1,7
        do i1=1,nz
           bndrlr(i1,i2,1)=p(i1+nb,nb-i2+1) !left
           bndrlr(i1,i2,2)=p(i1+nb,nb+nz+i2) !right
        enddo
     enddo
  else !boundary to p
     do i2=1,nx
        do i1=1,7
           p(nb-i1+1,i2+nb)=bndrtb(i1,i2,1) !top
           p(nb+nz+i1,i2+nb)=bndrtb(i1,i2,2)!bottom
        enddo
     enddo

     do i2=1,7
        do i1=1,nz
           p(i1+nb,nb-i2+1)=bndrlr(i1,i2,1) !left
           p(i1+nb,nb+nz+i2)=bndrlr(i1,i2,2) !right
        enddo
     enddo
  endif
end subroutine boundary_rw
