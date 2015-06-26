!  Forward modeling based upon Maxwell attenuation model and CPML ABC
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
program mexwell_cpml2
  use rsf
  implicit none

  logical :: order1
  integer :: ib, it, nt, nz, nx, nb, sx, sz, nxpad, nzpad
  real    :: dt, dz, dx, fm, tmp, idx, idz
  real, parameter::PI=3.14159265
  real, dimension (:),   allocatable :: wlt,bndr
  real, dimension (:,:), allocatable :: v0, vv, rho, eta
  real, dimension (:,:), allocatable :: p, vz, vx
  real, dimension(:,:,:),allocatable :: conv_pz,conv_px,conv_vz,conv_vx
  
  type(file) :: Fw, Fv, Frho, Feta  ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fv = rsf_input ("in")   ! source position 
  Fw = rsf_output("out")  ! output wavefield 
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
  call from_par("order1",order1,.true.) ! 1st order or 2nd order accuracy

  call to_par(Fw,"n1",nz)
  call to_par(Fw,"n2",nx)
  call to_par(Fw,"d1",dz)
  call to_par(Fw,"d2",dx)
  call to_par(Fw,"n3",nt)
  call to_par(Fw,"o3",0)
  call to_par(Fw,"d3",dt)

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
  allocate(conv_pz(nb,nxpad,2))
  allocate(conv_px(nzpad,nb,2))
  allocate(conv_vz(nb,nxpad,2))
  allocate(conv_vx(nzpad,nb,2))

  !generate ricker wavelet with a delay
  do it=1,nt  
     tmp=PI*fm*(it*dt-1.0/fm)
     tmp=tmp*tmp
     wlt(it)=(1.0-2.0*tmp)*exp(-tmp)
  enddo
  !generate coefficients for the absorbing boundary
  call cpmlcoeff_init(bndr,dx,nb)
  call rsf_read(Fv,v0)
  call expand2d(vv, v0, nz, nx, nb)
  call rsf_read(Frho,v0)
  call expand2d(rho, v0, nz, nx, nb)
  call rsf_read(Feta, v0)
  call expand2d(eta, v0, nz, nx, nb)
  p=0.
  vx=0.
  vz=0.
  conv_pz=0.
  conv_px=0.
  conv_vz=0.
  conv_vx=0.

  do it=1,nt
     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw,v0)

     if (order1) then ! scheme 1, 1st order accuracy, default
        call step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bndr,idz,idx,dt,nz,nx,nb)
        call step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_pzpx(p,vz,vx,conv_pz,conv_px,rho,vv,bndr,idz,idx,dt,nz,nx,nb)
        call add_attenuation(p, eta, rho, vv, dt, nzpad, nxpad)
     else 
        call add_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad)

        call step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bndr,idz,idx,dt,nz,nx,nb)
        call step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_pzpx(p,vz,vx,conv_pz,conv_px,rho,vv,bndr,idz,idx,dt,nz,nx,nb)

        call add_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad)
     endif

     call add_sources(p, eta, rho, vv, dt, wlt(it), sz, sx, nzpad, nxpad)
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
  deallocate(conv_pz)
  deallocate(conv_px)
  deallocate(conv_vz)
  deallocate(conv_vx)

  call exit(0)
end program mexwell_cpml2

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
subroutine step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2

  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

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
  return
end subroutine step_forward_v

!------------------------------------------------------------------------------
subroutine step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2

  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

  do i2=5,nxpad-3
     do i1=5,nzpad-3
        tmp=vv(i1,i2)
        tmp=rho(i1,i2)*tmp*tmp
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))&
             +c3*(vz(i1+2,i2)-vz(i1-3,i2))&
             +c4*(vz(i1+3,i2)-vz(i1-4,i2))
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))&
             +c3*(vx(i1,i2+2)-vx(i1,i2-3))&
             +c4*(vx(i1,i2+3)-vx(i1,i2-4))
        p(i1,i2)=p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2)
     enddo
  enddo
  return
end subroutine step_forward_p

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
     x=(ib-nb)*dx
     bndr(ib)=d0*x*x
  enddo
  return
end subroutine cpmlcoeff_init

!------------------------------------------------------------------------------
subroutine update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bndr,idz,idx,dt,nz,nx,nb)
  implicit none

  integer::nz,nx,nb
  real::idz,idx,dt
  real::bndr(nb)
  real,dimension(nz+2*nb,nx+2*nb)::p,vz,vx,rho,vv
  real::conv_pz(nb,nx+2*nb,2),conv_px(nz+2*nb,nb,2)

  integer::nzpad,nxpad,i1,i2,ib
  real::b,diff1,diff2

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

  !update conv_pz
  do i2=1,nxpad
     do i1=4,nb !top
        b=exp(-bndr(i1)*vv(i1,i2)*dt)
        diff1=c1*(p(i1+1,i2)-p(i1,i2)) &
             +c2*(p(i1+2,i2)-p(i1-1,i2)) &
             +c3*(p(i1+3,i2)-p(i1-2,i2)) &
             +c4*(p(i1+4,i2)-p(i1-3,i2))
        conv_pz(i1,i2,1)=b*conv_pz(i1,i2,1)+(b-1.)*diff1*idz
     enddo
     do i1=nz+nb+1,nzpad-4 !bottom
        ib=nzpad-i1+1
        b=exp(-bndr(ib)*vv(i1,i2)*dt)
        diff1=c1*(p(i1+1,i2)-p(i1,i2)) &
             +c2*(p(i1+2,i2)-p(i1-1,i2)) &
             +c3*(p(i1+3,i2)-p(i1-2,i2)) &
             +c4*(p(i1+4,i2)-p(i1-3,i2))
        conv_pz(ib,i2,2)=b*conv_pz(ib,i2,2)+(b-1.)*diff1*idz
     enddo
  enddo

  !update conv_px
  do i1=1,nzpad
     do i2=4,nb !left
        b=exp(-bndr(i2)*vv(i1,i2)*dt)
        diff2=c1*(p(i1,i2+1)-p(i1,i2))&
             +c2*(p(i1,i2+2)-p(i1,i2-1)) &
             +c3*(p(i1,i2+3)-p(i1,i2-2)) &
             +c4*(p(i1,i2+4)-p(i1,i2-3))
        conv_px(i1,i2,1)=b*conv_px(i1,i2,1)+(b-1.)*diff2*idx
     enddo
     do i2=nx+nb+1,nxpad-4 !right
        ib=nxpad-i1+1
        b=exp(-bndr(ib)*vv(i1,i2)*dt)
        diff2=c1*(p(i1,i2+1)-p(i1,i2)) &
             +c2*(p(i1,i2+2)-p(i1,i2-1)) &
             +c3*(p(i1,i2+3)-p(i1,i2-2)) &
             +c4*(p(i1,i2+4)-p(i1,i2-3))
        conv_px(i1,ib,2)=b*conv_px(i1,ib,2)+(b-1.)*diff2*idx
     enddo
  enddo

  !update vz
  do i2=1,nxpad
     do i1=1,nb !top
        vz(i1,i2)=vz(i1,i2)-dt*conv_pz(i1,i2,1)/rho(i1,i2)
     enddo
     do i1=nz+nb+1,nzpad !bottom
        ib=nzpad-i1+1
        vz(i1,i2)=vz(i1,i2)-dt*conv_pz(ib,i2,2)/rho(i1,i2)
     enddo
  enddo

  !update vx
  do i1=1,nzpad
     do i2=1,nb !left
        vx(i1,i2)=vx(i1,i2)-dt*conv_px(i1,i2,1)/rho(i1,i2)
     enddo
     do i2=nx+nb+1,nxpad !right
        ib=nxpad-i2+1
        vx(i1,i2)=vx(i1,i2)-dt*conv_px(i1,i2,2)/rho(i1,i2)
     enddo
  enddo

  return
end subroutine update_cpml_vzvx

!------------------------------------------------------------------------------
subroutine update_cpml_pzpx(p,vz,vx,conv_vz,conv_vx,rho,vv,bndr,idz,idx,dt,nz,nx,nb)
  implicit none
  
  integer::nz,nx,nb
  real::idz,idx,dt
  real,dimension(nb)::bndr
  real,dimension(nz+2*nb,nx+2*nb)::p,vz,vx,rho,vv
  real::conv_vz(nb,nx+2*nb,2),conv_vx(nz+2*nb,nb,2)

  integer::i1,i2,ib,nzpad,nxpad
  real::diff1,diff2,b,tmp

  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667
  real,parameter::c3=+0.009570312500000
  real,parameter::c4=-0.000697544642857

  nzpad=nz+2*nb
  nxpad=nx+2*nb

  !update conv_vz
  do i2=1,nxpad
     do i1=5,nb !top
        b=exp(-bndr(i1)*vv(i1,i2)*dt)
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))&
             +c3*(vz(i1+2,i2)-vz(i1-3,i2))&
             +c4*(vz(i1+3,i2)-vz(i1-4,i2))
        conv_vz(i1,i2,1)=b*conv_vz(i1,i2,1)+(b-1.)*diff1*idz
     enddo
     do i1=nz+nb+1,nzpad-3 !bottom
        ib=nzpad-i1+1
        b=exp(-bndr(ib)*vv(i1,i2)*dt)
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))&
             +c3*(vz(i1+2,i2)-vz(i1-3,i2))&
             +c4*(vz(i1+3,i2)-vz(i1-4,i2))
        conv_vz(ib,i2,2)=b*conv_vz(ib,i2,2)+(b-1.)*diff1*idz
     enddo
  enddo

  !update conv_vx
  do i1=1,nzpad
     do i2=5,nb !left
        b=exp(-bndr(i2)*vv(i1,i2)*dt)
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))&
             +c3*(vx(i1,i2+2)-vx(i1,i2-3))&
             +c4*(vx(i1,i2+3)-vx(i1,i2-4))
        conv_vx(i1,i2,1)=b*conv_vx(i1,i2,1)+(b-1.)*diff2*idx
     enddo
     do i2=nx+nb+1,nxpad-3 !right
        ib=nxpad-i2+1
        b=exp(-bndr(ib)*vv(i1,i2)*dt)
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))&
             +c3*(vx(i1,i2+2)-vx(i1,i2-3))&
             +c4*(vx(i1,i2+3)-vx(i1,i2-4))
        conv_vx(i1,ib,2)=b*conv_vx(i1,ib,2)+(b-1.)*diff2*idx
     enddo
  enddo

  !update pz
  do i2=1,nxpad
     do i1=1,nb !top
        tmp=vv(i1,i2);        tmp=rho(i1,i2)*tmp*tmp
        p(i1,i2)=p(i1,i2)-dt*tmp*conv_vz(i1,i2,1)
     enddo
     do i1=nz+nb+1,nzpad !bottom
        ib=nzpad-i1+1
        tmp=vv(i1,i2);        tmp=rho(i1,i2)*tmp*tmp
        p(i1,i2)=p(i1,i2)-dt*tmp*conv_vz(ib,i2,2)
     enddo
  enddo
  
  !update px
  do i1=1,nzpad
     do i2=1,nb !left
        tmp=vv(i1,i2);        tmp=rho(i1,i2)*tmp*tmp
        p(i1,i2)=p(i1,i2)-dt*tmp*conv_vx(i1,i2,1)
     enddo
     do i2=nx+nb+1,nxpad !right
        ib=nxpad-i2+1
        tmp=vv(i1,i2);        tmp=rho(i1,i2)*tmp*tmp
        p(i1,i2)=p(i1,i2)-dt*tmp*conv_vx(i1,ib,2)
     enddo
  enddo

  return
end subroutine update_cpml_pzpx

!------------------------------------------------------------------------------
subroutine add_attenuation(p, eta, rho, vv, dt, nzpad, nxpad)
  implicit none

  integer::i1,i2
  real::a,tau

  integer::nzpad,nxpad
  real::dt
  real, dimension(nzpad,nxpad)::p, eta, rho, vv

  do i2=1,nxpad
     do i1=1,nzpad
        a=rho(i1,i2)*vv(i1,i2)*vv(i1,i2)
        tau=eta(i1,i2)/a
        a=exp(-dt/tau)
        p(i1,i2)=a*p(i1,i2)
     enddo
  enddo

  return
end subroutine add_attenuation

!-------------------------------------------------------------------------------
subroutine add_sources(p, eta, rho, vv, dt, wlt, sz, sx, nzpad, nxpad)
  implicit none

  integer::sz,sx,nzpad, nxpad
  real::dt,wlt
  real,dimension(nzpad, nxpad)::p, eta, rho, vv

  real::a, tau

  a=rho(sz,sx)*vv(sz,sx)*vv(sz,sx)
  tau=eta(sz,sx)/a
  a=exp(-dt/tau)
  p(sz,sx)=p(sz,sx)+tau*(1.-a)*wlt

  return
end subroutine add_sources

