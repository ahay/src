!  Backward reconstruction based on DFT and inverse DFT
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
!!$
!!$  References: 
!!$  [1] Pengliang Yang, Romain Brossier, Ludovic Metivier and  Jean Virieux,
!!$      Reverse propagation of viscoacoustic forward wavefield with Maxwell 
!!$      attenuation using snapshots and saved boundaries, Technical report 
!!$	 No 83 - SEISCOPE project, University Joseph Fourier
!!$  [2] Pengliang Yang, Romain Brossier,and  Jean Virieux, Boundary reconstruction 
!!$      aftersignificant downsampling, Technical report No 84 - SEISCOPE project
!!$      University Joseph Fourier
program mexwell_dftinterp
  use rsf
  implicit none

  logical :: order1,attenuating
  integer :: ib, it, nt, mt, nz, nx, nb, sx, sz, nxpad, nzpad
  integer :: nsnap, isnap, snap_interval,r
  real  :: dt, dz, dx, fm, tmp, idx, idz, dtb
  real*8, parameter::PI=4.*atan(1.)
  real, dimension (:),   allocatable :: wlt,bb,aa,ef,eb, trf, trb
  real, dimension (:,:), allocatable :: v0, vv, rho, eta
  real, dimension (:,:), allocatable :: p, vz, vx
  real, dimension(:,:,:),allocatable :: conv_pz,conv_px,conv_vz,conv_vx
  complex, dimension(:,:,:,:),allocatable::bvz, bvx
  real, dimension(:,:,:),allocatable:: psnap, vzsnap, vxsnap
  type(file) :: Fv, Fw1, Fw2, Frho, Feta, Fef,Feb,Ftrf,Ftrb ! I/O files 

  call sf_init() ! initialize Madagascar

  ! setup I/O files 
  Fv = rsf_input ("in")   ! velocity model
  Fw1 =rsf_output("out")  ! output forward wavefield 
  Fw2 =rsf_output("back") ! output backward reconstructed wavefield
  Frho=rsf_input("rho")   ! density
  Feta=rsf_input("eta")   ! Pascal
  Fef=rsf_output("ef")    ! energy (kinematic+deformation) when forward modeling
  Feb=rsf_output("eb")    ! energy (kinematic+deformation) when back-reconstruction
  Ftrf=rsf_output("trf")  ! trf
  Ftrb=rsf_output("trb")  ! trb

  ! Read/Write axes
  call from_par(Fv,"n1",nz) ! velocity model: nz
  call from_par(Fv,"n2",nx) ! velocity model: nx
  call from_par(Fv,"d1",dz) ! velocity model: dz
  call from_par(Fv,"d2",dx) ! velocity model: dx

  call from_par("nb", nb, 30) ! thinkness of sponge ABC
  call from_par("nt", nt, 1000) !number of time steps
  call from_par("dt", dt, 0.001) ! time sampling interval
  call from_par("dtb", dtb, dt) ! dt for boundary saving
  call from_par("fm", fm, 20.) ! domainant frequency for ricker wavelet
  call from_par("nsnap", nsnap, 0) ! number of snapshots for re-initialization
  call from_par("order1",order1,.true.) ! 1st order or 2nd order accuracy
  call from_par("attenuating",attenuating,.true.) ! add attenuation or not

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

  call to_par(Ftrf,"n1",nt)
  call to_par(Ftrf,"o1",0)
  call to_par(Ftrf,"d1",dt)
  call to_par(Ftrf,"n2",1)

  idx=1./dx
  idz=1./dz
  nzpad=nz+2*nb
  nxpad=nx+2*nb
  sx=nxpad/2
  sz=nzpad/2
  r=nint(dtb/dt)
  if(mod(nt,r)==0) then
     write(0,*) 'store boundary every m step: m=', r
     mt=nt/r !N=nt;M=mt; M=N/m
  else
     write(0,*) 'dt for boundary saving is chosen inappropriately.'
     call exit(0)
  endif


  allocate(wlt(nt))
  allocate(bb(nb))
  allocate(aa(nb))
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
  allocate(bvz(3,nx,2,mt/2+1))
  allocate(bvx(nz,3,2,mt/2+1))
  allocate(ef(nt))
  allocate(eb(nt))
  if(nsnap>0) then !using snapshot strategy
     allocate(psnap(nzpad,nxpad,nsnap))
     allocate(vzsnap(nzpad,nxpad,nsnap))
     allocate(vxsnap(nzpad,nxpad,nsnap))
     snap_interval=int(nt/(nsnap+1))
  endif
  allocate(trf(nt))
  allocate(trb(nt))

  !generate ricker wavelet with a delay
  do it=1,nt  
     tmp=PI*fm*(it*dt-1.0/fm)
     tmp=tmp*tmp
     wlt(it)=1000.*(1.0-2.0*tmp)*exp(-tmp)
  enddo
  call rsf_read(Fv,v0)
  call check_sanity(maxval(v0),dt,dx,dz)
  call expand2d(vv, v0, nz, nx, nb)
  call rsf_read(Frho,v0)
  call expand2d(rho, v0, nz, nx, nb)
  call rsf_read(Feta, v0)
  call expand2d(eta, v0, nz, nx, nb)
  call cpmlcoeff_init(bb,aa,maxval(vv),dx,dt,fm,nb)
  p=0.
  vx=0.
  vz=0.
  conv_pz=0.
  conv_px=0.
  conv_vz=0.
  conv_vx=0.
  bvx=0.
  bvz=0.
  trf=0.

  !forward modeling
  if(nsnap>0) isnap=0
  do it=1,nt
     if (order1) then ! scheme 1, 1st order accuracy, default
        call step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
        call step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_pzpx(p,vz,vx,conv_vz,conv_vx,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
        if(attenuating) call apply_attenuation(p, eta, rho, vv, dt, nzpad, nxpad)
     else ! scheme 2, 2nd order accuracy
        if(attenuating) call apply_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad)
        call step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
        call step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
        call update_cpml_pzpx(p,vz,vx,conv_vz,conv_vx,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
        if(attenuating) call apply_attenuation(p, eta, rho, vv, 0.5*dt, nzpad, nxpad)
     endif
     call add_sources(p, dt, wlt(it), sz, sx, nzpad, nxpad)

     !boundary saving by computing the Fourier coefficients on the fly
     if(mod(it-1,r)==0) call boundary_interp(.true.,bvz,bvx,vz,vx,nz,nx,nb,it,nt,mt)

     if(nsnap>0) then !save snapshots for re-initialization in backward reconstruction
        if ((isnap<nsnap).and.(it==(isnap+1)*snap_interval)) then
           isnap=isnap+1
           !write(0,*)"isnap=",isnap
           psnap(:,:,isnap)=p
           vzsnap(:,:,isnap)=vz
           vxsnap(:,:,isnap)=vx       
        endif
     endif
     
     trf(it)=p(nb+nz/2,nb+nx+2)

     call window2d(v0, p, nz, nx, nb);
     call rsf_write(Fw1,v0)
     call compute_energy(ef(it),vz(nb+1:nb+nz,nb+1:nb+nx),&
          vx(nb+1:nb+nz,nb+1:nb+nx),p(nb+1:nb+nz,nb+1:nb+nx),&
          rho(nb+1:nb+nz,nb+1:nb+nx),vv(nb+1:nb+nz,nb+1:nb+nx),nz,nx)
  enddo
  call rsf_write(Fef,ef) ! store forward energy history
  call rsf_write(Ftrf,trf)

  !backward reconstruction
  do it=nt,1,-1
     !reconstruct the boundary wavefield by Fourier interpolation
     call boundary_interp(.false.,bvz,bvx,vz,vx,nz,nx,nb,it,nt,mt)

     trb(it)=p(nb+nz/2,nb+nx-2)

     call compute_energy(eb(it),vz(nb+1:nb+nz,nb+1:nb+nx),&
          vx(nb+1:nb+nz,nb+1:nb+nx),p(nb+1:nb+nz,nb+1:nb+nx),&
          rho(nb+1:nb+nz,nb+1:nb+nx),vv(nb+1:nb+nz,nb+1:nb+nx),nz,nx)

     call window2d(v0, p, nz, nx, nb)
     call rsf_write(Fw2,v0)
     if(nsnap>0) then !save snapshots for re-initialization in backward reconstruction
        if((isnap>0).and.(it==isnap*snap_interval)) then
           p=psnap(:,:,isnap)
           vz=vzsnap(:,:,isnap)
           vx=vxsnap(:,:,isnap)
           !write(0,*)"isnap=",isnap
           isnap=isnap-1
        endif
     endif
     call add_sources(p, -dt, wlt(it), sz, sx, nzpad, nxpad)
     if (order1) then ! scheme 1, 1st order accuracy, default
        if(attenuating) call apply_attenuation(p, eta, rho, vv, -dt, nzpad, nxpad)
        call step_forward_p(p, vz, vx, vv, rho, -dt, idz, idx, nzpad, nxpad)
        call step_forward_v(p, vz, vx, vv, rho, -dt, idz, idx, nzpad, nxpad)
     else ! scheme 2, 2nd order accuracy
        if(attenuating) call apply_attenuation(p, eta, rho, vv, -0.5*dt, nzpad, nxpad)
        call step_forward_p(p, vz, vx, vv, rho, -dt, idz, idx, nzpad, nxpad)
        call step_forward_v(p, vz, vx, vv, rho, -dt, idz, idx, nzpad, nxpad)
        if(attenuating) call apply_attenuation(p, eta, rho, vv,-0.5*dt, nzpad, nxpad)
     endif
  enddo
  call rsf_write(Feb,eb) !store backward energy history
  call rsf_write(Ftrb,trb)

  deallocate(wlt)
  deallocate(bb)
  deallocate(aa)
  deallocate(ef)
  deallocate(eb)
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
  deallocate(bvz)
  deallocate(bvx)
  if(nsnap>0) then
     deallocate(psnap)
     deallocate(vzsnap)
     deallocate(vxsnap)
  endif
  deallocate(trf)
  deallocate(trb)

  call exit(0)
end program mexwell_dftinterp

!------------------------------------------------------------------------------
! check the CFL/stability condition is satisfied or not
subroutine check_sanity(vpmax,dt,dx,dz)
  implicit none

  real::vpmax,dt,dx,dz
  real,parameter::c1=+1.196289062500000
  real,parameter::c2=-0.079752604166667

  real CFL,tmp

  tmp=c1-c2
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
end subroutine window2d

!-------------------------------------------------------------------------------
subroutine step_forward_v(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2

  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  do i2=2,nxpad-2
     do i1=2,nzpad-2
        diff1=c1*(p(i1+1,i2)-p(i1,i2))&
             +c2*(p(i1+2,i2)-p(i1-1,i2))
        diff2=c1*(p(i1,i2+1)-p(i1,i2))&
             +c2*(p(i1,i2+2)-p(i1,i2-1))
        vz(i1,i2)=vz(i1,i2)-dt*idz*diff1/rho(i1,i2)
        vx(i1,i2)=vx(i1,i2)-dt*idx*diff2/rho(i1,i2)
     enddo
  enddo
end subroutine step_forward_v

!------------------------------------------------------------------------------
subroutine step_forward_p(p, vz, vx, vv, rho, dt, idz, idx, nzpad, nxpad)
  implicit none

  integer::i1, i2
  real::tmp,diff1,diff2

  integer:: nzpad, nxpad
  real::idz,idx,dt
  real,dimension(nzpad,nxpad)::p, vz, vx, vv, rho

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  do i2=3,nxpad-1
     do i1=3,nzpad-1
        tmp=vv(i1,i2)
        tmp=rho(i1,i2)*tmp*tmp
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))
        p(i1,i2)=p(i1,i2)-dt*tmp*(idz*diff1+idx*diff2)
     enddo
  enddo
end subroutine step_forward_p

!------------------------------------------------------------------------------
! initialize PML damping profile
subroutine cpmlcoeff_init(bb,aa,vpmax,dx,dt,fm,nb)
  implicit none

  integer::nb
  real::dt,fm,dx,vpmax
  real,dimension(nb)::bb,aa

  integer::ib
  real::x,L,d0,d,alpha_max,alpha,r
  real,parameter::Rc=1.e-4

  L=nb*dx
  d0=-3.*log(Rc)*vpmax/(2.*L)
  alpha_max=4.*atan(1.)*fm

  do ib=1,nb
     x=(nb-ib+1)*dx     !x=1.-cos(0.5*(nb-ib)*PI/nb)
     r=x/L ! ratio of x over L
     d=d0*r**2; alpha=alpha_max*r
     bb(ib)=exp(-(d+alpha)*dt)
     aa(ib)=d*(bb(ib)-1.)/(d+alpha)
  enddo
end subroutine cpmlcoeff_init

!------------------------------------------------------------------------------
subroutine update_cpml_vzvx(p,vz,vx,conv_pz,conv_px,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
  implicit none

  integer::nz,nx,nb
  real::idz,idx,dt
  real,dimension(nb)::bb,aa
  real,dimension(nz+2*nb,nx+2*nb)::p,vz,vx,rho,vv
  real::conv_pz(nb,nx+2*nb,2),conv_px(nz+2*nb,nb,2)

  integer::nzpad,nxpad,i1,i2,ib
  real*8::b,diff1,diff2

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  !update conv_pz
  do i2=1,nxpad
     do i1=2,nb !top
        diff1=c1*(p(i1+1,i2)-p(i1,i2)) &
             +c2*(p(i1+2,i2)-p(i1-1,i2))
        conv_pz(i1,i2,1)=bb(i1)*conv_pz(i1,i2,1)+aa(i1)*diff1*idz
     enddo
     do i1=nz+nb+1,nzpad-2 !bottom
        ib=nzpad-i1+1
        diff1=c1*(p(i1+1,i2)-p(i1,i2)) &
             +c2*(p(i1+2,i2)-p(i1-1,i2))
        conv_pz(ib,i2,2)=bb(ib)*conv_pz(ib,i2,2)+aa(ib)*diff1*idz
     enddo
  enddo

  !update conv_px
  do i1=1,nzpad
     do i2=2,nb !left
        diff2=c1*(p(i1,i2+1)-p(i1,i2))&
             +c2*(p(i1,i2+2)-p(i1,i2-1)) 
        conv_px(i1,i2,1)=bb(i2)*conv_px(i1,i2,1)+aa(i2)*diff2*idx
     enddo
     do i2=nx+nb+1,nxpad-2 !right
        ib=nxpad-i2+1
        diff2=c1*(p(i1,i2+1)-p(i1,i2)) &
             +c2*(p(i1,i2+2)-p(i1,i2-1))
        conv_px(i1,ib,2)=bb(ib)*conv_px(i1,ib,2)+aa(ib)*diff2*idx
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
        vx(i1,i2)=vx(i1,i2)-dt*conv_px(i1,ib,2)/rho(i1,i2)
     enddo
  enddo
end subroutine update_cpml_vzvx

!------------------------------------------------------------------------------
subroutine update_cpml_pzpx(p,vz,vx,conv_vz,conv_vx,rho,vv,bb,aa,idz,idx,dt,nz,nx,nb)
  implicit none

  integer::nz,nx,nb
  real::idz,idx,dt
  real,dimension(nb)::bb,aa
  real,dimension(nz+2*nb,nx+2*nb)::p,vz,vx,rho,vv
  real::conv_vz(nb,nx+2*nb,2),conv_vx(nz+2*nb,nb,2)

  integer::i1,i2,ib,nzpad,nxpad
  real*8::diff1,diff2,b,tmp

  real,parameter::c1=1.125
  real,parameter::c2=-1./24.

  nzpad=nz+2*nb
  nxpad=nx+2*nb

  !update conv_vz
  do i2=1,nxpad
     do i1=3,nb !top
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))
        conv_vz(i1,i2,1)=bb(i1)*conv_vz(i1,i2,1)+aa(i1)*diff1*idz
     enddo
     do i1=nz+nb+1,nzpad-1 !bottom
        ib=nzpad-i1+1
        diff1=c1*(vz(i1,i2)-vz(i1-1,i2))&
             +c2*(vz(i1+1,i2)-vz(i1-2,i2))
        conv_vz(ib,i2,2)=bb(ib)*conv_vz(ib,i2,2)+aa(ib)*diff1*idz
     enddo
  enddo

  !update conv_vx
  do i1=1,nzpad
     do i2=3,nb !left
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))
        conv_vx(i1,i2,1)=bb(i2)*conv_vx(i1,i2,1)+aa(i2)*diff2*idx
     enddo
     do i2=nx+nb+1,nxpad-1 !right
        ib=nxpad-i2+1
        diff2=c1*(vx(i1,i2)-vx(i1,i2-1))&
             +c2*(vx(i1,i2+1)-vx(i1,i2-2))
        conv_vx(i1,ib,2)=bb(ib)*conv_vx(i1,ib,2)+aa(ib)*diff2*idx
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
end subroutine update_cpml_pzpx

!------------------------------------------------------------------------------
subroutine apply_attenuation(p, eta, rho, vv, dt, nzpad, nxpad)
  implicit none

  integer::i1,i2
  real*8::a,tau

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
end subroutine apply_attenuation

!-------------------------------------------------------------------------------
subroutine add_sources(p, dt, wlt, sz, sx, nzpad, nxpad)
  implicit none

  integer::sz,sx,nzpad, nxpad
  real::dt,wlt
  real,dimension(nzpad, nxpad)::p

  p(sz,sx)=p(sz,sx)+dt*wlt
end subroutine add_sources

!-------------------------------------------------------------------------------
subroutine compute_energy(e,vz,vx,p,rho,vv,nz,nx)
  !compute the sum of kinematic energy and deformation/potential energy
  implicit none

  integer::nz,nx
  real::e
  real,dimension(nz,nx)::rho,vv,p,vz,vx

  e=0.5*sum((vz*vz+vx*vx)*rho+p*p/(rho*vv*vv))
end subroutine compute_energy

!-------------------------------------------------------------------------------
subroutine boundary_interp(v2b,bvz,bvx,vz,vx,nz,nx,nb,it,nt,mt)
  !interpolate the boundary elements by DFT and its inverse
  implicit none

  logical::v2b !v to bourndary or reverse
  integer::nz,nx,nb,mt,nt,it
  real,dimension(nz+2*nb,nx+2*nb)::vz,vx
  complex::bvz(3,nx,2,mt/2+1),bvx(nz,3,2,mt/2+1)
  real*8::pi=4.*atan(1.)

  ! n=it; n'=itt;
  ! N=nt; N'=mt;
  integer::i1,i2,ik,r
  complex::factor(mt/2+1)
  r=nt/mt

  if(v2b) then !v to bourndary
     do ik=1,mt/2+1
        !complex exponential factor
        factor(ik)=exp(cmplx(0,-2.*pi*(ik-1)*(it-1)/nt))
     enddo
     do i2=1,nx
        do i1=1,3
           bvz(i1,i2,1,:)=bvz(i1,i2,1,:)+vz(i1+nb-2,i2+nb)*factor
           bvz(i1,i2,2,:)=bvz(i1,i2,2,:)+vz(i1+nz+nb-2,i2+nb)*factor
        enddo
     enddo
     do i1=1,nz
        do i2=1,3
           bvx(i1,i2,1,:)=bvx(i1,i2,1,:)+vx(i1+nb,i2+nb-2)*factor
           bvx(i1,i2,2,:)=bvx(i1,i2,2,:)+vx(i1+nb,i2+nx+nb-2)*factor
        enddo
     enddo
  else !boundary to v
     do ik=1,mt/2+1
        factor(ik)=exp(cmplx(0,2.*pi*(ik-1)*(it-1)/nt))/mt
     enddo
     factor(1)=0.5*factor(1)
     !handle the case mt is an even number
     if(mod(mt,2)==0) factor(mt/2+1)=0.5*factor(mt/2+1)
     do i2=1,nx
        do i1=1,3
           vz(i1+nb-2,i2+nb)=2.*real(sum(bvz(i1,i2,1,:)*factor))
           vz(i1+nz+nb-2,i2+nb)=2.*real(sum(bvz(i1,i2,2,:)*factor))
        enddo
     enddo
     do i1=1,nz
        do i2=1,3
           vx(i1+nb,i2+nb-2)=2.*real(sum(bvx(i1,i2,1,:)*factor))
           vx(i1+nb,i2+nx+nb-2)=2.*real(sum(bvx(i1,i2,2,:)*factor))
        enddo
     enddo
  endif
end subroutine boundary_interp
