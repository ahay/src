!Acoustic W.E 2d 4-2 order FDTD modeling, regular grids

!  Copyright (C) 2015 ZHENDONG Zhang @ KAUST  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


program fdac2d
use rsf
implicit none



type(file)                  :: Fslo, Fglo,Fwav, Fvel, Frec, Fwfd            !I/O files
type(axa)                   :: at, az, ax, asx, asz, agz, agx       !cube axes

logical                     :: verb
integer, allocatable        :: sz(:), sx(:), gz(:), gx(:)
real*4,  allocatable        :: slo(:,:), glo(:,:)
real*4,  allocatable        :: velmod(:,:), vel(:,:), wavelet(:)
real*4,  allocatable        :: seisout(:,:)


integer                     ::  nz, nx, nt, ns, ng, npml, jsnap
integer                     ::  iz, ix, it, is, ig
real*4                      ::  dz, dx, dt, c
real*4                      ::  d, d1



!for OpenMp choosing the cores
integer                     ::  Ncore
!*******************allocatable dynamic arraies**********************************
real*4,allocatable          ::  u_plus(:,:), u(:,:), u_minus(:,:)
real*4,allocatable          ::  p_plus(:,:), p(:,:), p_minus(:,:)
real*4,allocatable          ::  q_plus(:,:), q(:,:)
real*4,allocatable          ::  temp_plus(:,:), temp(:,:), temp_minus(:,:)
real*4,allocatable          ::  t_plus(:,:), t(:,:), t_minus(:,:)

call sf_init()
call from_par("verb", verb, .true.)


Fslo  = rsf_input("slo")
Fglo  = rsf_input("glo")
Fvel  = rsf_input("vel")
Fwav  = rsf_input("wavelet")

Frec  = rsf_output("frec")


! Read/Write axes
call iaxa(Fwav, at,  1)
call iaxa(Fvel, az,  1); call iaxa(Fvel, ax,  2)
call iaxa(Fslo, asz, 1); call iaxa(Fslo, asx, 2)
call iaxa(Fglo, agz, 1); call iaxa(Fglo, agx, 2)
call oaxa(Frec, at,  1); call oaxa(Frec, agz, 2)


nz     =  az%n
nx     =  ax%n
ns     =  asz%n
ng     =  agz%n
dx     =  ax%d
dz     =  az%d
dt     =  at%d
nt     =  at%n


call from_par("npml",     npml)          ! Grid points for PML
call from_par("nthreads", ncore)         ! Threads for OMP

allocate(slo(ns,2)); slo=0.0; call rsf_read(Fslo, slo)
allocate(glo(ng,2)); glo=0.0; call rsf_read(Fglo, glo)


!read velocity model
allocate(velmod(nz,nx)); velmod  = 0.0 ; call rsf_read(Fvel, velmod)
allocate(wavelet(nt))  ; wavelet = 0.0 ; call rsf_read(Fwav, wavelet)


 !********************************allocate******************************

!write(*,*)"allo start"
allocate(u_plus(nz+2*npml,nx+2*npml), u(nz+2*npml,nx+2*npml), u_minus(nz+2*npml,nx+2*npml))
allocate(p_plus(nz+2*npml,nx+2*npml), p(nz+2*npml,nx+2*npml), p_minus(nz+2*npml,nx+2*npml))
allocate(q_plus(nz+2*npml,nx+2*npml), q(nz+2*npml,nx+2*npml))
allocate(temp_plus(nz+2*npml,nx+2*npml), temp(nz+2*npml,nx+2*npml), temp_minus(nz+2*npml,nx+2*npml))
allocate(t_plus(nz+2*npml,nx+2*npml), t(nz+2*npml,nx+2*npml), t_minus(nz+2*npml,nx+2*npml))

allocate( vel(nz+2*npml,nx+2*npml), seisout(nt,ng))
allocate(sz(ns), sx(ns), gz(ng), gx(ng))

c          = 1/sqrt(real(2.0))

u          = 0.0
u_minus    = 0.0
u_plus     = 0.0

p          = 0.0
p_minus    = 0.0
p_plus     = 0.0

q          = 0.0
q_plus     = 0.0

temp_minus = 0.0
temp_plus  = 0.0
temp       = 0.0

t          = 0.0
t_minus    = 0.0
t_plus     = 0.0
vel        = 0.0
seisout    = 0.0


!******************initialize**************************


!----------------------------------------------Expand the velocity model-----------------------------------------------------------
do iz = 1, nz
   do ix = 1, nx
      vel(npml+iz, npml+ix) = velmod(iz,ix)
   end do
end do

do iz = 1,npml                      !top
   do ix = 1,nx
      vel(iz,npml+ix) = velmod(1,ix)
   end do
end do

do ix = 1,npml
   vel(1:ix,ix) = vel(1:ix,npml+1)
end do
do ix = nx+npml+1,nx+2*npml
   vel(1:(nx+2*npml-ix+1),ix) = vel(1:(nx+2*npml-ix+1),nx+npml)
end do


do iz = npml+nz+1,nz+2*npml                      !bottom
   do ix = 1,nx
      vel(iz,npml+ix) = velmod(nz,ix)
   end do
end do
do ix = 1,npml
   vel(nz+2*npml-ix+1:nz+2*npml,ix) = vel(nz+2*npml-ix+1:nz+2*npml,npml+1)
end do
do ix = nx+npml+1,nx+2*npml
   vel(nz+(ix-nx):nz+2*npml,ix) = vel(nz+(ix-nx):nz+2*npml,nx+npml)
end do


do ix = 1,npml                      !left
   do iz = 1,nz
      vel(npml+iz,ix) = velmod(iz,1)
   end do
end do
do iz = 1,npml
   vel(iz,1:iz) = vel(npml+1,1:iz)
end do
do iz = nz+npml+1,nz+2*npml
   vel(iz, 1:(nz+2*npml-iz+1)) = vel(npml+nz, 1:(nz+2*npml-iz+1))
end do


do ix = npml+nx+1,nx+2*npml                      !right
   do iz = 1,nz
      vel(npml+iz,ix) = velmod(iz,nx)
   end do
end do
do iz = 1,npml
   vel(iz,nx+2*npml-iz+1:nx+2*npml) = vel(npml+1,nx+2*npml-iz+1:nx+2*npml)
end do
do iz = nz+npml+1,nz+2*npml
   vel(iz, nx+npml+(iz-nz-npml):nx+2*npml) = vel(npml+nz, nx+npml+(iz-nz-npml):nx+2*npml)
end do


!write(0,*) vel(1:npml,1:npml), vel(1:npml,nx+npml+1:2*npml+nx),vel(nz+npml+1:2*npml+nz,1:npml)
!write(0,*) vel(nz+npml+1:2*npml+nz,nx+npml+1:2*npml+nx)
!----------------------------------------------------------------------------------------------------------------------------------
!open (97,file='vel.bin',status='replace',form='unformatted',access='direct',recl=(nz+2*npml)*(nx+2*npml)*4)
!write(97,rec=1) vel
!close(97)


do is=1,ns
   sz(is) = nint(slo(is,2) / (1.0*dx)) +npml
   sx(is) = nint(slo(is,1) / (1.0*dx)) +npml
end do
do ig=1,ng
   gz(ig) = nint(glo(ig,2) / (1.0*dx)) +npml
   gx(ig) = nint(glo(ig,1) / (1.0*dx)) +npml
end do

write(0,*) 'Start...'


u=0.0
u_minus=0.0
u_plus=0.0
p=0.0
p_minus=0.0
p_plus =0.0
q=0.0
q_plus=0.0

temp_minus=0.0
temp_plus=0.0
temp=0.0

t=0.0
t_minus =  0.0
t_plus =  0.0




do it=1,nt

!$OMP PARALLEL DO NUM_THREADS(Ncore), PRIVATE(iz,ix), SCHEDULE(DYNAMIC)

   do ix=1+npml,nx+npml

      do iz=1+npml,nz+npml

         u_plus(iz,ix)=2*u(iz,ix)-u_minus(iz,ix)+(vel(iz,ix)**2*dt**2/dx**2)*(u(iz,ix+1)&
              &-2*u(iz,ix)+u(iz,ix-1)+u(iz+1,ix)-2*u(iz,ix)+u(iz-1,ix))

	  !  u_plus(iz,ix)=2*u(iz,ix)-u_minus(iz,ix)+(vel(iz,ix)**2*dt**2/(12*dx**2))*(-u(iz,ix+2)+16*u(iz,ix+1)-&
	  !  &30*u(iz,ix)+16*u(iz,ix-1)-u(iz,ix-2)-u(iz+2,ix)+16*u(iz+1,ix)-30*u(iz,ix)+16*u(iz-1,ix)-u(iz-2,ix))

      end do
   end do


!$OMP END PARALLEL DO

   do is = 1,ns
     
      u(sz(is), sx(is))=u(sz(is),sx(is))-vel(sz(is),sx(is))**2*dt**2*wavelet(it)*1000.0
      
   end do




!TOP

  do ix=npml+1,nx+npml-1
     do iz=npml+1,2,-1
        d=9*(iz-npml-2)**2*vel(iz,ix)/(npml**3*dz)
        d1=18*(iz-npml-2)*vel(iz,ix)/(npml**3*dz**2)
       
     
        p_plus(iz,ix)=((U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
             &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
        
        temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
             &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
        q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
		
        t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
        u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
     
     end do
  end do


  do ix=2,npml
     do iz=npml,2,-1
        d=8*((iz-npml-2)*c)**2*vel(iz,ix)/(npml**3*dz)
        d1=16*c*((iz-npml-2)*c)*vel(iz,ix)/(npml**3*dz**2)
        p_plus(iz,ix)=((U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
             &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
		
        temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
             &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
        q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
				   	   
        t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
        u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)

     end do
  end do

  do ix=nx+npml,nx+2*npml-2
     do iz=npml,2,-1
        d=8*((iz-npml-2)*c)**2*vel(iz,ix)/(npml**3*dz)
        d1=16*c*((iz-npml-2)*c)*vel(iz,ix)/(npml**3*dz**2)
        p_plus(iz,ix)=((U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
             &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
		
        temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
             &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
        q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)		   
	   
        t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
        u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
   
     end do
 end do


!TOP end


!BOTTOM

   do ix=npml+1,nx+npml-1
      do iz=nz+npml-1,nz+2*npml-2
         d=9*(iz-nz-npml+2)**2*vel(iz,ix)/(npml**3*dz)
         d1=18*(iz-nz-npml+2)*vel(iz,ix)/(npml**3*dz**2)

         p_plus(iz,ix)=((U(iz+1,ix)-2*U(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
			 
         temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
         q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
             	      
         t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*U(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
	   
         u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
      end do
   end do	 

    

   do ix=2,npml	
      do iz=nz+npml-1,nz+2*npml-2
         d=8*((iz-nz-npml+2)*c)**2*vel(iz,ix)/(npml**3*dz)
         d1=16*c*((iz-nz-npml+2)*c)*vel(iz,ix)/(npml**3*dz**2)
   
         p_plus(iz,ix)=((U(iz+1,ix)-2*U(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
              &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
			
         temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
         q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
	   
         t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*U(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
	   
         u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)

      end do
   end do



   do ix=nx+npml,nx+2*npml-2	
      do iz=nz+npml-1,nz+2*npml-2
         d=8*((iz-nz-npml+2)*c)**2*vel(iz,ix)/(npml**3*dz)
         d1=16*c*((iz-nz-npml+2)*c)*vel(iz,ix)/(npml**3*dz**2)
  
         p_plus(iz,ix)=((U(iz+1,ix)-2*U(iz,ix)+U(iz-1,ix))*vel(iz,ix)**2/dz**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
			
         temp_plus(iz,ix)=(-0.5*d1*(U(iz+1,ix)-U(iz-1,ix))*vel(iz,ix)**2/dz-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
         q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
			
         t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz,ix+1)-2*U(iz,ix)+U(iz,ix-1))/dx**2+2*t(iz,ix)-t_minus(iz,ix)
	   
         u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)


      end do
  end do



!BOTTOM END





!LEFT


 do ix=npml+1,2,-1	
    d=9*(ix-npml-2)**2*vel(iz,ix)/(npml**3*dz)
    d1=18*(ix-npml-2)*vel(iz,ix)/(npml**3*dz**2)
    do iz=npml+1,nz+npml-1
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
	
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
	
       t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
   
    end do
 end do

 do ix=npml+1,2,-1	
    d=8*((ix-npml-2)*c)**2*vel(iz,ix)/(npml**3*dx)
    d1=16*c*((ix-npml-2)*c)*vel(iz,ix)/(npml**3*dx**2)

    do iz=2,npml
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
       
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
		
       T_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)

    end do
 end do



 do ix=npml+1,2,-1	
    d=8*((ix-npml-2)*c)**2*vel(iz,ix)/(npml**3*dx)
    d1=16*c*((ix-npml-2)*c)*vel(iz,ix)/(npml**3*dx**2)
    do iz=nz+npml,nz+2*npml-2
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
       
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
		   	   
       t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)

    end do
 end do


!LEFT end


!RIGHT


 do ix=nx+npml-1,nx+2*npml-2	
    d=9*(ix-nx-npml+2)**2*vel(iz,ix)/(npml**3*dx)
    d1=18*(ix-nx-npml+2)*vel(iz,ix)/(npml**3*dx**2)
    do iz=npml+1,nz+npml-1 
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
        
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)
 
       t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
    end do
 end do

 do ix=nx+npml-1,nx+2*npml-2	
    d=8*((ix-nx-npml+2)*c)**2*vel(iz,ix)/(npml**3*dx)
    d1=16*c*((ix-nx-npml+2)*c)*vel(iz,ix)/(npml**3*dx**2)

    do iz=2,npml 
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
    
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)	   
 	    
       t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)
     end do
 end do

 do ix=nx+npml-1,nx+2*npml-2	
    d=8*((ix-nx-npml+2)*c)**2*vel(iz,ix)/(npml**3*dx)
    d1=16*c*((ix-nx-npml+2)*c)*vel(iz,ix)/(npml**3*dx**2)
    do iz=nz+npml,nz+2*npml-2 
       p_plus(iz,ix)=((U(iz,ix+1)-2*u(iz,ix)+U(iz,ix-1))*vel(iz,ix)**2/dx**2-(-2*p(iz,ix)&
            &+p_minus(iz,ix))/dt**2+d*p_minus(iz,ix)/dt-d**2*p(iz,ix))/(1/dt**2+d/dt)
       
       temp_plus(iz,ix)=(-0.5*d1*(U(iz,ix+1)-U(iz,ix-1))*vel(iz,ix)**2/dx-(-2*temp(iz,ix)&
            &+temp_minus(iz,ix))/dt**2+d*temp_minus(iz,ix)/dt-d**2*temp(iz,ix))/(1/dt**2+d/dt)
       q_plus(iz,ix)=dt*temp(iz,ix)-dt*d*q(iz,ix)+q(iz,ix)

       t_plus(iz,ix)=vel(iz,ix)**2*dt**2*(U(iz+1,ix)-2*u(iz,ix)+U(iz-1,ix))/dz**2+2*t(iz,ix)-t_minus(iz,ix)
       u_plus(iz,ix)=p_plus(iz,ix)+q_plus(iz,ix)+t_plus(iz,ix)

    end do
  end do

  u_minus=u
  u=u_plus
 
  p_minus=p
  p=p_plus
 
  temp_minus=temp
  temp=temp_plus


  q=q_plus

  t_minus=t
  t=t_plus
!
  do ig = 1, ng
   
    seisout(it, ig) = u(gz(ig), gx(ig))
   
  end do


!seisforward(:,:, it) = u_plus(npml+1:nz+npml, npml+1:nx+npml)

end do

!end of time loop


!seis = u(npml+1:nz+npml, npml+1:nx+npml)

!print *,seisout

call rsf_write(Frec, seisout)

write(*,*)"Success"

deallocate(u_plus, u, u_minus)
deallocate(p_plus, p, p_minus)
deallocate(q_plus, q)
deallocate(temp_plus, temp, temp_minus)
deallocate(t_plus, t, t_minus)	
deallocate(wavelet,vel, velmod)	 
deallocate(slo,glo, seisout)  
deallocate(sz, sx, gz, gx)
write(0,*) 'End Successfully.'
call exit(0)
end program fdac2d
!subroutine end
