!!$  Copyright (C) 2007 Stanford University
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
module C2R2D
  use rsf
! use nonexistant
  implicit none

  integer, private :: nsx,nsz,norm,cnx,cnz,rnx,rnz
  real   , private  :: xmin,xmax,zmin,zmax,dxx,dzz,maxoff
  real, allocatable, dimension(:,:),private :: fold
  real, allocatable, dimension(:,:),private :: RC,CC,gnorm
  real, allocatable, dimension(:,:,:),private:: rays

contains

  !----------------------------------------------------------------  

 subroutine C2R_init(insx,insz,inorm,idxx,idzz,&
                      icnx,icnz,irnx,irnz,ixmin,izmin)
    integer :: insx,insz,inorm,icnx,icnz,irnx,irnz
    real    :: idxx,idzz,ixmin,izmin
    real    :: maxoffset
    nsx = insx;  nsz = insz;   norm= inorm
    xmin=ixmin;  zmin= izmin
    cnx = icnx;  cnz = icnz
    rnx = irnx;  rnz = irnz
    dxx = idxx;  dzz = idzz
    
 
    allocate(fold(cnz,cnx));    fold=1.
    allocate(CC(cnz,cnx),RC(rnz,rnx),gnorm(rnz,rnx))
    allocate(rays(rnz,rnx,2))
  end subroutine C2R_init

  !---------------------------------------------------------------- 

  subroutine C2R_run(RC,rays,CC,gnorm,adj)
    real    :: CC(:,:),RC(:,:),rays(:,:,:),gnorm(:,:)
    integer :: ix,iz,it,ir,im,jx,jz
    real    :: x,z,zz,xx,RCmax
    real    :: ddz,ddx,r,sincr,dr
    logical :: adj

    write(0,*) 'CART: nx,dx,ox' ,cnx,dxx,xmin
    write(0,*) 'CART: nz,dz,oz' ,cnz,dzz,zmin
    write(0,*) 'RAY : nr,minx,maxx',rnx ,minval(rays(:,:,1)),maxval(rays(:,:,1))
    write(0,*) 'RAY : nt,minz,maxz',rnz ,minval(rays(:,:,2)),maxval(rays(:,:,2))

    gnorm = gnorm / maxval(gnorm)
    dr = sqrt( dxx**2 + dzz**2 )
    do ir=1,rnx 
       do it=1,rnz   
         x=rays(it,ir,1)
          z=rays(it,ir,2)

          ix=floor((x-xmin)/dxx)+1
          if(ix<   1+nsx) cycle
          if(ix> cnx-nsx) cycle

          iz=floor((z-zmin)/dzz)+1
          if(iz<   1+nsz) cycle
          if(iz> cnz-nsz) cycle


          do jx=ix-nsx,ix+nsx
             xx=xmin+(jx-1)*dxx
             ddx=xx-x
             do jz=iz-nsz,iz+nsz
                zz=zmin+(jz-1)*dzz
                ddz=zz-z
                
                r=sqrt( ddz**2 + ddx**2 )

                if(abs(r)<epsilon(r)) then
                   sincr = 1.0
                else
                   r = r / dr
                   sincr = sin(r)/r
                end if

                if (adj) then
                   CC(jz,jx)=CC(jz,jx)+RC(it,ir)*sincr*gnorm(it,ir)
                   fold(jz,jx) = fold(jz,jx) + 1.
                else
                   RC(it,ir) = RC(it,ir) + CC(jz,jx) !*sincr
                end if

             end do
          end do

       end do !! End Tau
    end do !! End Mu

    if (adj) then
       do jx=1,cnx
          do jz=1,cnz
             CC(jz,jx)=CC(jz,jx)/(fold(jz,jx)+0.001)
          end do
       end do
    end if

    if (.not. adj) then
       RCmax = maxval(RC)
       where (RC .eq. 0.)
          RC = RCmax
       end where
    end if

    write(0,*) 'CC: min/max',minval(CC),maxval(CC)
    write(0,*) 'RC: min/max',minval(RC),maxval(RC)

  end subroutine C2R_run

  subroutine C2R_close()
    deallocate( CC,RC,gnorm,rays,fold ) 
  end subroutine C2R_close


end module C2R2D






