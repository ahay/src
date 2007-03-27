module C2R2D

  implicit none

  integer, private :: nsx,nsz,norm,cnx,cnz,rnx,rnz
  real   , private  :: xmin,xmax,zmin,zmax,dxx,dzz
  real, allocatable, dimension(:,:),private :: fold

contains

  !----------------------------------------------------------------  

  subroutine C2R_init(insx,insz,inorm,idxx,idzz,&
                      icnx,icnz,irnx,irnz,ixmin,izmin)
    integer :: insx,insz,inorm,icnx,icnz,irnx,irnz
    real    :: idxx,idzz,ixmin,izmin
    nsx = insx;  nsz = insz;   norm= inorm
    xmin=ixmin;  zmin= izmin
    cnx = icnx;  cnz = icnz
    rnx = irnx;  rnz = irnz
    dxx = idxx;  dzz = idzz
    write(0,*) 'NSX and NSZ:',nsx,nsz
    allocate(fold(cnz,cnx))

  end subroutine C2R_init

  !----------------------------------------------------------------  

  subroutine C2R_run(RC,rays,CC,gnorm)
    real    :: CC(:,:),RC(:,:),rays(:,:,:),gnorm(:,:)
    integer :: ix,iz,it,ig,im,jx,jz
    real    :: x,z,zz,xx
    real    :: ddz,ddx,r,sincr,dr

    write(0,*) 'MIN/MAX Image to interp',minval(RC),maxval(RC)
    write(0,*) ' RAY DIM',rnz,rnx
    write(0,*) 'CART DIM',cnz,cnx
    write(0,*) 'MIN/MAX Gnorm',minval(gnorm),maxval(gnorm)
    write(0,*) 'NSX and NSZ:',nsx,nsz
    dr = sqrt( dxx**2 + dzz**2 )
    nsx=nsz

    do ig=1,rnx
       do it=1,rnz   
          x=rays(it,ig,1)
          z=rays(it,ig,2)

          iz=floor((z-zmin)/dzz)+1
          ix=floor((x-xmin)/dxx)+1

          if(iz<   1+nsz) cycle
          if(iz> cnz-nsz) cycle
          if(ix<   1+nsx) cycle
          if(ix> cnx-nsx) cycle
          do jz=iz-nsz,iz+nsz
             zz=zmin+(jz-1)*dzz
             ddz=zz-z

             do jx=ix-nsx,ix+nsx
                xx=xmin+(jx-1)*dxx
                ddx=xx-x
                r=sqrt( ddz**2 + ddx**2 )

                if(abs(r)<epsilon(r)) then
                   sincr = 1.0
                else
                   r = r / dr
                   sincr = sin(r)/r
                end if
                CC(jz,jx)=CC(jz,jx)+RC(it,ig)*sincr*(gnorm(it,ig))**2
                fold(jz,jx) = fold(jz,jx) + 1.

             end do
          end do

       end do !! End Tau
    end do !! End Mu
    write(0,*) minval(CC),maxval(CC)
    !CC=CC !/(fold+0.001)
  end subroutine C2R_run

end module C2R2D






