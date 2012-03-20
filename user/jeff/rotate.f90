! Module to do 3-component rotation [Free-surface rotation matrix Kennett, 198x]
module rotate 
  use rsf

  implicit none
  integer, private              :: n1,n2,n3
  real   , private              :: alpha,beta,p1,p2,irot

contains

  !! Initialization subroutine
  subroutine rotatedata_init(n1_in,n2_in,n3_in,a_in,b_in,p1_in,p2_in,irot_in)
    integer                      :: n1_in,n2_in,n3_in
    real                         :: a_in,b_in,p1_in,p2_in,irot_in
    n1=n1_in; n2=n2_in; n3=n3_in; alpha=a_in; beta=b_in; p1=p1_in; p2=p2_in; irot=irot_in
  end subroutine rotatedata_init

  !! . . Rotation
  subroutine rotatedata(adj,add,data,modl)
    !! . . Declarations
    integer                               :: ii,jj
    real                                  :: Frp,Fzp,Frsv,Fzsv,qa,qb,delta,pi
    real                                  :: pmag,sina,cosa,slow,invfact,cosi,sini
    real, dimension(n1,n2,n3)             :: data,modl
    real, dimension(:,:,:),allocatable    :: tmp, tmp2
    logical                               :: adj,add

    allocate( tmp(n1,n2,n3), tmp2(n1,n2,n3) )
    if (p1/abs(p1)*p2/abs(p2)==1.0) then
       pmag=sqrt(p1**2+p2**2)   !! . . Magnitudes and Angles
    else
       pmag=-sqrt(p1**2+p2**2)  !! . . Magnitudes and Angles
    end if
    pi = acos(-1.)
    sini=sin(irot * (pi/180.))
    cosi=cos(irot * (pi/180.))
    sina=sin(atan2(p2,p1))
    cosa=cos(atan2(p2,p1))
    slow=pmag                !! . . Constant Definitions
    qa=sqrt(1./(alpha*alpha) - slow**2)
    qb=sqrt(1./(beta * beta) - slow**2)
    delta = 1.-4.*slow**2*beta**2 + 4.*slow**4*beta**4 + 4*beta**4*slow**2*qa*qb
    Frp = 4*alpha*beta**2*slow*qa*qb       /delta
    Fzp =-2*alpha*qa*(1-2.*slow**2*beta**2)/delta
    Frsv= 2*beta *qb*(1-2.*slow**2*beta**2)/delta
    Fzsv= 4.*beta**3*slow*qa*qb            /delta
    invfact=1./(-Frsv*Fzp+Frp*Fzsv)

    if (adj) then   !! . . Model to Data

       write(0,*) 'APPLYING ADJOINT OPERATOR'
      !! . . Data to Model
       do jj=1,n2
          do ii=1,n1            !! . . Rotation to go from [UN,UE,UD] to [u1,u2,u3]
             tmp(ii,jj,1)=  cosi*data(ii,jj,1) + sini*data(ii,jj,2)
             tmp(ii,jj,2)= -sini*data(ii,jj,1) + cosi*data(ii,jj,2)
             tmp(ii,jj,3)=       data(ii,jj,3)
          end do
       end do
       do jj=1,n2
          do ii=1,n1            !! . . Rotation to go from [u1,u2,u3] to [Ur,Ut,Uz]
             tmp2(ii,jj,1)=  cosa*tmp(ii,jj,1) + sina*tmp(ii,jj,2)
             tmp2(ii,jj,2)= -sina*tmp(ii,jj,1) + cosa*tmp(ii,jj,2)
             tmp2(ii,jj,3)=       tmp(ii,jj,3)
          end do
       end do
       do jj=1,n2       
          do ii=1,n1            !! . . From [Ur,Ut,Uz] to [SV,SH,P]
             modl(ii,jj,1)= -Fzp/invfact*tmp2(ii,jj,1) +  Frp/invfact*tmp2(ii,jj,3)
             modl(ii,jj,2)=        - 0.5*tmp2(ii,jj,2)
             modl(ii,jj,3)= Fzsv/invfact*tmp2(ii,jj,1) - Frsv/invfact*tmp2(ii,jj,3)
          end do
       end do
      
    else
  
      write(0,*) 'APPLYING FORWARD OPERATOR'
       do ii=1,n1            !! . . Decompose section from [SV,SH,P] to [R,T,Z]
          do jj=1,n2 
             tmp(ii,jj,1)= Frsv*modl(ii,jj,1) + Fzsv*modl(ii,jj,3)
             tmp(ii,jj,2)= -2.*  modl(ii,jj,2)
             tmp(ii,jj,3)= Frp *modl(ii,jj,1) + Fzp *modl(ii,jj,3)
          end do
       end do

       do ii=1,n1            !! . . Rotate Sections from [R,T,Z] to [u1,u2,u3]
          do jj=1,n2
             tmp2(ii,jj,1)= cosa*tmp(ii,jj,1) - sina*tmp(ii,jj,2)  !! Ur
             tmp2(ii,jj,2)= sina*tmp(ii,jj,1) + cosa*tmp(ii,jj,2)  !! Ut
             tmp2(ii,jj,3)=      tmp(ii,jj,3)                      !! Uz
          end do
       end do
       do ii=1,n1            !! . . Rotate Sections from [u1,u2,u3] to [UN,UE,UD]
          do jj=1,n2
             data(ii,jj,1)= cosi*tmp2(ii,jj,1) - sini*tmp2(ii,jj,2)
             data(ii,jj,2)= sini*tmp2(ii,jj,1) + cosi*tmp2(ii,jj,2)
             data(ii,jj,3)=      tmp2(ii,jj,3)
          end do
       end do
    end if
  end subroutine rotatedata

end module rotate


