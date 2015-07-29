program RWEab
  use rsf

  implicit none
  type (file) :: Fi, Fo, Fr, Fs

  complex, pointer :: rays(:,:)
  real,    allocatable :: x(:,:),z(:,:)
  real,    allocatable :: h1(:,:),h2(:,:)
  real,    allocatable :: ss(:,:),aa(:,:),bb(:,:)
  real,    allocatable :: ma(:,:),mb(:,:)

  type(axa)    ::at,ag,ar
  integer      ::it,ig
  real         ::gx,gz
  real         :: eps,peps

  logical          :: verb
  integer          :: naref, nbref
  real             :: mina,maxa,dela
  real             :: minb,maxb,delb
  real             :: a0,b0
  integer          :: ia,ib,ii
  complex, allocatable :: ab(:,:)
  real,    allocatable :: mm(:,:)
  real             :: tiny

  call sf_init()
  Fi  = rsf_input()      ! input rays
  Fs  = rsf_input("slo") ! input slowness (in RC)
  Fo  = rsf_output()     ! output variable  a,b,m
  Fr  = rsf_output("abr")! output reference a,b

  call from_par("verb",verb,.false.)

  call from_par("naref",naref,1) ! no of a refs
  call from_par("nbref",nbref,1) ! no of b refs
  if(verb) write(0,*) "naref=",naref
  if(verb) write(0,*) "nbref=",nbref

  call iaxa(Fi,at,1); if(verb) call raxa(at)
  call iaxa(Fi,ag,2); if(verb) call raxa(ag)

!  call to_par(Fo,"esize",4)
  call settype(Fo,sf_float)
  call to_par(Fo,"n3",5)

  call from_par("peps",peps,0.01)

  ar%n=naref*nbref
  ar%o=0.
  ar%d=1.
  call oaxa(Fr,at,1)
  call oaxa(Fr,ar,2)
  call to_par(Fr,"esize",8)

  allocate(rays(at%n,ag%n))
  allocate(   x(at%n,ag%n))
  allocate(   z(at%n,ag%n))
  allocate(  h1(at%n,ag%n))
  allocate(  h2(at%n,ag%n))

  allocate(  ss(at%n,ag%n))
  allocate(  aa(at%n,ag%n))
  allocate(  bb(at%n,ag%n))

  allocate(  ab(at%n,naref*nbref))
  allocate(  mm(at%n,ag%n))
  allocate(  ma(at%n,ag%n))
  allocate(  mb(at%n,ag%n))

  !----------------------------------------------------------------
  ! read rays and slowness
  call rsf_read(Fi,rays)
  z=aimag(rays)
  x= real(rays)
  call rsf_read(Fs,ss)
  !----------------------------------------------------------------

  !! h1 == alpha
  do it=1,at%n-1
     do ig=1,ag%n
        gx = (x(it+1,ig)-x(it,ig))/at%d
        gz = (z(it+1,ig)-z(it,ig))/at%d
        h1(it,ig) = sqrt(gx*gx+gz*gz)
     end do
  end do
  h1(at%n,:) = h1(at%n-1,:)

  !! h2 == J
  do it=1,at%n
     do ig=1+1,ag%n-1
        gx = (x(it,ig+1)-x(it,ig-1))/2./ag%d
        gz = (z(it,ig+1)-z(it,ig-1))/2./ag%d
        h2(it,ig) = sqrt(gx*gx+gz*gz)
     end do
  end do
  h2(:,1)    = h2(:,2)
  h2(:,ag%n) = h2(:,ag%n-1)
  h2(1,:)    = h2(2,:)

  eps=peps*(maxval(abs(h2))+minval(abs(h2)))/2.
  write(0,*) "eps=",eps
  do it=1,at%n
     do ig=1,ag%n
        if(abs(h2(it,ig))<eps) h2(it,ig)=eps
     end do
  end do
!  where(abs(h2)<eps)
!     h2=eps
!  end where

  aa=ss*h1
  bb=h1/h2
  mm=1.
  ma=1.
  mb=1.

  tiny=0.1

  !! compute reference a & b
  do it=1,at%n

     mina=minval(aa(it,:))
     maxa=maxval(aa(it,:))
     dela=(maxa-mina)/naref

     minb=minval(bb(it,:))
     maxb=maxval(bb(it,:))
     delb=(maxb-minb)/nbref

     !! reference medium (a0 and b0)
     ii=1
     do ia=1,naref
        a0    = mina + 0.5*dela + (ia-1)*dela
        do ib=1,nbref
           b0 = minb + 0.5*delb + (ib-1)*delb
           ab(it,ii) = cmplx(a0,b0)
           ii=ii+1
        end do
     end do

     !! mask for a
     do ia=1,naref
        a0    = mina + 0.5*dela + (ia-1)*dela
        where( a0-(0.5+tiny)*dela <= aa(it,:) .and. aa(it,:)<=a0+(0.5+tiny)*dela) 
           ma(it,:)=ia
        end where
     end do

     !! mask for b
     do ib=1,nbref
        b0    = minb + 0.5*delb + (ib-1)*delb
        where( b0-(0.5+tiny)*delb <= bb(it,:) .and. bb(it,:)<=b0+(0.5+tiny)*delb) 
           mb(it,:)=ib
        end where
     end do

  end do

  !! referece medium mask
  ii=1
  do ia=1,naref
     do ib=1,nbref
        
        do it=1,at%n
           do ig=1,ag%n
              if(ma(it,ig)==ia .and. mb(it,ig)==ib) mm(it,ig)=ii
           end do
        end do
!        where(ma==ia .and. mb==ib)
!           mm=ii
!        end where        
        ii=ii+1
     end do
  end do

  !! write output
  call rsf_write(Fo,aa)
  call rsf_write(Fo,bb)
  call rsf_write(Fo,mm)
  call rsf_write(Fo,ma)
  call rsf_write(Fo,mb)
  call rsf_write(Fr,ab)

  call exit(0)
end program RWEab
