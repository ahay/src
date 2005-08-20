program RWEaj
  use rsf

  implicit none
  type (file) :: Fi, Fo

  complex, pointer :: rays(:,:)
  real,    pointer :: x(:,:),z(:,:)
  real,    pointer :: h1(:,:),h2(:,:)

  type(axa)    ::at,ag
  integer      ::it,ig
  real         ::gx,gz

  real         :: eps,peps

  call sf_init()
  Fi = rsf_input()
  Fo = rsf_output()

  call iaxa(Fi,at,1)
  call iaxa(Fi,ag,2)

  call to_par(Fo,"esize",4)
  call to_par(Fo,"n3",2)

  call from_par("peps",peps,0.1)

  allocate(rays(at%n,ag%n))
  allocate(  h1(at%n,ag%n))
  allocate(  h2(at%n,ag%n))
  allocate(   x(at%n,ag%n))
  allocate(   z(at%n,ag%n))

  call rsf_read(Fi,rays)
  z=aimag(rays)
  x= real(rays)

  !! h1
  do it=1,at%n-1
     do ig=1,ag%n
        gx = (x(it+1,ig)-x(it,ig))/at%d
        gz = (z(it+1,ig)-z(it,ig))/at%d
        h1(it,ig) = sqrt(gx*gx+gz*gz)
     end do
  end do
  h1(at%n,:) = h1(at%n-1,:)

  !! h2
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
  where(abs(h2)<eps)
     h2=eps
  end where

  call rsf_write(Fo,h1)
  call rsf_write(Fo,h2)

  call exit(0)
end program RWEaj









