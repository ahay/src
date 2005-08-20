program CC2RC
  !! interpolates from cartesian coordinates to ray coordinates
  !! Paul Sava (paul@sep.stanford.edu)
  !! 08/15/03

  use rsf
  
  implicit none
  type (file) :: Fi, Fo, Fr

  logical            :: adj, verb

  real,    pointer :: mapCC(:,:)
  real,    pointer :: mapRC(:,:)
  complex, pointer ::  rays(:,:)

  type(axa)    :: ax,az,at,ag
  integer      :: ix,iz,it,ig, jx,jz
  real         :: x,z

  real         :: xmin,xmax,zmin,zmax
  real         :: fx,fz
  real         :: xo,zo
  real         :: w00,w01,w10,w11, sumw

  call sf_init()
  Fi = rsf_input()
  Fo = rsf_output()
  Fr = rsf_input("rays")

  call from_par("adj",adj,.false.)
  call from_par("verb",verb,.false.)

  call iaxa(Fr,at,1)
  call iaxa(Fr,ag,2)

  if(adj) then
     call from_par("azn",az%n,1)
     call from_par("azo",az%o,0.)
     call from_par("azd",az%d,1.)

     call from_par("axn",ax%n,1)
     call from_par("axo",ax%o,0.)
     call from_par("axd",ax%d,1.)

     call oaxa(Fo,az,1)
     call oaxa(Fo,ax,2)
  else
     call iaxa(Fi,az,1)
     call iaxa(Fi,ax,2)

     call oaxa(Fo,at,1)
     call oaxa(Fo,ag,2)
  end if

  allocate(mapCC(az%n,ax%n)); mapCC=0.
  allocate(mapRC(at%n,ag%n)); mapRC=0.
  allocate( rays(at%n,ag%n))

  call rsf_read(Fr,rays)
  if(adj) then
     call rsf_read(Fi,mapRC)
  else
     call rsf_read(Fi,mapCC)
  end if

  xmin = ax%o
  xmax = ax%o + (ax%n-1)*ax%d
  zmin = az%o
  zmax = az%o + (az%n-1)*az%d

  do ig=1,ag%n
     if(mod(ig,100)==0) write(0,*) ig,ag%n
     do it=1,at%n
        z=aimag(rays(it,ig))
        x= real(rays(it,ig))

        z = min(z,zmax)
        z = max(z,zmin)
        x = min(x,xmax)
        x = max(x,xmin)

        iz=floor((z-az%o)/az%d)+1
        ix=floor((x-ax%o)/ax%d)+1

        iz = min(iz  ,az%n)
        iz = max(iz  , 1  )
        jz = min(iz+1,az%n)

        ix = min(ix  ,ax%n)
        ix = max(ix  , 1  )
        jx = min(ix+1,ax%n)

        zo = az%o+(iz-1)*az%d
        xo = ax%o+(ix-1)*ax%d

        fz= (z-zo)/az%d
        fx= (x-xo)/ax%d

        fz = max(0., fz)
        fz = min(1., fz)
        fx = max(0., fx)
        fx = min(1., fx)

        !! 00   01       --> x
        !!    .         |
        !! 10   11      v
        !!               z
        w00=(1.-fx)*(1.-fz)
        w10=(   fx)*(1.-fz)
        w01=(1.-fx)*(   fz)
        w11=(   fx)*(   fz)

        !! normalize the sum of the weights
        sumw = w00+w10+w01+w11
        w00 = w00 / sumw
        w10 = w10 / sumw
        w01 = w01 / sumw
        w11 = w11 / sumw

        if(adj) then
           mapCC(iz,ix) = mapCC(iz,ix) + w00 * mapRC(it,ig)
           mapCC(jz,ix) = mapCC(jz,ix) + w10 * mapRC(it,ig)
           mapCC(iz,jx) = mapCC(iz,jx) + w01 * mapRC(it,ig)
           mapCC(jz,jx) = mapCC(jz,jx) + w11 * mapRC(it,ig)
        else
           mapRC(it,ig) = mapRC(it,ig) + w00 * mapCC(iz,ix)
           mapRC(it,ig) = mapRC(it,ig) + w10 * mapCC(jz,ix)
           mapRC(it,ig) = mapRC(it,ig) + w01 * mapCC(iz,jx)
           mapRC(it,ig) = mapRC(it,ig) + w11 * mapCC(jz,jx)
        end if

     end do
  end do

  if(adj) then
     call rsf_write(Fo,mapCC)
  else
     call rsf_write(Fo,mapRC)
  end if

  call exit(0)
end program CC2RC






