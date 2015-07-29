program C2Rsinc
  !! interpolates from cartesian coordinates to ray coordinates
  !! using sinc interpolation
  !! Paul Sava (paul@sep.stanford.edu)
  !! 08/15/03

  use rsf

  implicit none
  type (file) :: Fi, Fo, Fr

  logical            :: adj

  real,    pointer :: mapCC(:,:)
  real,    pointer :: mapRC(:,:)
  complex, pointer ::  rays(:,:)

  type(axa)    :: ax,az,at,ag
  integer      :: ix,iz,it,ig, jx,jz
  real         :: x,z

  integer      :: nsz,nsx
  real         :: zz,xx
  real         :: dz,dx
  real         :: r,sincr,dr

  call sf_init()
  Fi = rsf_input()
  Fo = rsf_output()
  Fr = rsf_input("rays")

  call from_par("adj",adj,.false.)
  call from_par("nsz",nsz,11)
  call from_par("nsx",nsx,11)

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

  dr = sqrt(az%d**2 + ax%d**2)

  do ig=1,ag%n
     if(mod(ig,100)==0) write(0,*) ig,ag%n
     do it=1,at%n   
        z=aimag(rays(it,ig))
        x= real(rays(it,ig))

        iz=floor((z-az%o)/az%d)+1
        ix=floor((x-ax%o)/ax%d)+1

        if(iz<   1+nsz) cycle
        if(iz>az%n-nsz) cycle
        if(ix<   1+nsx) cycle
        if(ix>ax%n-nsx) cycle

        do jz=iz-nsz,iz+nsz
           zz=az%o+(jz-1)*az%d
           dz=zz-z
           
           do jx=ix-nsx,ix+nsx
              xx=ax%o+(jx-1)*ax%d
              dx=xx-x
              
              r=sqrt(dz**2 + dx**2)
!              if(abs(r)<epsilon(r)) then
              if(abs(r) < 0.00001*dr) then
                 sincr = 1.0
              else
                 r = r / dr
                 sincr = sin(r)/r
              end if

              if(adj) then
                 mapCC(jz,jx) = mapCC(jz,jx) + mapRC(it,ig)*sincr
              else
                 mapRC(it,ig) = mapRC(it,ig) + mapCC(jz,jx)*sincr
              end if

           end do
        end do

     end do
  end do

  if(adj) then
     call rsf_write(Fo,mapCC)
  else
     call rsf_write(Fo,mapRC)
  end if

  call exit(0)
end program C2RSINC






