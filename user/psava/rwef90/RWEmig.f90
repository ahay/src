program RWEmig
  use rsf
  use wemig

  implicit none
  type (file) :: Fi, Fo, Fm, Fr

  type(axa) :: ax,az,aw,ar

  integer :: method
  integer :: iz

  logical :: verb

  real,    allocatable :: aa(:,:),bb(:,:),mm(:,:)
  complex, allocatable :: ab(:,:)
  real,    allocatable :: a0(:,:),b0(:,:)
  complex, pointer     :: dat(:,:) !! data
  real,    pointer     :: img(:,:) !! image

  call sf_init()
  Fi = rsf_input()
  Fo = rsf_output()
  Fm = rsf_input("abm")
  Fr = rsf_input("abr")

  call from_par("verb",verb,.false.)

  call from_par("method",method,0)

  call iaxa(Fi,ax,1)     !! x = 'space' axis (can be angle)
  call iaxa(Fi,aw,2)     !! w = frequency
  call iaxa(Fm,az,1)     !! z = extrapolation axis (can be time)
  call iaxa(Fr,ar,2)     !! references axis 
  if(method==0) ar%n=1   !! pure finite-differences

  if(verb) then
     write(0,*) ax%n,ax%o,ax%d
     write(0,*) az%n,az%o,az%d
     write(0,*) aw%n,aw%o,aw%d
     write(0,*) ar%n,ar%o,ar%d
  end if
  
  call oaxa(Fo,az,1)
  call oaxa(Fo,ax,2)
  call settype(Fo,sf_float)
!  call to_par(Fo,"esize",4)

  write(0,*) "read data"
  allocate(dat(ax%n,aw%n))
  call rsf_read(Fi,dat)

  !----------------------------------------------------------------
  write(0,*) "read ABM"
  allocate(aa(az%n,ax%n))
  allocate(bb(az%n,ax%n))
  allocate(mm(az%n,ax%n))
  call rsf_read(Fm,aa)
  call rsf_read(Fm,bb)
  call rsf_read(Fm,mm)
  if(verb) then
     write(0,*) "aa",minval(aa),maxval(aa)
     write(0,*) "bb",minval(bb),maxval(bb)
     write(0,*) "mm",minval(mm),maxval(mm)
  end if
  write(0,*) "read ABM ok"
  !----------------------------------------------------------------
  write(0,*) "read ABr"
  allocate(ab(az%n,ar%n))
  allocate(a0(az%n,ar%n))
  allocate(b0(az%n,ar%n))
  call rsf_read(Fr,ab)
  a0= real(ab)
  b0=aimag(ab)
  if(verb) then
     write(0,*) "a0",minval(a0),maxval(a0)
     write(0,*) "b0",minval(b0),maxval(b0)
  end if
  write(0,*) "read ABr ok"
  !----------------------------------------------------------------
  write(0,*) "allocate image"
  allocate(img(az%n,ax%n))
  img=0.0

  write(0,*) "init modeling"
  call wemig_init(ax,az,aw,ar,method)

  write(0,*) "run modeling"
  call wemig_driver(dat,img,aa,bb,mm,a0,b0)

  write(0,*) "write image"
  call rsf_write(Fo,img)
  !----------------------------------------------------------------

  deallocate(dat,img)
end program RWEmig
