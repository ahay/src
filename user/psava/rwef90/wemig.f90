module wemig
  use rsf
  use fft

  implicit none
  integer, private :: nz,nx,nw,nr
  real,    private :: dz,dx,dw
  real,    private :: oz,ox,ow

  real,    private :: pi
  real,    private :: c1,c2
  real,    private :: sixth
  integer, private :: nloop
  integer, private :: method

  complex, allocatable, private :: lk(:),rk(:),ck(:)
  complex, allocatable, private :: a(:),b(:),c(:),u(:)
  real,                 private :: kmu,    knu,    kro 
  real,    allocatable, private ::  mu(  :  ),nu(  :  ),ro(  :  )
  real,    allocatable, private ::  m0(:,:,:),n0(:,:,:),r0(:,:,:)
  real,                 private :: okx,dkx

  complex, allocatable, private :: wfl(:),wtt(:)
  real,    allocatable, private :: msk(:),mtt(:)

  real,    allocatable, private :: tap(:)
  integer                       :: ntap
contains
  !----------------------------------------------------------------  
  subroutine wemig_init(ax,az,aw,ar,method_in)
    type(axa)           :: ax,az,aw,ar
    integer, intent(in) :: method_in

    pi=acos(-1.)

    method=method_in

    nz=az%n
    nx=ax%n
    nw=aw%n

    dz=az%d
    dx=ax%d
    dw=aw%d * 2*pi

    oz=az%o
    ox=ax%o
    ow=aw%o * 2*pi

    nr=ar%n

    call from_par("c1",c1,0.47)
    call from_par("c2",c2,0.37)
    call from_par("sixth",sixth,0.122996)
    call from_par("nloop",nloop,3)

    allocate(lk(nx),rk(nx),ck(nx))
    allocate(a(nx),b(nx-1),c(nx-1),u(nx))

    kmu = 1./ dz         !! 1/ dz
    knu = 1 /( 2*dx**2)  !! 1/( 2 dx^2)
    kro = 1./(dz*dx**2)  !! 1/(dz dx^2)

    allocate(m0(nz,nx,nr),n0(nz,nx,nr),r0(nz,nx,nr))
    allocate(mu(   nx   ),nu(   nx   ),ro(   nx   ))

    okx=   -pi/    dx
    dkx=2.0*pi/(nx*dx)

    allocate(wfl(nx))
    allocate(wtt(nx))
    allocate(msk(nx)); msk=1.
    allocate(mtt(nx)); mtt=1.

    call from_par("ntap",ntap,1)
    allocate(tap(nx)); tap=1.
    call rwetaper_init()

  end subroutine wemig_init
  !----------------------------------------------------------------
  subroutine rwetaper_init()
    integer :: itap,jtap
    real    :: pi
    integer :: ix

    pi=acos(-1.)

    if(ntap>1) then
       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(   itap  ) = cos(pi/2.* jtap/(ntap-1))
!          tap(   itap  ) = 0.5*sin(pi* jtap/(ntap-1))
       end do
       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(nx-itap+1) = cos(pi/2.* jtap/(ntap-1))
!          tap(nx-itap+1) = 0.5*sin(pi* jtap/(ntap-1))
       end do
    end if
  end subroutine rwetaper_init
  subroutine rwetaper(dax)
    complex, pointer :: dax(:)
    dax = dax * tap
  end subroutine rwetaper
  !----------------------------------------------------------------
  !! WE migration driver
  subroutine wemig_driver(dat,img,aa,bb,mm,a0,b0)
    complex, pointer :: dat(:,:)
    real,    pointer :: img(:,:)
    real    :: aa(:,:),bb(:,:),mm(:,:)
    real    :: a0(:,:),b0(:,:)
    integer :: iw,iz,ix,ir
    real    :: w
    complex, pointer :: dax(:)

    if     (method==3) then 
       call wemig_psc_coef(aa,bb,a0,b0)
       call wemig_coef_report()
       do iw=1,nw;           w=ow+(iw-1)*dw; w=2*w
          write(0,*) "PSC",iw,nw,maxval(abs(img))
          dax => dat(:,iw)
          do iz=1,nz
             wtt=dax
             dax=0.
             do ir=1,nr
                wfl=wtt
                call wemig_phs(wfl,w,iz,ir,a0,b0)
                call wemig_ssf(wfl,w,iz,ir,aa,a0)
                call wemig_fds(wfl,w,iz,ir)
                call wemig_ref(wfl,  iz,ir,mm,dax)
             end do
             call rwetaper(dax)
             img(iz,:) = img(iz,:) + real(dax)
          end do
       end do
    else if(method==2) then
       call wemig_ffd_coef(aa,bb,a0,b0)
       call wemig_coef_report()
       do iw=1,nw;          w=ow+(iw-1)*dw; w=2*w
          write(0,*) "FFD",iw,nw,maxval(abs(img))
          dax => dat(:,iw)
          do iz=1,nz
             wtt=dax
             dax=0.
             do ir=1,nr
                wfl=wtt
                call wemig_phs(wfl,w,iz,ir,a0,b0)
                call wemig_ssf(wfl,w,iz,ir,aa,a0)
                call wemig_fds(wfl,w,iz,ir)
                call wemig_ref(wfl,  iz,ir,mm,dax)
             end do
             call rwetaper(dax)
             img(iz,:) = img(iz,:) + real(dax)
          end do
       end do
    else if(method==1) then
       do iw=1,nw;          w=ow+(iw-1)*dw; w=2*w
          write(0,*) "SSF",iw,nw,maxval(abs(img))
          dax => dat(:,iw)
          do iz=1,nz
             wtt=dax
             dax=0.
             do ir=1,nr
                wfl=wtt
                call wemig_phs(wfl,w,iz,ir,a0,b0)
                call wemig_ssf(wfl,w,iz,ir,aa,a0)
                call wemig_ref(wfl,  iz,ir,mm,dax)
             end do
             call rwetaper(dax)
             img(iz,:) = img(iz,:) + real(dax)
          end do
       end do
    else
       call wemig_xfd_coef(aa,bb)
       call wemig_coef_report()
       do iw=1,nw;          w=ow+(iw-1)*dw; w=2*w
          write(0,*) "XFD",iw,nw,maxval(abs(img))
          dax => dat(:,iw)
          do iz=1,nz
             call wemig_ssh(dax,w,iz,aa)
             call wemig_fds(dax,w,iz,1)
             call rwetaper(dax)
             img(iz,:) = img(iz,:) + real(dax)
          end do
       end do
    end if

  end subroutine wemig_driver
  !----------------------------------------------------------------
  subroutine wemig_ref(wfl,iz,ir,mm,dax)
    complex, pointer :: dax(:)
    complex          :: wfl(:)
    real             :: mm(:,:)
    integer          :: iz,ir,ix,ii

    where( nint(mm(iz,:))==ir) 
       mtt=1.
    elsewhere
       mtt=0.
    end where

    do ii=1,nloop
       msk(1)=mtt(1)
       do ix=1+1,nx-1
          msk(ix) = 0.5*( mtt(ix-1)+mtt(ix+1) )
       end do
       msk(nx)=mtt(nx-1)
       mtt=msk
    end do

    dax = dax + msk*wfl
  end subroutine wemig_ref
  !----------------------------------------------------------------  
  !! Fourier-domain phase shift
  subroutine wemig_phs(dax,w,iz,ir,a0,b0)
    complex :: dax(:)
    real    :: a0(:,:),b0(:,:)
    real    :: w,kx,arg
    integer :: ix,iz,ir
    complex :: ikz

    call fth(.false.,.false.,dax)
    do ix=1,nx
       kx=okx+(ix-1)*dkx
       arg=1.0- ( (b0(iz,ir)*kx) / (w*a0(iz,ir)) )**2
       if (arg < 0.0 ) then
          ikz = cmplx(   w*a0(iz,ir)*sqrt(-arg),0.)
       else
          ikz = cmplx(0.,w*a0(iz,ir)*sqrt(+arg)   )
       end if
       dax(ix) = dax(ix) * cexp( ikz*(-dz))
    end do
    call fth( .true.,.false.,dax)

  end subroutine wemig_phs
  !----------------------------------------------------------------
  !! space-domain phase shift
  subroutine wemig_ssh(dax,w,iz,aa)
    complex :: dax(:)
    real    :: aa(:,:)
    real    :: w
    integer :: iz

    dax = dax * cexp(cmplx(0.,-w*aa(iz,:) * dz))
  end subroutine wemig_ssh
  !----------------------------------------------------------------
  !! split-step Fourier
  subroutine wemig_ssf(dax,w,iz,ir,aa,a0) 
    complex :: dax(:)
    real    :: aa(:,:)
    real    :: a0(:,:)
    real    :: w
    integer :: iz,ir

    dax = dax * cexp(cmplx(0.,-w*(aa(iz,:)-a0(iz,ir)) * dz))
  end subroutine wemig_ssf
  !----------------------------------------------------------------
  subroutine wemig_xfd_coef(aa,bb)
    real    :: aa(:,:),bb(:,:)
    integer :: iz

    m0(:,:,1) = 1.
    n0(:,:,1) =-c1 * bb**2 / aa
    r0(:,:,1) = c2 * bb**2 / aa**2

    m0 = m0 * kmu
    n0 = n0 * knu
    r0 = r0 * kro

  end subroutine wemig_xfd_coef
  !----------------------------------------------------------------
  subroutine wemig_psc_coef(aa,bb,a0,b0)
    real    :: aa(:,:),bb(:,:)
    real    :: a0(:,:),b0(:,:)
    integer :: iz,ir
    real    :: tt
    real    :: boa

    write(0,*) "compute PSC coefficients"
    tt = abs(maxval(b0))

    m0 = 1 / tt**2
    do iz=1,nz
       do ir=1,nr
          boa = (b0(iz,ir)/tt)/a0(iz,ir)

          n0(iz,:,ir) =    c1*( aa(iz,:)/a0(iz,ir) - 1.) - ( bb(iz,:)/b0(iz,ir) - 1.)
          n0(iz,:,ir) = n0(iz,:,ir) * a0(iz,ir) * boa**2
          r0(iz,:,ir) = 3.*c2 * boa**2

       end do
    end do

    m0 = m0 * kmu
    n0 = n0 * knu
    r0 = r0 * kro
    write(0,*) "done PSC coefficients"

  end subroutine wemig_psc_coef
  !----------------------------------------------------------------
  subroutine wemig_ffd_coef(aa,bb,a0,b0)
    real    :: aa(:,:),bb(:,:)
    real    :: a0(:,:),b0(:,:)
    integer :: iz,ir
    real    :: tt
    real    :: boa
    real,    allocatable :: d1(:,:,:),d2(:,:,:)

    allocate(d1(nz,nx,nr),d2(nz,nx,nr))

    write(0,*) "compute FFD coefficients"
    tt = abs(maxval(b0))

    do iz=1,nz
       do ir=1,nr
          boa = (b0(iz,ir)/tt)/a0(iz,ir)
          d1(iz,:,ir) = a0(iz,ir) * boa**2 * ( (bb(iz,:)/b0(iz,ir))**2 * (a0(iz,ir)/aa(iz,:))     - 1. )
          d2(iz,:,ir) = a0(iz,ir) * boa**4 * ( (bb(iz,:)/b0(iz,ir))**4 * (a0(iz,ir)/aa(iz,:))**3  - 1. )
       end do
    end do
    where(abs(d1)<epsilon(d1)) 
       d1=epsilon(d1)
    end where

    m0 =      d1    / tt**2
    n0 =-c1 * d1**2
    r0 = c2 * d2

    m0 = m0 * kmu
    n0 = n0 * knu
    r0 = r0 * kro

    deallocate(d1,d2)
  end subroutine wemig_ffd_coef
  !----------------------------------------------------------------
  subroutine wemig_coef_report()

    write(0,*) minval(m0),"< mu <",maxval(m0)
    write(0,*) minval(n0),"< nu <",maxval(n0)
    write(0,*) minval(r0),"< ro <",maxval(r0)

  end subroutine wemig_coef_report
  !----------------------------------------------------------------  
  !! finite-differences solver
  subroutine wemig_fds(dax,w,iz,ir)
    real                  :: w
    integer               :: iz,ir
    complex, dimension(:) :: dax
    integer               :: ix

    mu = m0(iz,:,ir)
    nu =-n0(iz,:,ir) / w
    ro = r0(iz,:,ir) / w**2

    lk = cmplx(ro+sixth*mu, nu)
    rk = cmplx(ro+sixth*mu,-nu)
    ck = cmplx(         mu, 0.)

    !!             |                      |
    !!     b       |          a           |     c
    !!  ro + i nu  |  mu - 2(ro + i nu )  |  ro + i nu
    !!             |                      |
    a = ck - 2.*lk  
    b = lk(1  :nx-1)
    c = lk(1+1:nx  )

    do ix=2,nx-1
       u(ix)=           rk(ix) *(dax(ix-1)+dax(ix+1)) & 
       +    (ck(ix)-2.0*rk(ix))* dax(ix)
       !!             |                      | 
       !!  ro - i nu  |  mu - 2 (ro - i nu)  | ro - i nu
       !!             |                      | 
    end do
    ix= 1; u(ix)=(ck(ix)-2.0*rk(ix))*dax(ix) + rk(ix)*dax(ix+1)
    ix=nx; u(ix)=(ck(ix)-2.0*rk(ix))*dax(ix) + rk(ix)*dax(ix-1)

    call thr(a,b,c,u,nx)
    dax = u

  end subroutine wemig_fds
  !----------------------------------------------------------------
  !! tridiagonal solver
  subroutine thr(a,b,c,v,n)
    integer :: n
    complex :: a(n),b(n-1),c(n-1),v(n)
    integer :: i

    b(1)=b(1)/a(1)
    v(1)=v(1)/a(1)
    do i=2,n
       a(i)=a(i)-c(i-1)*b(i-1)
       if (i<n) then 
          b(i)=b(i)/a(i)
       end if
       v(i)=(v(i)-c(i-1)*v(i-1))/a(i)
    end do

    do i=n-1,1,-1
       v(i)=v(i)-b(i)*v(i+1)
    end do
  end subroutine thr
  !----------------------------------------------------------------
  subroutine attribute(x)
    complex :: x(:)
    write(0,*) minval( real(x)),maxval( real(x))
    write(0,*) minval(aimag(x)),maxval(aimag(x))
  end subroutine attribute
  !----------------------------------------------------------------
  subroutine cwrite(x)
    complex :: x
    write(0,*) '(',real(x),',',aimag(x),')'
  end subroutine cwrite
  !------------------------------------------------------------
end module wemig
