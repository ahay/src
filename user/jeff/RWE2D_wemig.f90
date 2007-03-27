module RWE2D_wemig
!! Migration driver for the RWE2D migration program
!! contains
!!		1) wemig_init		- wemig initialization routine
!!		2) wemig		- main migration driver routine
!!		3) wemig_ssfboth	- split-step fourier propagation
!!		4) wemig_phsboth	- phase-shift propagation
!!		5) wemig_refboth	- interpolation for various phase-shifts
!!		6) rwetaper_init	- a cosine-taper initialization routine
!!		7) rwetaper		- apply cosine-taper to reduce boundary artifacts
  use rsf		!! basic SEPlib routine I/O module
  use util		!! handle axes module
  use fft               !! Basic fft module

  implicit none
  integer, private :: nz,nx,nw,nr,nloop,nref,ntap
  real,    private :: dz,dx,dw,oz,ox,ow,okx,dkx,forward,eps,pi,sixth,c1,c2
  logical, private :: kin,verbose

	!! PRIVATE ARRAYS FOR PROPAGATION
  integer, allocatable,private :: allmask(:)
  real,    allocatable,private :: msk(:),mtt(:),tap(:),ikz2d(:)
  real,    allocatable,private :: allrefs(:,:),allfields(:,:)
  complex, allocatable,private :: swfl(:),swtt(:),rwfl(:),rwtt(:),rdax(:),sdax(:)

  complex, allocatable, private :: lk(:),rk(:),ck(:)
  complex, allocatable, private :: a(:),b(:),c(:),u(:)

  real,                 private :: kmu,   knu,    kro 
  real,    allocatable, private :: mu(  :),nu(  :),ro(  :)
  real,    allocatable, private :: m0(:),n0(:),r0(:)


contains

  !----------------------------------------------------------------
  !! Initialization routine - parameter passing and array allocation
  subroutine wemig_init(ax,az,aw,ar,forward_in,kin_in,nref_in,verbose_in)
    type(myaxis)     :: ax,az,aw,ar
    integer          ::                               nref_in
    real             ::             forward_in
    logical          ::                        kin_in,        verbose_in
    real             :: pi
    pi=acos(-1.)
    forward=forward_in
    nz=az%n;    nx=ax%n;    nw=aw%n
    dz=az%d;    dx=ax%d;    dw=aw%d * 2.*pi
    oz=az%o;    ox=ax%o;    ow=aw%o * 2.*pi

    nref=nref_in

    call from_par("ntap",ntap,5)  	!! Width of cosine taper
    call from_par("nloop",nloop,3)	!! # smoothing iterations for PSPI interpolation
    call from_par("c1",c1,0.5)
    call from_par("c2",c2,0.)
    call from_par("sixth",sixth,0.122996)
    okx=  -pi/    dx
    dkx=2.*pi/(nx*dx)

    kin=kin_in
    verbose=verbose_in;    
    eps=0.0000001
    allocate(lk(nx),rk(nx),ck(nx))
    allocate(a(nx),b(nx-1),c(nx-1),u(nx))

    kmu = 1./ dz         !! 1/ dz
    knu = 1./(2.*dx**2)  !! 1/( 2 dx^2)
    kro = 1./(dz*dx**2)  !! 1/(dz dx^2)

    allocate(m0(nx),n0(nx),r0(nx))
    allocate(mu(nx),nu(nx),ro(nx))

    allocate(rdax(nx),sdax(nx))
    allocate(ikz2d(nx),swfl(nx))
    allocate(swtt(nx),rwfl(nx))
    allocate(rwtt(nx),msk(nx));
    allocate(mtt(nx),tap(nx))
    allocate(allmask  (nx)   )
    allocate(allfields(nx,4) ) 
    allocate(allrefs(nref,4))
    msk=1.;mtt=1.;tap=1.;ikz2d=1.
    call rwetaper_init() !! Initialize the taper program

  end subroutine wemig_init

  !----------------------------------------------------------------
  !! WE migration driver
  !!
  !! Basic routine is:
  !!
  !!	for all frequencies
  !!		for all depths
  !!			Split-step Fourier Propagation
  !!			FFT wavefields
  !!			for all References
  !!				Phase-shift propagation
  !!				IFFT wavefields
  !! 				FD Correction
  !!				Wavefield Interpolation
  !!			end References
  !!			Imaging Condition
  !!		end depths
  !!	end frequencies
  !!
  subroutine wemig(rwf,swf,Pimg,iss,cutfields,refs,cutmask)
    integer          :: iw,iz,ix,ir,iss,cutmask(:,:)!,t1,t2,t3,t4,t5,t6
    real             :: rarg,saxmax,sarg,w
    real             :: Pimg(:,:),refs(:,:,:),cutfields(:,:,:)
    complex,pointer  :: rwf(:,:),swf(:,:),rax(:),sax(:)
!    logical          :: logic !! For timers

	!! Timer parameters
!!$    logic=init_sep_timers()
!!$    logic=setup_next_timer("WLoop",t1)
!!$    logic=setup_next_timer("PHS",t2)
!!$    logic=setup_next_timer("SSF",t3)
!!$    logic=setup_next_timer("REF",t4)

    !! . . Forward Scattering Option
    sarg=-1.; !! if ( forward .ne. 0. ) sarg=1.
    rarg= 1.
    write(0,*) 'STARTING MIGRATION'
    write(0,*) 'USING REFERENCE MAX OF :',nref
	
    do iw=1,nw  !! FREQUENCY LOOP         
!       call start_timer_num(t1)

       w=ow+(iw-1)*dw 		!! Assign frequency
       rax => rwf(:,iw)			!! Choose receiver wavefield frequency
       sax => swf(:,iw)			!! Choose source wavefield frequency
       saxmax=maxval(cabs(sax)) !! Max for normalization

       do iz=2,nz		!! DEPTH LOOP

	  !! Apply Split-step Fourier propagation
!          call start_timer_num(t3)
          call wemig_ssfboth(sax,rax,w,cutmask(:,iz),cutfields(:,:,iz),&
                                       refs(:,:,iz),rarg,sarg)
!          call stop_timer_num(t3)
			
	  !! Assign temp wavefields
          rwtt=rax;rax=0. 
          swtt=sax;sax=0.

	  !! FFT both wavefields
          call fth(.false.,.false.,swtt)
          call fth(.false.,.false.,rwtt)

          do ir=1,nref  !! REFERENCE VELOCITY LOOP
             if (refs(ir,2,iz) .lt. 0.000000001) cycle 
				!! Assign temp wavefields
             rwfl=rwtt
             swfl=swtt

	     !! Phase-shift wavefields
!             call start_timer_num(t2)
             call wemig_phsboth(swfl,rwfl,w,iz,ir,refs(:,:,iz),rarg,sarg)
!             call stop_timer_num(t2)
	
	     !! IFFT wavefields
             call fth( .true.,.false.,swfl)
             call fth( .true.,.false.,rwfl)

             !! Screen
!             call start_timer_num(t5)
             call wemig_psc_coef(cutfields(:,2,iz),cutfields(:,3,iz),&
                                     refs(ir,2,iz),    refs(ir,3,iz))
             call wemig_fds(rwfl,w,iz,rarg)
             call wemig_fds(swfl,w,iz,sarg)
!             call stop_timer_num(t5)

             !! Interpolate wavefields for a PSPI approach
!             call start_timer_num(t4)
             call wemig_refboth(swfl,rwfl,iz,ir,cutmask(:,iz),sax,rax)
!             call stop_timer_num(t4)

          end do !! END Reference loop
			
	  !! Apply taper to wavefields to handle boundary artifaces
          call rwetaper(rax) 
          call rwetaper(sax)

	  !! Apply correlation imaging condition
	  !! (a deconvolution imaging condition is zeroed out)
          Pimg(iz,:)=Pimg(iz,:)+real(rax*conjg(sax))!/&
          !                                  max( maxval(real(sax*conjg(sax))) ,0.01 )

       end do !! END Depth loop

       write(0,*) "SSF",iw,nw,iss,maxval(abs(Pimg(:,:)))
!       call stop_timer_num(t1)

    end do  !! END Frequency Loop
    write(0,*) 'FINISHED OFF SHOT!'

!    call print_timers() !! Output benchmarking statistics
  end subroutine wemig

  !----------------------------------------------------------------
  !!
  !! Split-step Fourier Propagation subroutine 
  !!	
  subroutine wemig_ssfboth(sdax,rdax,w,allmask,allfields,allrefs,rdir,sdir) 
    integer :: ix,allmask(:)
    real    :: w,sdir,rdir,allfields(:,:),allrefs(:,:)
    complex :: sdax(:),rdax(:)

!! Calculate w-x domain correction step
    do ix=1,nx
       ikz2d(ix)=-dz*(w**2*allrefs(allmask(ix),2 )*&
       (allfields(ix,2 )-allrefs(allmask(ix),2))-&
       allrefs(allmask(ix),4)*(allfields(ix,4)-allrefs(allmask(ix),4)))/&
       sqrt(abs(allrefs(allmask(ix),2)**2*w**2-allrefs(allmask(ix),4)**2)+eps ) 
    end do

!! Apply to each wavefield.  Note rdir and sdir control the direction of propagation
!! Should be oppposite for backscattering and same for forward-scattering
    rdax = rdax*cmplx(cos(ikz2d),rdir*sin(ikz2d))
    sdax = sdax*cmplx(cos(ikz2d),sdir*sin(ikz2d))

  end subroutine wemig_ssfboth
  !----------------------------------------------------------------  
  !!
  !! Fourier-domain phase shift
  !!
 
  subroutine wemig_phsboth(sdax,rdax,w,iz,ir,allrefs,rdir,sdir)
    integer :: ix,iz,ir
    real    :: w,kx,rdir,sdir,carg,carg1,allrefs(:,:),pkz
    complex :: kz,rdax(:),sdax(:)

    carg1= (allrefs(ir,2)* w)**2 - allrefs(ir,4)**2

    do ix=1,nx  !! For each Wavenumber
       kx=okx+(ix-1)*dkx  !! Calculate wavenumber
       carg = carg1 - (allrefs(ir,3)*kx)**2  !! Beneath square-root discriminant

       if ( carg .lt. 0) then  !! If discriminant is negative - apply exponential decay
          pkz = -dz*abs(sqrt(-carg))   !! Wavenumber
          rdax(ix)=rdax(ix)*exp(pkz)   !! Apply operator to rwf
          sdax(ix)=sdax(ix)*exp(pkz)   !! Apply operator to swf

       else !! If discriminant is positive - then do propagation

          pkz =-dz*sqrt(carg) !! Wave-number
          rdax(ix)=rdax(ix)*cmplx(cos(pkz),rdir*sin(pkz)) !! Apply operator to rwf
          sdax(ix)=sdax(ix)*cmplx(cos(pkz),sdir*sin(pkz)) !! Apply operator to swf

       end if

       pkz = dz*kx*allrefs(ir,1) !! Non-orthogonal correction term
       rdax(ix)=rdax(ix)*cmplx(cos(pkz),rdir*sin(pkz))
       sdax(ix)=sdax(ix)*cmplx(cos(pkz),sdir*sin(pkz))

    end do !! END Wavenumber 

  end subroutine wemig_phsboth


!!----------------------------------------------------------------
!!
!! Interpolation routine for different phase-shifts
!!
  subroutine wemig_refboth(swfl,rwfl,iz,ir,allmask,sdax,rdax)
    complex :: sdax(:),rdax(:),swfl(:),rwfl(:)
    integer :: iz,ir,ix,ii,allmask(:)

	!! Find out locations where want to take wavefield for given
	!! Reference values
    mtt=0.
    where( allmask .eq. ir) 
       mtt=1.
    end where

    do ii=1,nloop !! Smoothing Iterations
       msk(1 ) = mtt(1   )
       msk(nx) = mtt(nx-1)
       do ix=2,nx-1
          msk(ix)=( mtt(ix-1)+mtt(ix+1) )/2.
       end do
       mtt=msk
    end do
	!! Interpolate Wavefields
	
    sdax = sdax + msk*swfl
    rdax = rdax + msk*rwfl

  end subroutine wemig_refboth

  !----------------------------------------------------------------
!! Taper Initialization Program
  subroutine rwetaper_init()
    integer :: itap,jtap,ix
    real :: pi
    pi=acos(-1.)

    if (ntap .gt. 1) then

       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(   itap) = cos(pi/2.* jtap/(ntap-1))
       end do

       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(nx-itap+1) = cos(pi/2.* jtap/(ntap-1))
       end do

    end if

  end subroutine rwetaper_init

!!----------------------------------------------------------------
!!
!! Apply Taper to wavefield
!!
  subroutine rwetaper(dax)
    complex, pointer  :: dax(:)
    dax = dax * tap
  end subroutine rwetaper

!----------------------------------------------------------------
!! . . REINSERTED FD STUFF 
!----------------------------------------------------------------
  subroutine wemig_psc_coef(aa,bb,a0,b0)
    real    :: aa(:),bb(:)
    real    :: a0,  b0
    integer :: iz

    !write(0,*) "compute PSC coefficients"
    !    n0 = a0*(c1*(aa/a0-1.)-(bb/b0-1.))*(b0/a0)**2
    n0 = -c1 * (b0/a0)**2 *(aa-a0)+ b0/a0   *(bb-b0)
    r0 = 3*c2* (b0/a0)**2 
    m0 = 1.

    m0 = m0 * kmu
    n0 = n0 * knu
    r0 = r0 * kro
    !   write(0,*) 'COEFFs m0',minval(m0),maxval(m0)
    !   write(0,*) 'COEFFs n0',minval(n0),maxval(n0)
    !   write(0,*) 'COEFFs r0',minval(r0),maxval(r0)


  end subroutine wemig_psc_coef
 !----------------------------------------------------------------

  subroutine wemig_ffd_coef(aa,bb,a0,b0)
    real    :: aa(:),bb(:)
    real    :: a0,b0
    integer :: iz
    real,    allocatable :: d1(:),d2(:)
    real    :: tt

    allocate(d1(nx),d2(nx))

    tt = b0
 
!    d1=aa*(bb/aa)**2-a0*(b0/a0)**2
!    d2=aa*(bb/aa)**2-a0*(b0/a0)**4

    d1=a0*((b0/tt)/a0)**2 * ( ( bb/b0)**2* a0/aa    - 1. )
    d2=a0*((b0/tt)/a0)**4 * ( ( bb/b0)**4*(a0/aa)**3- 1. )

    where(abs(d1)<eps ) 
       d1=eps
    end where

    m0 =      1    
    n0 = c1 * d1 ! * tt**2
    r0 = c2 !* d2/d1 * tt**2

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
  subroutine wemig_fds(dax,w,iz,dir)
    real                  :: w,dir
    integer               :: iz
    complex, dimension(:) :: dax
    integer               :: ix

    mu = m0
    nu = n0 / w
    ro = r0 / w**2

    lk = cmplx(ro+sixth*mu, dir*nu)
    rk = cmplx(ro+sixth*mu,-dir*nu)
    ck = cmplx(         mu, 0.)
!    lk = cmplx(ro+sixth*mu, nu)
!    rk = cmplx(ro+sixth*mu,-nu)
!    ck = cmplx(         mu, 0.)

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
       b(i)=b(i)/a(i)
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

end module  RWE2D_wemig
