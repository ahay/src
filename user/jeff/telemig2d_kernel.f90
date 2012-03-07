module telemig2d_kernel
  
	implicit none
  
  	!! . . Define kernel parameter type
  	type ker_par_type
  	 	integer :: nx
     	real    :: dx,dz,fkx,dkx,ax1,bx1,ax2,bx2,tx ! for finite difference
  	end type ker_par_type
  	type(ker_par_type), private:: ker_par

  	!! . . Define arrays
  	real,   allocatable :: cax(:,:),cbx(:,:)
  	complex,allocatable :: dax(:,:),rbx(:,:),rcx(:,:),ux(:,:),vx(:,:)
  	real   ,allocatable :: KKX(:),xav(:,:)
   
contains

  	!! . . Initialization routine for Kernel propagation
	subroutine shot_ker_init(nth,nx,dx,dz,a,b,trick)
    	integer :: nx,nth,i,j
    	real    :: dx,dz,a,b,trick,pi,a1,a2,b1,b2

    	!! . . Coefficients for two-step propagation
    	a1=0.040315157; b1=0.873981642
    	a2=0.457289566; b2=0.222691983

    	pi=acos(-1.)
    	ker_par%nx=nx
    	ker_par%dx=abs(dx)
    	ker_par%dz=dz
    	ker_par%fkx=-pi/ker_par%dx      
    	ker_par%dkx=2.0*pi/(nx*ker_par%dx)
    	ker_par%ax1=a1*dz/(2.0*ker_par%dx*ker_par%dx)
    	ker_par%bx1=b1/(ker_par%dx*ker_par%dx)
    	ker_par%ax2=a2*dz/(2.0*ker_par%dx*ker_par%dx)
    	ker_par%bx2=b2/(ker_par%dx*ker_par%dx)
    	ker_par%tx=trick   
    
    	allocate ( dax(nx,nth),rbx(nx,nth),rcx(nx,nth),ux(nx,nth),vx(nx,nth))
    	allocate ( cax(nx,nth),cbx(nx,nth),KKX(nx),xav(nx,nth) )
	
   	 	call klinemap(KKX,nx,dx,0)

  	end subroutine shot_ker_init


	!----------------------------------------------------------------
	!! . . Propagate a wavefield one extrapolation step
  	subroutine shot_ker_onestep(id,wld,vel,w,signn)
    	integer :: ix,id
    	real    :: vel(:),w,signn
    	complex :: wld(:,:)

    	!! . . Define RWE coefficients 
    	do ix=1,ker_par%nx
    	   xav(ix,id)=vel(ix) / w
    	end do

    	!! . . Split-step propagation
    	call shot_ker_ssf(id,wld,w,vel,ker_par,signn)

    	!! . . First FFD step
    	call shot_ker_fd (id,wld,xav,xav,ker_par,-signn*ker_par%ax1,ker_par%bx1)

    	!! . . Second FFD step
    	call shot_ker_fd (id,wld,xav,xav,ker_par,-signn*ker_par%ax2,ker_par%bx2)
  	end subroutine shot_ker_onestep

	!---------------------------------------------------------------
	!! . . The split-step extrapolation part
  	subroutine shot_ker_ssf(id,wld,w,vel,ker_par,signn)
    	type(ker_par_type) :: ker_par
    	integer :: ix,id
    	real    :: w,vel(:),signn,phsft
    	complex :: wld(:,:)

    	!! . . Apply appropriate w-x phase-shift
    	do ix=1,ker_par%nx
    	   phsft = -signn * ker_par%dz * w / vel(ix)
    	   wld(ix,id) = wld(ix,id) * cmplx(cos(phsft),sin(phsft))
    	end do

  	end subroutine shot_ker_ssf
  
	!--------------------------------------------------------------
 	 !! . . Define which parameters we're looking at
  	subroutine shot_ker_fd(id,wld,ca,cb,ker_par,ax,bx)
    	type(ker_par_type) ::ker_par
    	integer :: nx,ny,ix,iy,id
    	real    :: ca(:,:),cb(:,:),tx,ax,bx
    	complex :: wld(:,:)

    	nx=ker_par%nx
    	tx=ker_par%tx
	
    	!! . . Inline propagation (x)
    	!! . . Loop over all crossline (y) indices
    	ux (:,id) = wld(:,id)
    	cax(:,id) =  ca(:,id)
    	cbx(:,id) =  cb(:,id)
    	call fd(ux(:,id),vx(:,id),tx,ax,bx,cax(:,id),cbx(:,id),nx,dax(:,id),rbx(:,id),rcx(:,id))
    	wld(:,id)=vx(:,id)

  	end subroutine shot_ker_fd

	!-------------------------------------------------------------
  	!! . . Setup the standard finite difference caller
	subroutine fd(u,v,tr,a,b,ca,cb,nx,ra,rb,rc)
    	integer :: nx,ix
    	real    :: ca(:),cb(:),a,b,tr,signn
    	complex :: url,urr,ror,rol,u(:),v(:),ra(:),rb(:),rc(:)

    	!! Case where ix=1
    	url   = 0. 
    	urr   = u(2)
    	ror   = cmplx( tr+b*(cb(1))**2 , a*ca(1) )
    	rol   = conjg(ror)
    	v (1) = u(1)+ror*(url+urr-2.*u(1))
    	ra(1) = 1.-2.*rol
    	rb(1) = rol
    	rc(1) = 0.

    	do ix=2,nx-1
       		url      = u(ix-1)
      		urr      = u(ix+1)
       		ror      = cmplx( tr+b*(cb(ix))**2 , a*ca(ix) )
       		rol      = conjg(ror)
       		v (ix  ) = u(ix)+ror*(url+urr-2.*u(ix))
       		ra(ix  ) = 1.-2.*rol
       		rb(ix  ) = rol
       		rc(ix-1) = rol
    	end do

    	!! Case where ix=nx
    	url      = u(nx-1)
    	urr      = 0.
    	ror      = cmplx( tr+b*(cb(nx))**2 , a*ca(nx))
    	rol      = conjg(ror)
    	v (nx  ) = u(nx)+ror*(url+urr-2.*u(nx))
    	ra(nx  ) = 1.-2.*rol
    	rb(nx  ) = 0.
    	rc(nx-1) = rol

    	call thryj(ra,rb,rc,v,nx)
  	end subroutine fd

	!-------------------------------------------------------------
	!! . . Tri-diagonal solver
  	subroutine thryj(a,b,c,v,n)
    	complex :: a(:),b(:),c(:),v(:)
    	integer :: ii,n

    	!! . . for i=1
    	b(1) = b(1)/a(1)
    	v(1) = v(1)/a(1)

    	!! for i=2,n-1
    	do ii=2,n-1
     		 a(ii)= a(ii)-b(ii-1)*c(ii-1)
       		v(ii)=(v(ii)-c(ii-1)*v(ii-1))/a(ii)
      		b(ii)= b(ii)/a(ii)
    	end do

    	!! for i=n
    	a(n)= a(n)-b(n-1)*c(n-1)
    	v(n)=(v(n)-c(n-1)*v(n-1))/a(n)

   	 	!! backsubstitution
    	do ii=n-1,1,-1
       		v(ii)=v(ii)-b(ii)*v(ii+1)
    	end do

  	end subroutine thryj

	!------------------------------------------------------------
	!! . . Fourier-domain phase shift correction
 	 subroutine wemig_phs_correction(id,dax,w,vmin,darg)
    	integer :: ikx,iky,id,nx,ny
    	real    :: w,kr2,darg,vmin,kx,ky,avX,avY
    	real    :: tx,bx,ty,by,xav,yav
    	complex :: dax(:,:)
    	real    :: w_v0,sr2,phsft,dz,a1,a2,b1,b2,wv2

    	!! . . Coefficients for filtering
    	a1=-0.040315157; b1=-0.873981642
    	a2=-0.457289566; b2=-0.222691983
    	dz=ker_par%dz
    	nx=ker_par%nx

    	w_v0=w/vmin
    	wv2=w_v0*w_v0

    	!! E1 correction factors
    	avX = 1/w_v0

    	!! . . Loop over all wave-numbers
    	do ikx=1,nx

       		sr2 = ( (KKX(ikx)*xav)**2 )/wv2

       		!! . . Filter out evanescent energy
       		if ( sr2 .gt. 1) then
          		dax(ikx,id) = 0.
       		else

          		!! . . Coefficients
          		tx = (KKX(ikx)*avX)**2
          		bx = (KKX(ikx)*avX)**2

         	 	phsft=-dz*w_v0*(sqrt(1.-sr2)-(1.+a1*tx/(1.+b1*bx)+a2*tx/(1.+b2*bx)))

          		dax(ikx,id)=dax(ikx,id)*cmplx(cos(phsft),darg*sin(phsft))

       		end if
    	end do

  	end subroutine wemig_phs_correction

	!------------------------------------------------------------

  	subroutine klinemap(k,n,d,shift)
    	!! . . Set the correct wavenumbers
    	real, dimension(:) :: k
    	integer, intent(in):: n,shift
    	real,    intent(in):: d
    	integer            :: i,m

    	k=0.
    	if(n/=1) then
       		m = n/2
       		do i=1,m+1
          		k(i+shift)=i-1
       		end do
       		do i=m+2,n
          		k(i+shift)=i-n-1
       		end do
    	else
       		k=0. !! make sure k=0 if n=1
    	end if

    	k = k / real(n*d)       !! 1/length
    	k = 2. * acos(-1.) * k  !! radians

  	end subroutine klinemap

	!----------------------------------------------------------------
  	!! . . Release the kernel memory
  	subroutine shot_ker_release
    	deallocate(dax,rbx,rcx,cax,cbx,ux,vx,xav,KKX)
 	end subroutine shot_ker_release

end module telemig2d_kernel
