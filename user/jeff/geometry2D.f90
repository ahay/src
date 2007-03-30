module geometry2D
!! ---------------------------------------------
!!
!! Geometry calculation module for RWE2D
!! Calculates:
!!		1) Reference Parmeters
!!		2) Coefficient fields
!!		3) Mask for most appropriate reference
!!
!! Written Jeff Shragge May 16, 2006
!!			email: jshragge@gmail.com
!!
!! ---------------------------------------------
  use rsf
  use util	          !! Regular SEP Utilities
  use lloyd		  !! Lloyds module

  implicit none
  integer, private :: nz,nx
  real,    private :: dz,dx
  integer, private :: nref
  real   , private :: eps
  logical, private :: kin,verbose

  !! GEOMETRY
  real,    allocatable, private :: ray(:,:,:)
  real,    allocatable,dimension(:), private  :: g11,g12,g13,g22,g23,g33
  real,    allocatable,dimension(:), private  :: iG11,iG12,iG13,iG22,iG23,iG33
  real,    allocatable,dimension(:), private  :: m11,m12,m13,m22,m23,m33
  real,    allocatable,dimension(:), private  :: gdetr,n1,n2,n3

  integer, allocatable,dimension(:  ), private  :: gmask
  real,    allocatable,dimension(:,:), private  :: gfields

  !! For Reference Calculation
  integer, allocatable,private :: iregion(:),nstart(:)
  real,    allocatable,private :: scaly(:),have_zeros(:)
  real,    allocatable,private :: vslice(:),grefs(:,:),slice(:,:)

  type(file), private :: ref, a
contains

  !----------------------------------------------------------------  
  !!
  !! Geometry initialization routine
  !!
  subroutine geometry_init(nz_in,nx_in,dz_in,dx_in,kin_in,nref_in,verbose_in)
    integer          ::    nz_in,nx_in,                   nref_in, i
    real             ::                dz_in,dx_in
    logical          ::                            kin_in,        verbose_in
    character (len=128)      :: buffer

    nz=nz_in;  nx=nx_in
    dz=dz_in;  dx=dx_in
    nref=nref_in;    kin=kin_in
    verbose=verbose_in;    eps=0.0001;

    !! For Calculating Geometry
    allocate(ray(3,nx,2) )
    allocate(g11(nx),g12(nx),g13(nx),g22(nx),g23(nx),g33(nx))
    allocate(iG11(nx),iG12(nx),iG13(nx),iG22(nx),iG23(nx),iG33(nx))
    allocate(m11(nx),m12(nx),m13(nx),m22(nx),m23(nx),m33(nx))
    allocate(gdetr(nx),n1(nx),n2(nx),n3(nx),vslice(nx) )
    g11=0.;g12=0.;g13=0.;g22=0.;g23=0.; g33=0.;
    iG11=0.;iG12=0.;iG13=0.;iG22=0.;iG23=0.;iG33=0.;
    m11=0.;m12=0.;m13=0.;m22=0.;m23=0.;m33=0.;
    gdetr=0.;n1=0.;n2=0.;n3=0.; vslice=0.;

    allocate( gfields(nx,4), grefs(nref,4), gmask(nx) )
    allocate( slice(nx,4) )
    gfields=0.; grefs=0.;gmask=0;

    !! For Calculating Reference Parameters
    allocate(have_zeros(4));have_zeros=0.;
    allocate(scaly(4)); scaly=0.
    allocate(iregion(nx)); iregion=0;
    allocate(nstart(4));nstart=0;

    call init_lloyd()

    ref = rsf_output("ref.H")
    call to_par(ref,"n1",nref)
    call to_par(ref,"n3",nz-2)
    call to_par(ref,"n2",4)
    call to_par(ref,"o1",0.)
    call to_par(ref,"o2",0.)
    call to_par(ref,"d1",1.)
    call to_par(ref,"d2",1.)
    call settype(ref,sf_float)

    a = rsf_output("a.H")
    call to_par(a,"n1",nx)
    call to_par(a,"n3",nz-2)
    call to_par(a,"n2",4)
    call to_par(a,"o1",0.)
    call to_par(a,"o2",0.)
    call to_par(a,"d1",1.)
    call to_par(a,"d2",1.)
    call settype(a,sf_float)
    
  end subroutine geometry_init

!!----------------------------------------------------------------  
!! 
!! Geometry calcuation routine
!!
!! Outputs:
!!		1) fields - coefficient fields
!!		2) refs	  - references
!!		3) mask	  - mask for approprite reference
!!		4) gnorm  - normalization array
  subroutine geometry_calc(rays,Vel,fields,refs,mask,gnorm)
    !! Main subroutine dealing with geometry
    integer :: iz,inn,mask(:,:),ix
    integer, allocatable :: masky(:,:)
    real    :: rays(:,:,:),Vel(:,:),fields(:,:,:),refs(:,:,:),gnorm(:,:)

    allocate( masky(nx,nz) );    masky=0

    !! 2D Constants
    m12=0.
    m23=0. 
    g12=0.; 
    g22=1.; 
    g23=0.; 

    iG12 = 0. 
    iG22 = 1. 
    iG23 = 0. 

    do iz=2,nz-1 !! For all extrapolation steps

       !! Metric tensor calculation
       call calc_g11_2d(rays(iz-1:iz+1,:,:),g11);
       call calc_g13_2d(rays(iz-1:iz+1,:,:),g13); 
       if ( maxval(abs(g13)) .lt. 0.01) then
          g13=0.
       end if
       call calc_g33_2d(rays(iz-1:iz+1,:,:),g33);
       call calc_determinant(g11,g13,g33,gdetr)

	!! Inverse metric tensor calculation
       do ix=1,nx
          iG11(ix) = g33(ix) / gdetr(ix)**2
       end do
       do ix=1,nx
          iG13(ix) =-g13(ix) / gdetr(ix)**2
       end do
       do ix=1,nx
          iG33(ix) = g11(ix) / gdetr(ix)**2
       end do

	!! Weighted metric tensor calculation
       do ix=1,nx
          m11(ix)=gdetr(ix)*iG11(ix)
       end do
       do ix=1,nx
          m13(ix)=gdetr(ix)*iG13(ix)
       end do
       do ix=1,nx
          m22(ix)=gdetr(ix)
       end do
       do ix=1,nx
          m33(ix)=gdetr(ix)*iG33(ix)
       end do
       
!        call geom_report()

	!! Stretch/rotate parameter calculation
       call calc_n_2d(m11,m12,m13,n1)
       call calc_n_2d(m12,m22,m23,n2)
       call calc_n_2d(m13,m23,m33,n3)

       !! Non-stationary coefficient fields
       do ix=1,nx
          gfields(ix,1) = iG13(ix)/iG33(ix)
       end do
       do ix=1,nx
             gfields(ix,2) = sqrt(1./(Vel(iz,ix)**2*iG33(ix)))
       end do
       do ix=1,nx
          gfields(ix,3) = sqrt(abs(iG11(ix)/iG33(ix) - (iG13(ix)/iG33(ix))**2))
       end do
       do ix=1,nx
          gfields(ix,4)= n3(ix)/m33(ix)
       end do
       !! Compute Reference Fields and Masks using N-D Lloyds algorithm
       call get_references(gfields,grefs,gmask)


       !! Store all fields in arrays
       do ix=1,nx
          gnorm (iz,ix   ) = gdetr(ix)
       end do

       do ix=1,nx
          fields(ix,1:4,iz) = gfields(ix,1:4)
       end do

       do ix=1,nx
          masky (ix,iz) = gmask(ix)
       end do

       do ix=1,nref
          refs  (ix,1:4,iz) = grefs(ix,1:4)
       end do

    end do

    !! Treat ends
    do ix=1,nx
       gnorm(1 ,ix)=gnorm(2   ,ix)
    end do
    do ix=1,nx
       gnorm(nz,ix)=gnorm(nz-1,ix)
    end do

    do ix=1,nx
       fields(ix,1:4,1 )=fields(ix,1:4,2); 
    end do
    do ix=1,nx
       fields(ix,1:4,nz)=fields(ix,1:4,nz-1)
    end do
    
    do ix=1,nref
        refs(ix,1:4,1 )=refs(ix,1:4,2   ) 
    end do
    do ix=1,nref
       refs(ix,1:4,nz)=refs(ix,1:4,nz-1)
    end do

    do ix=1,nx
       masky (ix,1 )=masky (ix,2   ) 
    end do
    do ix=1,nx
       masky (ix,nz)=masky (ix,nz-1)
    end do
    
    do iz=1,nz
       do ix=1,nx
          mask(ix,iz)=masky(ix,iz)
       end do
    end do

    write(0,*) 'Computed reference Parameters and Geometry'
    call coeff_report(fields,refs)
    deallocate ( masky )

  end subroutine geometry_calc
 
 !----------------------------------------------------------------

  subroutine calc_g11_2d(ray,g11)
    integer          :: ix
    real,intent(in ) :: ray(:,:,:)
    real,intent(out) :: g11(:)

    do ix=2,nx-1
    g11(ix) = &
         ( ( ray(2,ix+1,1)-ray(2,ix-1,1) ) / (2.*dx) )**2 + &
         ( ( ray(2,ix+1,2)-ray(2,ix-1,2) ) / (2.*dx) )**2 
     end do
    g11(1 ) = g11(   2) 
    g11(nx) = g11(nx-1)

  end subroutine calc_g11_2d

   !----------------------------------------------------------------

  subroutine calc_g33_2d(ray,g33)
    integer          :: ix
    real,intent(in ) :: ray(:,:,:)
    real,intent(out) :: g33(:)

    do ix=1,nx 
       g33(ix) = &
           ( ( ray(3,ix,1)-ray(1,ix,1) ) / (2.*dz) )**2 + &
           ( ( ray(3,ix,2)-ray(1,ix,2) ) / (2.*dz) )**2 
    end do
  end subroutine calc_g33_2d

  !----------------------------------------------------------------

   subroutine calc_g13_2d(ray,g13)
     integer          :: ix
     real,intent(in ) :: ray(:,:,:)
     real,intent(out) :: g13(:)
     do ix=2,nx-1
        g13(ix) = &
             ( ray(2,ix+1,1)-ray(2,ix-1,1) ) * &
             ( ray(3,ix  ,1)-ray(1,ix  ,1) ) / (4.*dx*dz) + &
             ( ray(2,ix+1,2)-ray(2,ix-1,2) ) * &
             ( ray(3,ix  ,2)-ray(1,ix  ,2) ) / (4.*dx*dz) 
     end do
     g13(1 ) = g13(   2) 
     g13(nx) = g13(nx-1)

   end subroutine calc_g13_2d

  !----------------------------------------------------------------

  subroutine calc_determinant(g11,g13,g33,gdetr)
    integer                        :: ix
    real, dimension(:),intent(in ) :: g11,g13,g33
    real, dimension(:),intent(out) :: gdetr

    do ix=1,nx
       gdetr(ix) = sqrt( abs(g11(ix)*g33(ix)-g13(ix)**2) )	
    end do
  end subroutine calc_determinant

  !----------------------------------------------------------------

  subroutine calc_inverse(g11,g13,g33,iG11,iG12,iG13,iG22,iG23,iG33,gdetr)
    integer           :: ix
    real,dimension(:) :: g11,g13,g33,gdetr
    real,dimension(:) :: iG11,iG12,iG13,iG22,iG23,iG33

    do ix=1,nx
       iG11(ix) = g33(ix) / gdetr(ix)**2
    end do
    do ix=1,nx
       iG12(ix) = 0. 
    end do
    do ix=1,nx
       iG13(ix) =-g13(ix) / gdetr(ix)**2
    end do
    do ix=1,nx
       iG22(ix) = 1. 
    end do
    do ix=1,nx
       iG23(ix) = 0. 
    end do
    do ix=1,nx
       iG33(ix) = g11(ix) / gdetr(ix)**2
    end do

  end subroutine calc_inverse

  !----------------------------------------------------------------

  subroutine calc_n_2d(m11,m12,m13,n)
    integer                        :: ix
    real, intent(in) :: m11(:),m12(:),m13(:)
    real, intent(out):: n(:)
    do ix=2,nx-1
       n(ix) = ( m11(ix+1) - m11(ix-1) )/ (2.*dx)
    end do
    n(1 ) = n(   2) 
    n(nx) = n(nx-1)

    end subroutine calc_n_2d

  !----------------------------------------------------------------  

    subroutine get_references(gfields,grefs,gmask)
      integer :: iref,ii,ix
      real   :: gfields(:,:)
      real ,intent(out)  :: grefs(:,:)
      integer,intent(out):: gmask(:)

      !! gfields (in)  - Input fields to be used (nx,4)
      !! scaly (tmp)    - Scaling function (4)
      !! slice (tmp)     - 1D array used to get refs (nx,4)
      !! iregion (tmp)   - What region has the best attribute set (nx)
      !! nstart (tmp)    - Number of regions (nref)
      !! have_zeros(tmp) - zero mask(4)
      !! gmask(out)    - Output mask (nx)
      !! grefs (out)   - Reference values in (nref,4)

      scaly=0.
!      nstart = (/ 2,8,4,4 /)
      nstart = (/ 1,6,6,1 /)
      
      scaly(1) = 1.
      scaly(2) = (minval(gfields(:,2)))
      scaly(3) = (minval(gfields(:,3))) 
      scaly(4) = 1.

      where (scaly .eq. 0.) 
         have_zeros = 0.
      elsewhere
         have_zeros = 1. 
      end where
      
      slice = gfields !reshape( gfields, (/nx,4/) )

      if (maxval(abs(gfields(:,1)))    .lt. 0.0001) slice(:,1)=1.
      if (maxval(abs(gfields(:,3))-1.) .lt. 0.0001) slice(:,3)=1.
      if (maxval(abs(gfields(:,4)))    .lt. 0.0001) slice(:,4)=1.
  
      do ii=1,4 !! Rescale Data
         if ( have_zeros(ii) .eq. 1.) then
            slice(:,ii)=slice(:,ii)/scaly(ii)
         else
            slice(:,ii)=0. 
         end if
      end do

      iref=lloyd_go(slice,grefs,nref,iregion,nstart) 

      do ii=1,4  
         if ( have_zeros(ii) .eq. 1.) then
            grefs(1:nref,ii)=(grefs(1:nref,ii))*scaly(ii)
         else
            grefs(1:nref,ii)=0. 
         end if
      end do
      if (maxval(abs(gfields(:,1)))   .lt. 0.0001) grefs(:,1)=0.
      if (maxval(abs(gfields(:,3))-1.).lt. 0.0001) grefs(:,3)=1.
      if (maxval(abs(gfields(:,4)))   .lt. 0.0001) grefs(:,4)=0.

      call rsf_write(ref,grefs)
      call rsf_write(a,gfields)

      if(iref < nref)  grefs(iref+1:nref,:)=-1.

      gmask=reshape(iregion,(/nx/) ) !! Mask

    end subroutine get_references

  !----------------------------------------------------------------

  subroutine coeff_report(fields,refs)
    real  ::  fields(:,:,:),refs(:,:,:)
      write(0,*) "COEFFICIENT RANGES"
    write(0,*) "a1  min/max", minval(fields(2:,1,:)),maxval(fields(2:,1,:))
    write(0,*) "a4  min/max", minval(fields(2:,2,:)),maxval(fields(2:,2,:))
    write(0,*) "a5  min/max", minval(fields(2:,3,:)),maxval(fields(2:,3,:))
    write(0,*) "a10 min/max", minval(fields(2:,4,:)),maxval(fields(2:,4,:))
    write(0,*) " "
    write(0,*) "REFERENCE PARAMETERS"
    write(0,*) "b1  min/max", minval(refs(:,1,:)),maxval(refs(:,1,:))
    write(0,*) "b4  min/max", minval(refs(:,2,:)),maxval(refs(:,2,:))
    write(0,*) "b5  min/max", minval(refs(:,3,:)),maxval(refs(:,3,:))
    write(0,*) "b10 min/max", minval(refs(:,4,:)),maxval(refs(:,4,:))

  end subroutine coeff_report
  !----------------------------------------------------------------

  subroutine geom_report( )
    write(0,*) 'METRIC COEFFICIENTS'
    write(0,*) "g11 min/max", minval(g11(2:)),maxval(g11(2:))
    write(0,*) "g13 min/max", minval(g13(2:)),maxval(g13(2:))
    write(0,*) "g33 min/max", minval(g33(2:)),maxval(g33(2:))
    write(0,*) "gdetr min/max", minval(gdetr(2:)),maxval(gdetr(2:))
    write(0,*) " "
    write(0,*) 'INVERSE METRIC COEFFICIENTS'
    write(0,*) "G11 min/max", minval(iG11(2:)),maxval(iG11(2:))
    write(0,*) "G13 min/max", minval(iG13(2:)),maxval(iG13(2:))
    write(0,*) "G22 min/max", minval(iG22(2:)),maxval(iG22(2:))
    write(0,*) "G33 min/max", minval(iG33(2:)),maxval(iG33(2:))
    write(0,*) " "
 
  end subroutine geom_report

  !----------------------------------------------------------------

  subroutine geometry_close()

   ! deallocate(g11,g12,g13,g23,g33,iG11,iG12,iG13,iG22,iG23,iG33,m11,m12,m13)
   ! deallocate(m22,m23,m33,gdetr,n1,n2,n3,gfields,grefs,gmask)
   ! deallocate(have_zeros,scaly,slice,iregion,nstart,vslice)

  end subroutine geometry_close

  !----------------------------------------------------------------  

end module geometry2D
