!
!
! Lloyds algorthm for velocity selection
!
! Reference velocity selection by a generalized Lloyd method(Clapp 2002)
! Multi-d reference selection (Clapp 2005-6)
!
module lloyd
  use rsf
  use helixcartmod
  implicit none
  real, private, save :: min_r_pct,min_pct,minme,maxme,min_dist
  real, private :: bound2(1000)
  integer,private,save :: ndiv,niter,verb_level
  logical,private,save :: perc_start
  contains
  !INIT FUNCTION


  subroutine init_lloyd()
    !MINIMUM ALLOWED PCT OF ALL SLOWNESS BEFORE GROUP IS DESTROYED
    call from_par("min_region_pct",min_r_pct,1.)
    !MINIMUM ALLOWED DEVIATION BETWEEN CENTROIDS
    call from_par("perc_start",perc_start,.false.)
    call from_par("del_dist",min_dist,.1)
    call from_par("verb",verb_level,0)
    !NUMBER OF ITERATIONS
    call from_par("niter_lloyd",niter,20);
  end subroutine 


  subroutine start_quantile(array,npts,center)
    real :: array(:,:)
    integer :: npts(:)
    real    :: center(:,:)
    real, allocatable :: pts(:,:)
    integer,allocatable :: ipts(:)
    integer :: idim,i
    real    ::mmin,mmax,d

    allocate(pts(maxval(npts),size(npts)))
    allocate(ipts(size(npts)))

    !evenly space for each variable
    do idim=1,size(npts) 
      mmax=maxval(array(:,idim))
      mmin=minval(array(:,idim))
      d=(mmax-mmin)/(npts(idim)-1)
      do i=1,npts(idim)
        pts(i,idim)=mmin+(i-1)*d
      end do
    end do
    do i=1,product(npts)
      call helix2cart(npts,i,ipts)
       do idim=1,size(npts)
         center(i,idim)=pts(ipts(idim),idim)
       end do 
    end do
    do i=1,product(npts)
      if(verb_level>2) write(0,*) "Creating center",center(i,:)
    end do
    deallocate(pts,ipts)
  end subroutine



  integer function lloyd_go(array, center, nn, iregion,nstart)
    real    :: array(:,:),center(:,:)
    integer :: nn,nstart(:),iregion(:)
    real    :: dist_near,distort_local(nn)
    real :: smin,smax,ds,dist_old
    complex :: j(nn)
    integer :: block(size(array)),region_count(nn)
    integer :: i

                                                                                
    call start_quantile(array,nstart,center)
    ndiv=product(nstart)
    call find_region(array,center,iregion,ndiv)
    call compute_center(array,center,iregion,ndiv)
    do i=1,niter
      if(verb_level>1) write(0,*) "Iteration",i," npts ",ndiv
      call del_close(array,center,iregion,ndiv)
      call del_small_num(array,center,iregion,ndiv)
      if(i < niter/2) call split_region(array,center,iregion,ndiv)
    end do
    lloyd_go=ndiv
    if(verb_level>0) write(0,*)  product(nstart),"=IN OUT=",ndiv
  end function

  subroutine split_region(array,near,iregion,ndiv)
    integer   :: ndiv,nmin
    real      :: array(:,:),near(:,:)
    integer   :: iregion(:)
    logical   :: split(ndiv)
    integer   :: ncount(ndiv),i,imax,nmax,ilocs(2),idiv,ir
    real      :: rand(2),maxd
    real      :: dev(ndiv),dist

    if(size(near,1) > ndiv) then
      ncount=0
      dev=0.
      !calc statistics on regions
      do i =1,size(array,1)
        ir=iregion(i)
        ncount(ir)=ncount(ir)+1
        do idiv=1,size(array,2)
          dev(ir)=dev(ir)+(array(i,idiv)-near(ir,idiv))**2
        end do
      end do
      imax=1; nmax=1;maxd=0.
      nmin=2.*min_dist/100.*size(array,1)
      do i=1,ndiv
        dev(i)=dev(i)/ncount(i)
        if(ncount(i)>nmin .and. dev(i) > maxd) then
          maxd=dev(i)
          imax=i
          nmax=ncount(i)
        end if
      end do
      if(verb_level>2) write(0,*) "  Spliting region 1",near(imax,:)
      call random_number(rand)
      ilocs=nint(max(1.,min(nmax*1.,rand*nmax)))
      nmax=0
      do i =1,size(array,1)
        if(iregion(i)==imax) then
          nmax=nmax+1
          if(ilocs(1)==nmax) then
            near(imax,:)=array(i,:)
            if(verb_level>2) write(0,*) "    New center 1",near(imax,:)
          end if
          if(ilocs(2)==nmax) then
            near(ndiv+1,:)=array(i,:)
            if(verb_level>2) write(0,*) "    New center 2",near(ndiv+1,:)
          end if
        end if
      end do
      ndiv=ndiv+1
      call find_region(array,near,iregion,ndiv)
      call compute_center(array,near,iregion,ndiv)
    end if
  end subroutine

  subroutine del_close(array,near,iregion,ndiv)
    integer   :: ndiv
    real      :: array(:,:),near(:,:)
    integer   :: iregion(:)
    integer   :: idiv,ncount,nmin,i,ia
    logical   :: del(ndiv)
    real      :: dist

    nmin=min_r_pct/100.*size(array,1)
    del=.false.
    do  idiv=1,ndiv
      do i=idiv+1,ndiv
        if(.not.del(i)) then
          dist=0
          do ia=1,size(near,2)
            dist=dist+(near(i,ia)-near(idiv,ia))**2
          end do
          dist=sqrt(dist)
          if(dist < min_dist) then 
            del(i)=.true.
            if(verb_level>2) write(0,*) "  Delete close",near(idiv,:),near(i,:)
          end if
        end if 
      end do
    end do

    i=0
    do idiv=1,ndiv
      if(.not. del(idiv)) then
        i=i+1
        near(i,:)=near(idiv,:)
      end if
    end do
    ndiv=i
    call find_region(array,near,iregion,ndiv)
    call compute_center(array,near,iregion,ndiv)
  end subroutine

  subroutine del_small_num(array,near,iregion,ndiv)
    integer   :: ndiv
    real      :: array(:,:),near(:,:)
    integer   :: iregion(:)
    integer   :: idiv,ncount(ndiv),nmin,i
    logical   :: del(ndiv)

    nmin=min_r_pct/100.*size(array,1)
    del=.false.
    ncount=0
    do  i=1,size(array,1)
      ncount(iregion(i))=ncount(iregion(i))+1
    end do
    do idiv=1,ndiv
      if(ncount(idiv) < nmin) then
        del(idiv)=.true.
        if(verb_level>2) write(0,*) "  Delete few",near(idiv,:),ncount(idiv)
      end if
    end do
    i=0
    do idiv=1,ndiv
      if(.not. del(idiv)) then
        i=i+1
        near(i,:)=near(idiv,:)
      end if
    end do
    ndiv=i
    call find_region(array,near,iregion,ndiv)
    call compute_center(array,near,iregion,ndiv)
  end subroutine



  subroutine compute_center(array,near,iregion,ndiv)
    real    :: array(:,:),near(:,:)
    integer :: iregion(:),ndiv
    real    :: sums(ndiv,size(array,2))
    integer :: npts(ndiv)
    integer :: i

    sums=0
    npts=0
    do i=1,size(array,1)
      sums(iregion(i),:)=sums(iregion(i),:)+array(i,:) 
      npts(iregion(i))=npts(iregion(i))+1
    end do
    near=0.
    do i=1,ndiv
      if(npts(i)>0) near(i,:)=sums(i,:)/npts(i)
    end do
  end subroutine

  subroutine calc_dist(array,near,dist)
    real     :: array(:,:)
    real    :: near(:),dist(:)
    integer :: i,id
    dist=0
    do i=1,size(array,1)  
      do id=1,size(array,2)
        dist(i)=dist(i)+(array(i,id)-near(id))**2
      end do
    end do
  end subroutine

  !!this should be done with a n-d tree, but I dont need speed (for now)
  subroutine find_region(array,near,iregion,ndiv)
    real     :: array(:,:),near(:,:)
    integer  :: iregion(:),ndiv
    real     :: dist(size(array,1))
    real     :: dmin(size(array,1))
    integer  :: idiv,i

    call calc_dist(array,near(1,:),dist)
    iregion=1
    dmin=dist
    do idiv=2,ndiv
      call calc_dist(array,near(idiv,:),dist)
      do i=1,size(array,1)
        if(dist(i) < dmin(i)) then
           iregion(i)=idiv
           dmin(i)=dist(i)
        end if
      end do
    end do
  end subroutine

end module lloyd
