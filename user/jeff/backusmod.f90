! Module to implement backus averaging from sfbackus
module backusmod
  use rsf

  implicit none
  integer, private :: navg,nz

contains

  subroutine dobackus_init(navg_in,nz_in)
    integer :: navg_in,nz_in

    navg = navg_in;
    nz   = nz_in;
  end subroutine dobackus_init

!!***********************************************************************
!! . . copy first valid backus value on each end for a distance of navg/2
!!***********************************************************************
  subroutine handlEnds(p)
    integer :: ii,mm
    real :: gt,gb
    real :: p(:)
    
    !! halfwidth of rms window
    mm = (navg-1)/2+1

    !! good value at top and bottom
    gt = p(mm)
    gb = p(nz-mm)
    
    !! fix top
    do ii=1,mm
       p(ii) = gt
    end do

    !! fix bottom
    do ii=nz-mm,nz
       p(ii)=gb
    end do
    
  end subroutine handlEnds

!*******************************************************
!Do Backus averaging of parameter vector p[nz] using a 
!centered window of navg samples.
!
!avg_j = (1/m) ( sum_{i=j-m/s}^{j+m/2} p_i )
!
!where m=navg
!*******************************************************
  subroutine dobackus(p,avg)
    integer :: ii,jj,mm
    real :: p(:),avg(:)
    real :: val

    !! . . Half-width of rms window
    mm = (navg-1)/2
 
    !! . . Loop over output times
    do ii=1,nz
       val=0.
       !! . . Ensure on grid
       if (ii-mm .ge. 1 .and. ii+mm .le. nz) then
          do jj=ii-mm,ii+mm
             val = val+p(jj)
          end do
          avg(ii) = val/navg
       end if

       !! . . Fix end effects
       if (ii .le. mm) avg(ii) = p(ii)
       if (ii .ge. nz-mm) avg(ii) = p(ii)
    end do

  end subroutine dobackus

end module backusmod
