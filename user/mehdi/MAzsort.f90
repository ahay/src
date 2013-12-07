!  Extract a common offset section with variable azimuth from data in t,x,y domain and do stacking in azimuth directtion. 
!
!
!  Parameters: 
!
!
!!$  Copyright (C) 2013 UT Austin
!!$  
!!$  This program is free software; you can redistribute it and/or modify
!!$  it under the terms of the GNU General Public License as published by
!!$  the Free Software Foundation; either version 2 of the License, or
!!$  (at your option) any later version.
!!$  
!!$  This program is distributed in the hope that it will be useful,
!!$  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!$  GNU General Public License for more details.
!!$  
!!$  You should have received a copy of the GNU General Public License
!!$  along with this program; if not, write to the Free Software
!!$  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


program MAzsort
  use rsf
 
  implicit none
  type (file)                      :: two, two2, three,sectionf,commonoffset, filtered
  type(axa)  :: at,az,ax     ! cube axes 
  integer                          :: n1, n2,dimen,offset, i,j, radius
  real, allocatable :: vol (:,:,:), slice (:,:,:), cmofst (:,:), zeros1(:,:,:,:)
  real, allocatable :: zeros2(:,:,:,:),section(:,:)
  integer, allocatable :: is(:), js(:)
  real              :: win
  !radius=8
  call sf_init()            ! initialize RSF
  three = rsf_input("in")	    ! two is the input file name 
  two = rsf_output("out") ! three is the output file name 
  two2 = rsf_output("out2") ! three is the output file name 
  commonoffset = rsf_output("common") ! common_offset is the output file name for common offset section   
  sectionf = rsf_output("sectio") ! common_offset is the output file name for common offset section   
  filtered = rsf_output("filt") ! common_offset is the output file name for common offset section   

  if (sf_float /= gettype(three)) call sf_error("Need float type")
  call from_par("dimen",dimen)   ! command-line parameter 
  call from_par("win",win)   ! command-line parameter 
  call from_par("radius",radius)   ! command-line parameter 
  call from_par(three,"n1",n1) ! get the lenth of each trace (sample number)



  call iaxa(three,az,1);       ! read the axis info for the first axis (time) from input file 
  call oaxa(two,az,1);     ! write the axis info for the first axis (time) for the output file
  call oaxa(two2,az,1);     ! write the axis info for the first axis (time) for the output file
  call oaxa(commonoffset,az,1);     ! write the axis info for the first axis (time) for the output file


  call oaxa(filtered,az,1);     ! write the axis info for the first axis (time) for the output file
  !call oaxa(filtered,at,2);     ! write the axis info for the first axis (time) for the output file
 ! call oaxa(filtered,ax,3);     ! write the axis info for the first axis (time) for the output file

  call to_par(commonoffset,"n2",dimen)     ! write the sample number for the 2nd axis for the output file 
  call to_par(commonoffset,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(commonoffset,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 



  call to_par(sectionf,"n1",2001)     ! write the sample number for the 2nd axis for the output file 
  call to_par(sectionf,"o1",0)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(sectionf,"d1",.002)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 

  call to_par(sectionf,"n2",24)     ! write the sample number for the 2nd axis for the output file 
  call to_par(sectionf,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(sectionf,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 


  call to_par(two,"n2",dimen)  ! write the sample number for the 2nd axis for the output file 
  call to_par(two,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(two,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(two,"n3",4)     ! write the sample number for the 3rd axis for the output file 
  call to_par(two,"o3",1)      ! write the starting point for the 3rd axis (here, starting value is 0) for the output file 
  call to_par(two,"d3",1)      ! write the increament for the 3rd axis (here, increament is 1) for the output file 

  call to_par(two2,"n2",dimen)  ! write the sample number for the 2nd axis for the output file 
  call to_par(two2,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(two2,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(two2,"n3",4)     ! write the sample number for the 3rd axis for the output file 
  call to_par(two2,"o3",1)      ! write the starting point for the 3rd axis (here, starting value is 0) for the output file 
  call to_par(two2,"d3",1)      ! write the increament for the 3rd axis (here, increament is 1) for the output file 

  
  call to_par(filtered,"n2",dimen)  ! write the sample number for the 2nd axis for the output file 
  call to_par(filtered,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(filtered,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(filtered,"n3",dimen)     ! write the sample number for the 3rd axis for the output file 
  call to_par(filtered,"o3",1)      ! write the starting point for the 3rd axis (here, starting value is 0) for the output file 
  call to_par(filtered,"d3",1)      ! write the increament for the 3rd axis (here, increament is 1) for the output file 

  allocate (slice (n1,dimen,dimen)) ; slice=0.; call rsf_read(three,slice)  ! allocate the 2D array for reading the input file,ile into array "slice"
  allocate (zeros1 (n1,4,dimen,dimen/2)) ; zeros1=0.
  allocate (zeros2 (n1,4,dimen,dimen/2)) ; zeros1=0.
  allocate (cmofst (n1,dimen)) ; cmofst=0.; 			 ! allocate the 2D array for the output file (common offset section)
  allocate (section (n1,24)) ; section=0.; 			 ! allocate the 2D array for the output file (common offset section)
  allocate (vol (n1,dimen,dimen))  			 ! allocate the 3D array for 3D output data
  allocate (is (24))  			 ! allocate the 3D array for 3D output data
  allocate (js (24))  			 ! allocate the 3D array for 3D output data
  vol=slice


!######## MANUALLY PICK THE DESIRED OFFSETS AND DO STACKING IN AZIMUTHAL DIRECTION (GIVES BETTER S/N THAN PICKING AROUND A COMPUTED CIRCLE) #### 
  is=(/10,9,8,7,6,5,4,4,4,4,5,6,7,8,9,10,11,12,13,13,13,13,12,11/)
  js=(/4,4,4,4,5,6,7,8,9,10,11,12,13,13,13,13,12,11,10,9,8,7,6,5/)


  do i=1,24
	if (6<is(i) .and. is(i)<11) then
	    section(:,i)=(slice(:,is(i),js(i))+slice(:,is(i)-1,js(i))+&
		     slice(:,is(i)+1,js(i)))
	elseif (6<js(i) .and. js(i)<11) then
	    section(:,i)=(slice(:,is(i),js(i))+slice(:,is(i),js(i)-1)+&
		     slice(:,is(i),js(i)+1))
        elseif (((10<is(i) .and. is(i)<13) .and. (4<js(i) .and. js(i)<7)) .or. &
                ((4<is(i) .and. is(i)<7) .and. (10<js(i) .and. js(i)<13))) then

	    section(:,i)=(slice(:,is(i),js(i))+slice(:,is(i)-1,js(i)-1)+&
		     slice(:,is(i)+1,js(i)+1))
        else 
	    section(:,i)=(slice(:,is(i),js(i))+slice(:,is(i)+1,js(i)-1)+&
		     slice(:,is(i)-1,js(i)+1))
        end if
  end do 

!#####  GO AROUND A CIRCLE AND SELECT DESIRED OFFSETS #################################################
!#############   top  ##############################################
  do j=1,dimen 
    do i=dimen/2,1,-1
     if (sqrt((i-dimen/2)**2.+(j-dimen/2)**2.)>radius-win .and. &
             sqrt((i-dimen/2)**2.+(j-dimen/2)**2.)<radius+win) then
        
        zeros1(:,1,j,i)=slice(:,i,j)
	zeros1(:,2,j,i)=slice(:,i-1,j)
	zeros1(:,3,j,i)=slice(:,i,j-1)
	zeros1(:,4,j,i)=slice(:,i-1,j-1)
     else 
        vol(:,i,j)=0.
     end if
    end do
  end do
  
!#############   bottom   ##############################################
  do j=1,dimen
    do i=dimen/2+1,dimen
     if (sqrt((i-dimen/2)**2.+(j-dimen/2)**2.)>radius-win .and. &
             sqrt((i-dimen/2)**2.+(j-dimen/2)**2.)<radius+win) then
        
        zeros2(:,1,j,i-dimen/2)=slice(:,i,j)
	zeros2(:,2,j,i-dimen/2)=slice(:,i-1,j)
	zeros2(:,3,j,i-dimen/2)=slice(:,i,j-1)
	zeros2(:,4,j,i-dimen/2)=slice(:,i-1,j-1)
     else 
        vol(:,i,j)=0.
     end if
    end do
  end do








  cmofst=sum(sum(zeros1,4),2)


  call rsf_write(commonoffset,cmofst) 	! write out the common offset section
  call rsf_write(sectionf,section) 	! write out the common offset section
  call rsf_write(two,zeros1) 	! write out the 3D volume
  call rsf_write(two2,zeros2) 	! write out the 3D volume
  call rsf_write(filtered,vol) 	! write out the 3D volume

 call exit(0)

end program MAzsort





