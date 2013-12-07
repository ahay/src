!  Thickness attribute from reflectors picked by painting
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


program Mthickat
  use rsf
 
  implicit none
  type (file)                      :: reference, reflectors, thickness
  type(axa)  :: at,az     ! cube axes 
  integer                          :: n1, n2,i,j, n1refl
  real, allocatable :: slice (:,:), reflec (:,:), thick (:,:), col1 (:),reflecsorted(:,:)
  integer , dimension(1)                         :: mxloc
  integer , allocatable                         ::  col1sort (:), reflint(:,:)
  real                          :: o1, d1
  call sf_init()            ! initialize RSF

  reflectors = rsf_input("in")	    ! two is the input file name 
  reference = rsf_input("refseis")	    ! two is the input file name 
  thickness = rsf_output("out") ! common_offset is the output file name for common offset section   


  if (sf_float /= gettype(reflectors)) call sf_error("Need float type")
  call from_par(reference,"n1",n1) ! get the lenth of each trace (sample number)
  call from_par(reference,"o1",o1) ! get the lenth of each trace (sample number)
  call from_par(reference,"d1",d1) ! get the lenth of each trace (sample number)
  call from_par(reflectors,"n1",n1refl) ! get the lenth of each trace (sample number)
  n2 = filesize(reference,1)       ! get the number of traces
  d1=.002


  call iaxa(reference,az,1);       ! read the axis info for the first axis (time) from input file 
!  call oaxa(thickness,az,1);     ! write the axis info for the first axis (time) for the output file
!  call oaxa(commonoffset,az,1);     ! write the axis info for the first axis (time) for the output file

  call to_par(thickness,"n1",n1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(thickness,"o1",0)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(thickness,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 


  call to_par(thickness,"n2",n2)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(thickness,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(thickness,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 



  allocate (reflec (n1refl,n2)) ; reflec=0.; call rsf_read(reflectors,reflec)  ! allocate the 2D array for reading the input file,ile into array "slice"
  allocate (reflecsorted (n1refl,n2))
  allocate (reflint (n1refl,n2))
  allocate (slice (n1,n2)) ; slice=0.; call rsf_read(reference,slice)  ! allocate the 2D array for reading the input file,ile into array "slice"
  allocate (thick (n1,n2)) ; thick=0.; 			 ! allocate the 2D array for the output file (common offset section)
  allocate (col1 (n1refl)) 
  allocate (col1sort (n1refl)) 
!  allocate (vol (n1,n2/cut,cut)) ; vol=0.; 			 ! allocate the 3D array for 3D output data
  reflec=reflec-o1
  col1=reflec(:,1)
  !write (0,*) 'col1', col1

  do i=1,n1refl
     mxloc=minloc(col1)
     !mxloc=maxloc(col1) !for sorting in descending order 
     col1sort(i)=mxloc(1)
     col1(mxloc(1))=9999999.
     !col1(mxloc(1))=0.    !for sorting in descending order 
     reflecsorted(i,:)=reflec(col1sort(i),:)
  end do   
   
  !write (0,*) 'col1sort', col1sort

  reflint=(2*int(int(reflecsorted*1000.)/2))/2
  !reflecsorted=real(2*int(int(reflecsorted*1000.)/2))/1000.

  !write (0,*) 'reflecsorted', reflint(:,2)
!  write (0,*) 'reflecsorted', reflint(:,1),'n1',n1,'max reflint',maxval(reflint)
!  write (0,*) 'list',reflint(10,10),'list',reflint(10+1,10)
!  write (0,*) 'dim thick', shape(thick),'min reflint',minval(reflint)


  do i=1,n1refl-1
     do j=1,n2
	 if (n1<reflint(i+1,j)) then
         thick(reflint(i,j):n1,j)=real(reflint(i+1,j)-reflint(i,j))*.002
         else 
         thick(reflint(i,j):reflint(i+1,j),j)=real(reflint(i+1,j)-reflint(i,j))*.002
         end if
 !          write (0,*) 'reflint(i,j)', reflint(i,j),'reflint(i+1,j)', reflint(i+1,j)
     end do
  end do

  
   ! write (0,*) 'thick', thick(100,:)

!  write (0,*) 'list',reflint(i,j),'list',reflint(i+1,j)
  
 ! write (0,*) 'dim thick', shape(thick)




  call rsf_write(thickness,thick) 	! write out the common offset section
!    write (0,*) 'thick 1', thick(1,:),'thick 700', thick(700,:)



 call exit(0)

end program Mthickat





