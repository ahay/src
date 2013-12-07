!  Compute anisotropic Traveltime
!
!
!  
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


program MTTAni
  use rsf
 
  implicit none
  type (file)                      :: two, three,  commonoffset
  type(axa)  :: at,az,ax     ! cube axes 
  integer                          :: n1, n2,cut, offset
  real, allocatable :: vol (:,:,:), slice (:,:), cmofst (:,:)

  call sf_init()            ! initialize RSF
  two = rsf_input("in")	    ! two is the input file name 
  three = rsf_output("out") ! three is the output file name 
  commonoffset = rsf_output("common") ! common_offset is the output file name for common offset section   


  if (sf_float /= gettype(two)) call sf_error("Need float type")
  call from_par("cut",cut)   ! command-line parameter 
  call from_par("offset",offset)   ! command-line parameter 
  call from_par(two,"n1",n1) ! get the lenth of each trace (sample number)
  n2 = filesize(two,1)       ! get the number of traces


  call iaxa(two,az,1);       ! read the axis info for the first axis (time) from input file 
  call oaxa(three,az,1);     ! write the axis info for the first axis (time) for the output file
  call oaxa(commonoffset,az,1);     ! write the axis info for the first axis (time) for the output file

  call to_par(commonoffset,"n2",cut)     ! write the sample number for the 2nd axis for the output file 
  call to_par(commonoffset,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(commonoffset,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 


  call to_par(three,"n2",n2/cut)  ! write the sample number for the 2nd axis for the output file 
  call to_par(three,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file 
  call to_par(three,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file 
  call to_par(three,"n3",cut)     ! write the sample number for the 3rd axis for the output file 
  call to_par(three,"o3",1)      ! write the starting point for the 3rd axis (here, starting value is 0) for the output file 
  call to_par(three,"d3",1)      ! write the increament for the 3rd axis (here, increament is 1) for the output file 



  allocate (slice (n1,n2)) ; slice=0.; call rsf_read(two,slice)  ! allocate the 2D array for reading the input file,ile into array "slice"
  allocate (cmofst (n1,cut)) ; cmofst=0.; 			 ! allocate the 2D array for the output file (common offset section)
  allocate (vol (n1,n2/cut,cut)) ; vol=0.; 			 ! allocate the 3D array for 3D output data


  vol=reshape(slice,(/ n1,n2/cut,cut /)) 			! reshape 2D to 3D

  cmofst=vol(:,offset,:)

!##############################################################################################################################

1/3 (n1^2 + (Sqrt[
      3] (-n1^2 (C55 n1^2 + 2 C45 n1 n2 + C44 n2^2 + C33 n3^2) - 
        n1^2 (C66 n1^2 + 2 C26 n1 n2 + C22 n2^2 + C44 n3^2) + 
        2/3 n1^2 (C11 n1^2 + C55 n1^2 + C66 n1^2 + 2 C16 n1 n2 + 
           2 C26 n1 n2 + 2 C45 n1 n2 + C22 n2^2 + C44 n2^2 + 
           C66 n2^2 + (C33 + C44 + C55) n3^2)))/(\[Sqrt](((C36 + 
             C45) n1 n3 + (C23 + C44) n2 n3)^2 + ((C13 + 
             C55) n1 n3 + (C36 + C45) n2 n3)^2 - (C55 n1^2 + 
           2 C45 n1 n2 + C44 n2^2 + C33 n3^2) (C66 n1^2 + 
           2 C26 n1 n2 + C22 n2^2 + C44 n3^2) + (C16 n1^2 + 
          n2 ((C12 + C66) n1 + C26 n2) + 
          C45 n3^2)^2 - (C55 n1^2 + 2 C45 n1 n2 + C44 n2^2 + 
           C33 n3^2) (C11 n1^2 + 2 C16 n1 n2 + C66 n2^2 + 
           C55 n3^2) - (C66 n1^2 + 2 C26 n1 n2 + C22 n2^2 + 
           C44 n3^2) (C11 n1^2 + 2 C16 n1 n2 + C66 n2^2 + C55 n3^2) + 
        1/3 (C11 n1^2 + C55 n1^2 + C66 n1^2 + 2 C16 n1 n2 + 
           2 C26 n1 n2 + 2 C45 n1 n2 + C22 n2^2 + C44 n2^2 + 
           C66 n2^2 + (C33 + C44 + C55) n3^2)^2)))



1/3 (n1^2+(Sqrt[3] (-n1^2 (C55 n1^2+2 C45 n1 n2+C44 n2^2+C33 n3^2)-n1^2 (C66 n1^2+2 C26 n1 n2+C22 n2^2+C44 n3^2)+2/3 n1^2 (C11 n1^2+C55 n1^2+C66 n1^2+2 C16 n1 n2+2 C26 n1 n2+2 C45 n1 n2+C22 n2^2+C44 n2^2+C66 n2^2+(C33+C44+C55) n3^2)))/(\[Sqrt](((C36+C45) n1 n3+(C23+C44) n2 n3)^2+((C13+C55) n1 n3+(C36+C45) n2 n3)^2-(C55 n1^2+2 C45 n1 n2+C44 n2^2+C33 n3^2) (C66 n1^2+2 C26 n1 n2+C22 n2^2+C44 n3^2)+(C16 n1^2+n2 ((C12+C66) n1+C26 n2)+C45 n3^2)^2-(C55 n1^2+2 C45 n1 n2+C44 n2^2+C33 n3^2) (C11 n1^2+2 C16 n1 n2+C66 n2^2+C55 n3^2)-(C66 n1^2+2 C26 n1 n2+C22 n2^2+C44 n3^2) (C11 n1^2+2 C16 n1 n2+C66 n2^2+C55 n3^2)+1/3 (C11 n1^2+C55 n1^2+C66 n1^2+2 C16 n1 n2+2 C26 n1 n2+2 C45 n1 n2+C22 n2^2+C44 n2^2+C66 n2^2+(C33+C44+C55) n3^2)^2)))




!##############################################################################################################################


  
  write (0,*) 'dim cmof', shape(cmofst),'            dim vol', shape(vol)



  call rsf_write(commonoffset,cmofst) 	! write out the common offset section
  call rsf_write(three,vol) 	! write out the 3D volume


 call exit(0)

end program MR2to3

