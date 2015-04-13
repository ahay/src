!Find the maximum value of the data output is file contains one value

!  Copyright (C) 2015 ZHENDONG Zhang @ KAUST  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

program findmaximum
use rsf

implicit none

type(file)           :: in, out
integer              :: n1, n2, i1, i2
real*4               :: d1, temp
real*4,allocatable   :: fin(:,:), fout(:,:)


call sf_init()

in    = rsf_input("in")
out   = rsf_output("out") 

if (sf_float /= gettype(in)) call sf_error("Need float type")

call from_par(in, "n1", n1)  !Extract parameter n1
n2  = filesize(in, 1)        !Fast way to get n2. It means n2*n3*n4....


call from_par(in, "d1", d1)

allocate (fin(n1,n2), fout(1,1))
fin   =  0.0
fout  =  0.0

call rsf_read(in, fin)

temp = fin(1,1)

do i2 = 1, n2
   
   do i1 = 1, n1

      if (fin(i1,i2) .GT. temp) then

         temp = fin(i1,i2)
      else
         temp = temp
         
      end if
      
   end do
end do

fout(1,1)  =  temp

call to_par(out, "n1", 1)
call to_par(out, "n2", 1)
call to_par(out, "o1", 0)
call to_par(out, "d1", 1)

call rsf_write(out, fout)

deallocate(fin, fout)

end program findmaximum
