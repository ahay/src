!Convert interval velocity(Vertical axis is time) to average velocity

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

program intvel2avg
use rsf

implicit none

type(file)           :: in, out
integer              :: n1, n2, i1, i2
real*4               :: d1, temp
real*4,allocatable   :: velint(:,:), velavg(:,:)

call sf_init()
in    = rsf_input("in")
out   = rsf_output("out") 

if (sf_float /= gettype(in)) call sf_error("Need float type")

call from_par(in, "n1", n1)  !Extract parameter n1
n2  = filesize(in, 1)        !Fast way to get n2. It means n2*n3*n4....

call from_par(in, "d1", d1)

allocate (velint(n1,n2), velavg(n1,n2))
velint  =  0.0
velavg  =  0.0

call rsf_read(in, velint)


do i2 = 1, n2
   temp = 0.0
   do i1 = 1, n1

      temp = temp + velint(i1, i2)
      velavg(i1,i2) = temp / (1.0*i1)
      
   end do
end do

call rsf_write(out, velavg)

deallocate(velint, velavg)

end program intvel2avg
