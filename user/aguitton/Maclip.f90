! Test program

!!$  Copyright (C) 2006 University of Texas at Austin
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
program Clipit
  use rsf

  implicit none
  type (file)                      :: in, out
  integer                          :: n1, n2, i1, i2
  real                             :: clip
  real, dimension (:), allocatable :: trace

  call sf_init()            ! initialize RSF
  in = rsf_input()
  out = rsf_output()

  if (sf_float /= gettype(in)) call sf_error("Need float type")

  call from_par(in,"n1",n1)
  n2 = filesize(in,1)

  call from_par("clip",clip) ! command-line parameter 

  allocate (trace (n1))

  do i2=1, n2                ! loop over traces
     call rsf_read(in,trace)
     
     where (trace >  clip) trace =  clip
     where (trace < -clip) trace = -clip

     call rsf_write(out,trace)
  end do
end program Clipit
