/* Testing forward interpolation in 1-D. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
 
#include <rsf.h>

#include "int1.h"
#include "interp.h"
#include "interp_hale.h"
  use interp_cube
  use interp_sinc
  use interp_spline
  use prefilter
  use postfilter 

  implicit none
  integer                          :: n, n2, nd, stat, i, nw, i2
  real, dimension (:), allocatable :: mm, tmp, coord, z
  real                             :: o, oo, d, dd, kai
  character (len=128)              :: intp
  logical                          :: spl, adj

  call sep_init ()
  call from_par ("adj",adj,.false.)

  if (adj) then
     call from_history (nd,n2) 
     call from_aux ("model", "n1", n)              
     call to_history ("n1",n)
     call from_history ("d1",dd) 
     call from_aux ("model","d1",d)  
     call to_history ("d1",d)
     call from_history ("o1",oo)
     call from_aux ("model","o1",o)  
     call to_history ("o1",o)
  else
     call from_history (n,n2) 
     call from_aux ("coord", "n1", nd)              
     call to_history ("n1",nd) 
     call from_history ("d1",d) 
     call from_aux ("coord","d1",dd)  
     call to_history ("d1",dd) 
     call from_history ("o1",o)
     call from_aux ("coord","o1",oo)  
     call to_history ("o1",oo) 
  end if

  ! initialize interpolation
  call from_par ("interp",intp)
  call from_par ("nw",nw)
  call from_par ("kai",kai,4.0)
  call sep_close ()

  
  allocate (coord (nd))
  call sep_read (coord,"coord")

  spl = .false.
  select case (intp(1:3))
  case ("lag")
     call int1_init (coord, o, d, n, lg_int, nd, nw)
  case ("cub")
     call int1_init (coord, o, d, n, cube_int, nd, nw)
  case ("mui")
     call int1_init (coord, o, d, n, muir_int, nd, nw)
  case ("hal")
     call hale_init (nw)
     call int1_init (coord, o, d, n, hale_int, nd, nw)
     call hale_close()
  case ("kai")
     call kaiser_init (kai)
     call int1_init (coord, o, d, n, kaiser_int, nd, nw)
  case ("lan")
     call int1_init (coord, o, d, n, lanczos_int, nd, nw)
  case ("cos")
     call int1_init (coord, o, d, n, cosine_int, nd, nw)
  case ("wel")
     call int1_init (coord, o, d, n, welch_int, nd, nw)
  case ("mom")
     spl = .true.
     call prefilter_init (-nw, n)     
     call int1_init (coord, o, d, n, mom_int, nd, nw) 
  case ("spl")
     spl = .true.
     if (adj) then
        call postfilter_init (nw, .false.)
        allocate (tmp (n))
     else
        call prefilter_init (nw, n)
     end if
     call int1_init (coord, o, d, n, spline_int, nd, nw)
  end select
  
  allocate (z (nd), mm (n))

  do i2 = 1, n2
     if (adj) then
        call sep_read (z)
     else     
        call sep_read (mm)
        if (spl) call prefilter_apply (mm)
     end if

     stat = int1_lop (adj,.false.,mm,z)
     
     if (adj) then
        if (spl) then
           call postfilter_apply (mm,tmp)
           call sep_write (tmp)
        else
           call sep_write (mm)
        end if
     else
        call sep_write (z)
     end if
  end do

  call int1_close ()
  if (spl) then
     if (adj) then
        deallocate (tmp)
     else
        call prefilter_close ()
     end if
  end if
  deallocate (z, mm, coord)
  call exit (0)
end program INTERPOLATION



