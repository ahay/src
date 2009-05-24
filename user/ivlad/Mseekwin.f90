! Test for the seek wrapper function in the F90 API
! Works like sfwindow in 1-D. Extracts consecutive sequence of values in N-d.
! The following commands should output the integer sequence from 10 to 19:
! sfmath n1=20 o1=0 d1=1 type=float output=x1 > junk.rsf
! <junk.rsf sfseekwin | sfdisfil col=10 number=n
! Cannot take input from a pipe because stdin cannot be seeked through

!  Copyright (C) 2009 Ioan Vlad
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

program seekwin

  use rsf

  implicit none

  type(file) :: in, out
  integer(kind=OFFKIND) :: nseek
  integer :: n_elem_seek, whence, nread, nmax, esize
  real :: o1, d1
  real, dimension(:), allocatable :: val

  call sf_init()
  in  = rsf_input()
  out = rsf_output()

  ! Test is done on floating point data because it is most common
  ! Operator could be overloaded easily
  if( sf_float /= gettype(in) ) call sf_error("Need float input")

  ! Output is 1-D, with d1 inherited from input and o1 computed
  call from_par(in, "o1"   , o1   )
  call from_par(in, "d1"   , d1   )
  call from_par(in, "esize", esize)

  nmax = filesize(in, 0)

  call from_par("whence", whence     , sf_seek_set)
  call from_par("nseek" , n_elem_seek, 10         )
  call from_par("nread" , nread      , 10         )

  if( nread < 1 ) call sf_error(" please read at least one element")
  if( nseek + nread > nmax ) call sf_error("gone over file limit")

  call to_par(out, "n1", nread            )
  call to_par(out, "o1", o1+d1*n_elem_seek)
  call to_par(out, "n2", 1                )

  nseek = n_elem_seek * esize

  ! Argument 2 in sf_seek must be of type integer(kind=OFFKIND)
  ! Argument 3 must be either:
  ! 0(=sf_seek_set): seek from beginning of file
  ! 1(=sf_seek_cur): seek from current position
  ! 2(=sf_seek_end): seek from the end of the file
  call sf_seek(in, nseek, whence)

  allocate(val(nread))

  call rsf_read(in, val)

  call rsf_write(out, val)

end program seekwin
