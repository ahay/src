#include <rsf.h>

#include "helix.h"
#include "mis2.h"
  use bound

  implicit none
  integer                                     :: dim, na,ia, niter, n1, n2
  integer                                     :: maxlag, padin, padout, p1, p2
  real                                        :: a0
  integer                                     :: prec
  integer, dimension (:), allocatable         :: n, m, a
  real,    dimension (:), allocatable, target :: mm, kk
  real,    dimension (:),             pointer :: pmm, pkk
  logical, dimension (:), allocatable, target :: known
  logical, dimension (:),             pointer :: pknown
  type( filter)                               :: aa

  call sep_init (SOURCE)
  dim = sep_dimension ()
  allocate (n (dim), m (dim), a (dim))
  call to_history ("dim",dim)
  call from_history (n)
  call from_par ("prec", prec, 2)
  call from_par ("niter", niter, 100)
  call from_par ("padin", padin, 0)   
  call from_par ("padout", padout, 0)
  n (dim) = n (dim) + padin + padout
  call from_aux ("filt", "n1", na)
  call allocatehelix (aa, na)

  call from_aux ("filt","a0", a0, 1.)
  call from_aux ("filt","n",m)
  a = 1
  call from_aux ("filt","a",a, a)
  call from_aux ("filt","lag", aa%lag, (/ (ia, ia = 1, na) /))
  call boundn (m, n, a, aa)
  call sep_read (aa%flt, "filt"); aa%flt = aa%flt/a0
  call to_history ("lag", aa%lag)
  call sep_close ()

  n2 = product (n)
  if (dim > 1) then
     n1 = product (n(:dim-1))
  else
     n1 = 1
  end if

  p1 = padin*n1+1; p2 = n2 - padout*n1
  allocate (mm (n2), kk (n2), known (n2))
  mm = 0.; kk = 0.; known = .false.
  pmm => mm (p1:p2)
  pkk => kk (p1:p2)
  pknown => known (p1:p2)
  call sep_read (pmm)

  if (exist_file ("mask")) then
     call sep_read (pkk, "mask")
     pknown = (pkk /= 0.)
  else
     pknown = (pmm /= 0.)
  end if
  maxlag = maxval (aa%lag)
  if (padout*n1 >= maxlag) known (n2 - maxlag + 1:n2) = .true.

  if (prec.eq.2) where (pknown) pkk = pmm
  call mis2 (niter, mm, aa, known, prec)
  if (prec.eq.2) where (pknown) pmm = pkk
  call sep_write (pmm)

  deallocate (n, m, a, mm, kk, known)
  call deallocatehelix( aa)
  call exit (0)
end program Miss
