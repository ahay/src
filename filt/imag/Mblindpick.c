program Picking
  use sep
  use tridiagonal

  implicit none
  type (tris)                         :: slv
  integer                             :: nt, ns, nx, dim, ix, is, it
  real,    dimension (:), allocatable :: trace, pick, ampl
  real,    dimension (:), pointer     :: diag, offd
  integer, dimension (:), allocatable :: n, ipick
  real                                :: s0, ds, eps, lam

  call sep_init ()
  dim = sep_dimension ()
  allocate (n (dim))
  call from_history (n)
  nt = n (1)
  ns = n (2)
  nx = product (n (3:))
  deallocate (n)

  call from_history ("o2",s0) 
  call from_history ("d2",ds)
  call to_history ("n2",1)
  call from_par ("eps",eps,0.01)
  call from_par ("lam",lam,0.01)
  call sep_close ()

  allocate (trace (nt), pick (nt), ipick (nt), offd (nt), diag (nt), ampl (nt))
  call tridiagonal_init (nt, slv)

  offd = -eps
  pick = 0.
  do ix = 1, nx
     ampl = 0.
     ipick = 1.
     do is = 1, ns
        call sep_read (trace)
        trace = trace*trace
        where (trace > ampl)
           ampl = trace
           ipick = is
        end where
     end do

     ampl = ampl*ampl
     ampl = ampl/sum(ampl)
     diag (2:nt-1) = 2.*eps + ampl (2:nt-1)
     diag (1) = eps + ampl (1)
     diag (nt) = eps + ampl (nt)
     if (ix > 1) diag = diag + lam
     call tridiagonal_define (slv, diag, offd)     
     pick = lam*pick + ampl*(s0 + (ipick-1)*ds)
     call tridiagonal_solve (slv, pick)     
     call sep_write (pick)
  end do

  call tridiagonal_close(slv)
  deallocate (trace, pick, ipick, offd, diag, ampl)
  call exit (0)
end program Picking

