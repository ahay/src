program FinDif
  use cband
  use sep

  implicit none
  integer                             :: nw,nh,nx, iw,ix,ih, band, k
  real                                :: w0,dw, h0,dh,dx, pi,w,w2, h,h2,dh2,hdh
  complex                             :: diag, diag2
  complex, dimension (:), allocatable :: in, out, offd, offd2

  call sep_init ()
  call from_history (nx,nw,esize=8)
  call from_history ("d1",dx)
  pi = acos (-1.)
  call from_history ("o2",w0); w0 = w0*2.*pi
  call from_history ("d2",dw); dw = dw*2.*pi
  call from_par ("nh",nh)
  call from_par ("dh",dh); dh = dh/dx; dh2 = dh*dh
  call from_par ("h0",h0); h0 = h0/dx
  call from_par ("band",band,1)
  call sep_close ()

  allocate (in (nx), out (nx), offd (band), offd2 (band))
  call cband_init (nx, band)

  out = 0.
  call sep_write (out) ! nyquist

  do iw = 2, nw
     write (0,*) iw, nw
     w = w0 + (iw-1)*dw; w2 = w*w

     call sep_read (out)
     do ih = 1, nh
        in = out

        h = h0 + (ih-1)*dh; h2 = h*h; hdh = 2.*dh*h

        select case (band)
        case (1)
           diag  = 2.*cmplx(-w*(12.*(h2 - dh2 - hdh) + 5.*w2), &
           -27.*(dh2 + hdh + 4.*h2) - 3.*(5. + dh2 + hdh)*w2)       
           diag2 = 2.*cmplx(-w*(12.*(h2 + 2.*dh2 + 2.*hdh) + 5.*w2), &
           -27.*(3.*dh2 + 3.*hdh + 4.*h2) + 3.*(-5. + dh2 + hdh)*w2)
           
           offd  = cmplx(-w*(12.*(dh2 + hdh - h2) + w2), &
           27.*(dh2 + hdh + 4.*h2) + 3.*(-1. + dh2 + hdh)*w2)
           offd2 = cmplx(-w*(-12.*(2.*dh2 + 2.*hdh + h2) + w2), &
           27.*(3.*dh2 + 3.*hdh + 4.*h2) - 3.*(1. + dh2 + hdh)*w2)
        end select

        call cband_define (diag2, offd2)

        out(1) = diag * in(1) + sum (offd(1:band) * in (2:1+band))
        do k = 2, band
           out(k) = diag * in(k) &
           + sum (offd(1:band) * in (k+1:k+band)) &
           + sum (offd(1:k-1) * in (k-1:1:-1))
        end do
        do k = band + 1, nx - band
           out(k) = diag * in(k) &
           + sum (offd(1:band) * (in (k+1:k+band) +in (k-1:k-band:-1)))
        end do
        do k = nx - band + 1, nx - 1
           out(k) = diag * in(k) &
           + sum (offd(1:nx-k) * in (k+1:nx)) &
           + sum (offd(1:band) * in (k-1:k-band:-1))
        end do
        out (nx) = diag * in (nx) + sum (offd(1:band) * in (nx-1:nx-band:-1))

        call cband_solve (out)
     end do
     call sep_write (out)
  end do
  
  call cband_close ()

  deallocate (in, out, offd, offd2)
  call exit (0)
end program FinDif


