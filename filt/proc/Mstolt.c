#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"

int main(int argc, char* argv[])
{
    int inv, nw,nt,nx,ny, iw,ix,iy, nf;
    float dw,dx,dy, vel, x,y, w,st,sq, *str, *out, *trace;
    sf_file in, out;

    sf_init (argc,argv);
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) nx=1;
    if (!sf_histint(in,"n3",&ny)) ny=1;

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");

  call from_par ("pad",nw,nt)
  call from_history ("d1",dw); dw = pi/(dw*nw)


  call from_history ("d2",dx); dx = dx*pi* abs (vel) * 0.5
  call from_history ("d3",dy); dy = dy*pi* abs (vel) * 0.5
  call from_either ("stretch", st, 1.); if (vel < 0) st = 2.-st
  call from_par ("nf",nf,2)
  call sep_close ()

  allocate (out (nw), trace (nw), str (nw))

  call prefilter_init (nf, 3*nw)
  do iy = 1, ny
     y = (iy-1)*dy
     y = y*y
     do ix = 1, nx
        x = (ix-1)*dx
        x = st*(x*x + y)  
        do iw = 1, nw
           w = (iw-1)*dw
           sq = w*w + sign (x, vel)           
           if (sq > 0.) then              
              str (iw) = w*(1.-1./st) + sqrt (sq)/st
           else ! evanescent
              str (iw) = - 2.*dw
           end if
        end do
        
        call int1_init (str, 0., dw, nw, spline_int, nw, nf)

        call sep_read (trace (:nt))
        if (nw > nt) trace (nt+1:nw) = 0.
        call cosft (nw, trace, .true.)

        call prefilter_apply (trace)
        stat = int1_lop (.false.,.false.,trace,out)

        call cosft (nw, out, .false.); out = out*2./(nw-1)        
        call sep_write (out (:nt))
     end do
  end do
  call int1_close ()
  call prefilter_close ()
           
  deallocate (out, str, trace)
  call exit (0)
end program Stolt
