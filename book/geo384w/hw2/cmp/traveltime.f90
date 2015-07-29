program Traveltime
! Compute traveltime in a V(z) model. 
use rsf

character(len=FSTRLEN) :: type
integer :: ih, nh, it, nt, ir, nr, iter, niter
real    :: h, dh, h0, dt, t0, t2, h2, v2, s, p, hp, tp
integer, allocatable, dimension (:) :: r
real,    allocatable, dimension (:) :: v, t
type (file) :: vel, tim

call sf_init() ! initialize Madagascar

! input and output 
vel = rsf_input("in")
tim = rsf_output("out")
    
! time axis from input 
call from_par(vel,"n1",nt)
call from_par(vel,"d1",dt)

! offset axis from command line 

call from_par("nh",nh,1)    ! number of offsets 
call from_par("dh",dh,0.01) ! offset sampling
call from_par("h0",h0,0.0)  ! first offset

! get reflectors 

call from_par("nr",nr,1) ! number of reflectors

allocate (r(nr))

call from_par("r",r)

call from_par("type",type,"hyperbolic")
! traveltime computation type 

call from_par("niter",niter,10)
! maximum number of shooting iterations 

! put dimensions in output 
call to_par(tim,"n1",nh)
call to_par(tim,"d1",dh)
call to_par(tim,"o1",h0)
call to_par(tim,"n2",nr)

! read velocity 
allocate (v(nt))
call rsf_read(vel,v)

! convert to velocity squared 
v = v*v

allocate(t(nh))

do ir=1, nr
   nt = r(ir)
   t0 = nt*dt ! zero-offset time 
   t2 = t0*t0

   p = 0.0;

   do ih=1, nh
      h = h0+(ih-1)*dh ! offset
      h2 = h*h

      select case (type(1:1))
      case ("h") ! hyperbolic approximation 
         v2 = 0.0
         do it=1, nt
            v2 = v2 + v(it)
         end do
         v2 = v2/nt

         t(ih) = sqrt(t2+h2/v2)
      case("s") ! shifted hyperbola 
         
!!! MODIFY BELOW !!! 

         s = 0.0
         v2 = 0.0
         do it=1, nt
            v2 = v2 + v(it)
         end do
         v2 = v2/nt

         t(ih) = sqrt(t2+h2/v2)
      case("e") ! exact 
		    
!!! MODIFY BELOW !!! 

         do iter=1, niter 
            hp = 0.0
            do it=1,nt
               v2 = v(it)
               hp = hp + v2/sqrt(1.0-p*p*v2)
            end do
            hp = hp*p*dt

!!! SOLVE h(p)=h !!! 
	
         end do

         tp = 0.0
         do it=1, nt
            v2 = v(it)
            tp = tp + dt/sqrt(1.0-p*p*v2)
         end do

         t(ih) = tp
      case default
         call sf_error("Unknown type")
      end select
   end do

   call rsf_write(tim,t)
end do

call exit(0)
end program Traveltime
