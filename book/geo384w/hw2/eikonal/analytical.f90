program Analytical
! Analytical first-arrival traveltimes. */
use rsf

character(len=FSTRLEN) :: type
integer :: n1, n2, i1, i2
real    :: d1,d2, g1,g2, s,v0, x1,x2,gsq,g,s2,z,d
real, dimension (:), allocatable :: time
type (file) :: in, out

call sf_init() ! initialize Madagascar
in = rsf_input("in")
out = rsf_output("out")

! Get grid dimensions 
call from_par(in,"n1",n1)
call from_par(in,"n2",n2)
call from_par(in,"d1",d1)
call from_par(in,"d2",d2)

call from_par("g1",g1,0.) ! vertical gradient 
call from_par("g2",g2,0.) ! horizontal gradient 

gsq = g1*g1+g2*g2
g = sqrt(gsq)

call from_par("v0",v0)
! initial velocity or slowness squared 

call from_par("s",s,0.)
! shot location at the surface 

call from_par("type",type,"constant")
! case of velocity distribution 

if (0.0 == g1 .and. 0.0 == g2) type="const"
    
allocate (time(n1))

do i2 = 1, n2
   x2 = (i2-1)*d2
   do i1 = 1, n1
      x1 = (i1-1)*d1
      d = x1*x1+(x2-s)*(x2-s)
      
      select case (type(1:1))
      case("s") ! slowness squared 
         s2 = v0+g1*x1+g2*x2
         z = 2.0*d/(s2+sqrt(s2*s2-gsq*d))
         time(i1) = (s2-gsq*z/6.0)*sqrt(z)
      case("v") ! velocity 
         s2 = 2.0*v0*(v0+g1*x1+g2*x2)

!!! CHANGE BELOW !!! 
         time(i1) = hypot(x2-s,x1)/v0
      case("c") ! constant velocity 
      case default
         time(i1) = hypot(x2-s,x1)/v0
      end select
   end do
   call rsf_write(out,time)
end do

call exit (0)
end program Analytical
