#
#%
real wide
wide=3.
call initpar()
call vpfile ('junk.v')

#call vpuclip ( -wide-dx/2 ,H-tmax-dt/2+.5, wide+dx/2, H+.05)

call vpuorig ( -.5-wide, 0.);		call doit( wide, 1 )
call vpuorig ( -.5-wide-2.2*wide, 0.);	call doit( wide, 2 )
call vpendplot (); call exit(0); end


subroutine doit(wide, job)
real v,t,z,sqrt,x,H,dt, tmax, wide, zbot,dx,dy
integer	labelsize, iz, job
labelsize = 12
H=9.			# height of top of frame
v=1.;	dx=.04;	dt=.4;	tmax=8.; dy=.3
v= 1.3 * wide/tmax
call vpfat ( 6)
call vpumove (-wide,H-0.); 	call vpudraw (wide,H-0)
call vpumove (0,H-0.);		call vpudraw (0.,H-tmax)
call vpfat ( 3)
call vputext( wide-.3, H-       0.4, labelsize, 0, 'x')	
if( job == 1)
call vputext(     0.15,H- tmax+.03 , labelsize, 0, 't')	
if( job == 2)
call vputext(     0.15,H- tmax+.03 , labelsize, 0, 't+z/v')	

do iz=1,3 {
	z   = .2*iz*tmax
	zbot= .2* 3*tmax
	call vppenup ()
	#for( x=-wide+dx/2.; x<wide; x=x+dx) {
	x=-wide+dx/2.
	while( x<wide ) {
		t = sqrt(z**2+ (x/v)**2 )
		if( t< sqrt(2.)*z) {
			if( job == 2)
				t = t + zbot-z
			call vpupendn (x, H-t)
			}
		x=x+dx
		}
	}
return; end

