#
#%
subroutine goodie()
integer outfd, output
integer it,nt, ix,nx, i, itlast, iv, nv, iz, iy, ny
real dt,dx, xmin, xmax, tmax, x,t,z, x0,t0, v(2)
call initpar()
from either:    integer nv=1
call hclose()
outfd = output()
call vpfilep( outfd)
nx=18
nt=48
tmax = 9.5;     xmax = 9.; xmin=-2.
dt =  tmax / (nt-1)
dx =  xmax / (nx-1)
t0 = 0;         x0 = dx/2
v(1) = 1.01 * xmax/ sqrt( tmax*tmax - z*z)
#v(2) = 1.4 * v(1)
call vpuorig( -1.+xmin, 0.)
call vpumove( xmin, 9.5);       call vpudraw( xmax, 9.5)
call vpumove(   0., 9.5);       call vpudraw( 0., 9.5-tmax)
call vpcolor(5)
call vputext( xmax-.35, 9.5-.45, 12, 0, 'x')
call vputext( .25     ,  .2    , 12, 0, 't')
#do iv=1,nv {
call vpcolor(6)
#do iz=1,3 {
iz=2
do iy=0,1 {
        z = (-.10 + .33*iz) * tmax + dt/2
        itlast = 1 + z / dt
        do ix=1,nx {
                x = x0 + (ix-1)*dx
                t = sqrt( z*z + (x/v(1))**2)
                x = x0 + (ix-1)*dx  - iy * (xmax - dx)
                it = 1 + t / dt
                t = t0+(it-1)*dt
                do i = min0(itlast+1,it), it {
                        t = t0+(i-1)*dt
                        if( t < tmax ) {
                                if( iy == 0 ) {
                                        call onebox(  x, 9.5-t, dx/2, dt/2 )
                                        if( -x > xmin )
                                       call onebox( -x, 9.5-t, dx/2, dt/2 )
                                        }
                                else
                                        call diamond( -x, 9.5-t, dx/2, dt/2 )
                                }
                        }
                itlast = it
                }
        }

call vpendplot(); call exit(0); end

subroutine onebox( x, y, dx, dy )
real x,y, dx,dy, vx(4),vy(4), edge
edge = .95
vx(1) = x-edge*dx;      vy(1) = y-edge*dy;
vx(2) = x+edge*dx;      vy(2) = y-edge*dy;
vx(3) = x+edge*dx;      vy(3) = y+edge*dy;
vx(4) = x-edge*dx;      vy(4) = y+edge*dy;
call vpumove( vx(4), vy(4))
call vpudraw( vx(1), vy(1))
call vpudraw( vx(2), vy(2))
call vpudraw( vx(3), vy(3))
call vpudraw( vx(4), vy(4))
return; end

subroutine diamond( x, y, dx, dy )
real x,y, dx,dy, vx(4),vy(4), edge
edge = .95
vx(1) = x        ;      vy(1) = y-edge*dy;
vx(2) = x+edge*dx;      vy(2) = y        ;
vx(3) = x        ;      vy(3) = y+edge*dy;
vx(4) = x-edge*dx;      vy(4) = y        ;
call vpumove( vx(4), vy(4))
call vpudraw( vx(1), vy(1))
call vpudraw( vx(2), vy(2))
call vpudraw( vx(3), vy(3))
call vpudraw( vx(4), vy(4))
return; end
