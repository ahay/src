#include <rsf.h>

#include "ctridiagonal.h"

int main(int argc, char* argv[])
{
    bool inv;
    int nw,nz,nx, iw,ix,iz;
    float dw,dz,dx, vel0,pi, eps, den, beta;
    float complex w, a, *ctime, *tt, *diag1, *diag2, *offd1, *offd2;
    float **depth, **vel, **voff, *time;
    ctris slv;
    sf_file in, out, velocity;

    sf_init ();
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getfloat("beta",&beta)) beta=1./12.;

    if (NULL != sf_getstring ("velocity")) {
	velocity = sf_input("velocity");
    } else {
	velocity = NULL;
    }

    if (inv) { /* modeling */
	if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

	if (!sf_getint("nt",&nw)) sf_error("Need nt=");
	if (!sf_getfloat("dt",&dw)) sf_error("Need dt=");
	dw = 1./(dw*(nw-1));

	sf_putint (out,"n2",nw); 
	sf_putfloat (out,"d2",dw);
	sf_putint (out,"n1",nx); 
	sf_putfloat (out,"d1",dx);
    } else { /* migration */
	f (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");

	if (NULL != velocity) {
	    if (!sf_histint(velocity,"n1",&nz)) 
		sf_error("No n1= in velocity");
	    if (!sf_hisfloat(velocity,"d1",&dz)) 
		sf_error("No d1= in velocity");
	} else {
	    if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	    if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	}

	sf_putint (out,"n2",nx); 
	sf_putfloat (out,"d2",dx);
	sf_putint (out,"n1",nz); 
	sf_putfloat (out,"d1",dz);
    }

    vel = sf_floatalloc2(nz,nx);
    voff = sf_floatalloc2(nz,nx);

    if (NULL != velocity) {
	sf_read (vel[0],sizeof(float),nz*nx,velocity);
    } else { /* constant velocity */
	if (!sf_getfloat ("vel", &vel0)) sf_error("Need vel0=");
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		vel[ix][iz] = vel0;
	    }
	}
    }

    dw *= SF_PI;

    for (ix=0; ix < nx-1; ix++) {
	for (iz=0; iz < nz; iz++) {
	    vel[ix][iz] *= 0.5; /* post-stack exploding reflector */
	    voff[ix][iz] = sqrtf(vel[ix][iz]*vel[ix+1][iz]);
	}
    }
    for (iz=0; iz < nz; iz++) {
	vel[nx-1][iz] *= 0.5;
    }
    dx = 0.25/(dx*dx);

    depth = sf_floatalloc2(nz,nx);
    time = sf_floatalloc(nx);
    ctime = sf_complexalloc (nx);
    tt = sf_complexalloc (nx);
    diag1 = sf_complexalloc (nx);
    diag2 = sf_complexalloc (nx);
    offd1 = sf_complexalloc (nx);
    offd2 = sf_complexalloc (nx);

    if (inv) {
	sf_read(depth[0],sizeof(float),nz*nx,in);
    } else {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		depth[ix][iz] = 0.;
	    }
	}
    }
    
    slv = ctridiagonal_init (nx);

/*           (1. + k2*0.25*(1.-w*dz)/(w*w))/ & */
/*           (1. + k2*0.25*(1.+w*dz)/(w*w))    */

    for (iw = 0; iw < nw; iw++) {
	sf_warning("frequency %d of %d",iw+1, nw);

	if (inv) { /* modeling */
	    w = dw*(eps+I*(iw-1));

	    for (ix=0; ix < nx; ix++) {
		ctime[ix] = depth[ix][nz-1];
	    }

	    for (iz = nz-2; iz >= 0; iz--) { /* step up */

		for (ix=0; ix < nx; ix++) {
		    vel0 = vel[ix][iz];
		    diag1[ix] =   -2.*(beta - (vel0/w-dz)*vel0*dx/w);
		    diag2[ix] = 1.-2.*(beta - (vel0/w+dz)*vel0*dx/w);
		}           

		vel0 = vel[0][iz];
		a = cexpf(-0.5*w/(vel0*sqrtf(dx)));
		
		diag1[0] =    (a-2.)*(beta - (vel0/w-dz)*vel0*dx/w);
		diag2[0] = 1.+(a-2.)*(beta - (vel0/w+dz)*vel0*dx/w);

		vel0 = vel[nx-1][iz];
		a = cexpf(-0.5*w/(vel0*sqrtf(dx)));

		diag1[nx-1] =    (a-2.)*(beta - (vel0/w-dz)*vel0*dx/w);
		diag2[nx-1] = 1.+(a-2.)*(beta - (vel0/w+dz)*vel0*dx/w);

		for (ix=0; ix < nx-1; ix++) {
		    vel0 = voff[ix][iz];
		    offd1[ix] = beta - (vel0/w-dz)*vel0*dx/w;
		    offd2[ix] = beta - (vel0/w+dz)*vel0*dx/w;
		}

		tt[0] = diag1[0]*ctime[0] + offd1[0]*ctime[1];
		for (ix=1; ix < nx-1; ix++) {
		    tt[ix] = 
			offd1[ix-1]*ctime[ix-1] +
			diag1[ix]*ctime[ix] + &
			offd1[ix]*ctime[ix+1];
		}
		tt[nx-1] = offd1[nx-2]*ctime[nx-2] + diag1[nx-1]*ctime[nx-1];
 
		for (ix=0; ix < nx-1; ix++) {
		    ctime[ix] += tt[ix];
		}

		ctridiagonal_define (slv, diag2, offd2);
		ctridiagonal_solve (slv, ctime);

		for (ix=0; ix < nx-1; ix++) {
		    vel0 = vel[ix][iz];
		    ctime[ix] = ctime[ix] * cexpf(-w*dz/vel0) + depth[ix][iz];
		}
	    }

	    for (ix=0; ix < nx-1; ix++) {
		time[ix] = crealf (ctime[ix]);
	    }
	    sf_write (time,sizeof(float),nx,out);
	} else { /* migration */
        w = cmplx(eps*dw,-(iw-1)*dw)
        call sep_read (time)
        ctime = time
        do iz = 1, nz
           depth (iz,:) = depth(iz,:) + real (ctime)

           diag1 =   -2.*(beta - (vel(iz,:)/w-dz)*vel(iz,:)*dx/w)
           diag2 = 1.-2.*(beta - (vel(iz,:)/w+dz)*vel(iz,:)*dx/w)

           a = cexp(-0.5*w/(vel(iz,1)*sqrt(dx)))

           diag1(1) =    (a-2.)*(beta - (vel(iz,1)/w-dz)*vel(iz,1)*dx/w)
           diag2(1) = 1.+(a-2.)*(beta - (vel(iz,1)/w+dz)*vel(iz,1)*dx/w)

           a = cexp(-0.5*w/(vel(iz,nx)*sqrt(dx)))

           diag1(nx) =    (a-2.)*(beta - (vel(iz,nx)/w-dz)*vel(iz,nx)*dx/w)
           diag2(nx) = 1.+(a-2.)*(beta - (vel(iz,nx)/w+dz)*vel(iz,nx)*dx/w)

           offd1 = beta - (voff(iz,:)/w-dz)*voff(iz,:)*dx/w
           offd2 = beta - (voff(iz,:)/w+dz)*voff(iz,:)*dx/w

           tt(1) = diag1(1)*ctime(1) + offd1(1)*ctime(2)
           tt(2:nx-1) = &
                offd1(1:nx-2)*ctime(1:nx-2) + &
                diag1(2:nx-1)*ctime(2:nx-1) + &
                offd1(2:nx-1)*ctime(3:nx)
           tt(nx) = offd1(nx-1)*ctime(nx-1) + diag1(nx)*ctime(nx)
  
           ctime = ctime + tt

           call ctridiagonal_define (slv, diag2, offd2)
           call ctridiagonal_solve (slv, ctime)

           ctime = ctime * cexp(-w*dz/vel(iz,:)) 
        end do
     end if
  end do
  call ctridiagonal_close (slv)
  if (.not. inv) call sep_write (depth)

  deallocate (depth, time, ctime, vel, voff, tt)
  deallocate (diag1, offd1, diag2, offd2)
  call exit (0)
end program FinDif


