#include <rsf.h>

#include "fint1.h"

int main (int argc, char *argv[])
{
    int nt, nz, nx, iz, it, ix, nw;
    bool slow;
    float t0, dt, z0, dz, t, z;
    float *time, *depth, *vel;
    fint1 fnt;
    sf_file in, out, velocity;

    sf_init(argc, argv);
    in = sf_input("in");
    velocity = sf_input("velocity");
    out = sf_output("out");

    if (!sf_histint (in,"n1",&nt)) sf_error ("No n1= in input");
    if (!sf_histfloat (in,"d1",&dt)) sf_error ("No d1= in input");
    if (!sf_histfloat (in,"o1",&t0)) t0 = 0.;

    nx = sf_leftsize(in,1);

    if (!sf_getint ("nz",&nz)) {
	nz = nt; 
    } else {
	sf_putfloat(out,"n1",nz);
    }
    if (!sf_getfloat ("dz",&dz)) {
	dz = dt; 
    } else {
	sf_putfloat(out,"d1",dz);
    }
    if (!sf_getfloat ("z0",&z0)) {
	z0 = 0.; 
    } 
    sf_putfloat(out,"o1",z0);

    if (!sf_getint ("nw",&nw)) nw = 4;
    if (!sf_getbool ("slow",&slow)) slow = false;

    fnt = fint1_init (nw, nt);

    time = sf_floatalloc (nt);
    depth = sf_floatalloc (nz);
    vel = sf_floatalloc (nz);

    for (ix = 0; ix < nx; ix++) {
	sf_read (time,sizeof(float),nt,in);
	sf_read (vel,sizeof(float),nz,velocity);

	fint1_set (fnt, time);

	for (iz = 0; iz < nz; iz++) {
	    if (iz == 0) {
		z =  slow? z0*vel[0]/dz: z0/(dz*vel[0]);
	    } else {
		z += slow? vel[iz-1]: 1./vel[iz-1];
	    }
  
	    t = (2.*z*dz-t0)/dt;
	    it = t;
	    t = t - it;

	    if (it < -nw/2) {
		depth[iz] = 0.;
	    } else if (it > nt + nw/2) {
		depth[iz] = time[nt-1];
	    } else {
		depth[iz] = fint1_apply (fnt, it, t, 0);
	    }
	}

	sf_write (depth,sizeof(float),nz,out);
    }

    exit (0);
}

