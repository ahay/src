/* Second-order cell ray tracing with locally parabolic rays.

Takes: < velocity.rsf > rays.rsf

Rays and wavefronts can be displayed with sfplotrays program.
*/

#include <math.h>

#include <rsf.h>

#include "celltrace.h"

int main(int argc, char* argv[])
{
    bool velocity;
    int is, nz, nx, im, nm, order, nshot, ndim, nsr;
    int nt, nr, ir, it;
    float da=0., a0, amax, t;
    float x[2], p[2], dz, dx, z0, x0, **traj, *slow, **s, *a;
    celltrace ct;
    sf_file shots, vel, angles;

    sf_init (argc,argv);
    vel = sf_input("in");

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",&nz))   sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&x0)) x0=0.;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    if(!sf_getint("order",&order)) order=4;

    if (!sf_getint("nt",&nt)) nt=nx*nz;

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
	shots = sf_input("shotfile");
	if (!sf_histint(shots,"n1",&ndim) || 2 != ndim) 
	    sf_error("Must have n1=2 in shotfile");
	if (!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_read(s[0],sizeof(float),ndim*nshot,shots);
	sf_fileclose (shots);
    } else {
	nshot = 1;
	ndim = 2;

	s = sf_floatalloc2 (ndim,nshot);

	if (!sf_getfloat("zshot",s[0]))   s[0][0]=0.; 
	if (!sf_getfloat("yshot",s[0]+1)) s[0][1]=x0 + 0.5*(nx-1)*dx;
	
	sf_warning("Shooting from z=%f, x=%f",s[0][0],s[0][1]);
    }

    if (NULL != sf_getstring("anglefile")) {
	angles = sf_input("anglefile");

	if (!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
    } else {
	angles = NULL;

	if (!sf_getint("nr",&nr)) sf_error("Need nr=");
	if (!sf_getfloat("a0",&a0)) a0 = 0.; 
	if (!sf_getfloat("amax",&amax)) amax=360.;

	/* convert degrees to radians */
	a0 = a0*SF_PI/180.;
	amax = amax*SF_PI/180.;

	/* figure out angle spacing */
	da = (nr > 1)? (amax - a0)/(nr-1) : 0.;
    }

    a = sf_floatalloc(nr);
 
    /* specify output dimensions */
    nsr = nr*nshot;
    fwrite(&nsr,sizeof(int),1,stdout);
	    
    /* get slowness */
    nm = nz*nx;
    slow = sf_floatalloc(nm);

    sf_read(slow,sizeof(float),nm,vel);

    if (vel) { /* convert to slowness */
	for(im = 0; im < nm; im++){
	    slow[im] = 1./slow[im];
	}
    }

    /* initialize ray tracing object */
    ct = celltrace_init (order, nt, nz, nx, dz, dx, z0, x0, slow);
    
    free (slow);

    traj = sf_floatalloc2 (ndim,nt);

    for( is = 0; is < nshot; is++) { /* loop over shots */
	/* initialize angles */
	if (NULL != angles) {
	    sf_read(a,sizeof(float),nr,angles);
	} else {
	    for (ir = 0; ir < nr; ir++) {
		a[ir] = a0+da*ir;
	    }
	}
	for (ir = 0; ir < nr; ir++) { /* loop over rays */
	    /* initialize position */
	    x[0] = s[is][0]; 
	    x[1] = s[is][1];

	    /* initialize direction */
	    p[0] = -cosf(a[ir]);
	    p[1] = sinf(a[ir]);

	    t = cell_trace (ct, x, p, &it, traj);
	    if (it < 0) it = -it; /* keep side-exiting rays */
	    fwrite(&it,sizeof(int),1,stdout);
	    fwrite(traj[0],sizeof(float),(it+1)*ndim,stdout);
	}
    }

    exit (0);
}
