#include <math.h>
#include <float.h>

#include <rsf.h>

#include "celltrace.h"
#include "fzero.h"

static celltrace ct;
static float xs[2], xr, t;

static float shooting(float a);

int main(int argc, char* argv[])
{
    bool velocity;
    int is, nz, nx, im, nm, order, nshot, ndim, nt, nr, ir, i1, i2, ia;
    float dt, a, r0, dr, rmax, xmax, a0, amax, da, pi, r1, r2;
    float dz, dx, z0, x0, *angle, *slow, **shot, *time;
    sf_file shots, vel, out;

    sf_init (argc,argv);
    vel = sf_input("in");
    out = sf_output("out");

    /* get 2-D grid parameters */
    if (!sf_histint(vel,"n1",&nz))   sf_error("No n1= in input");
    if (!sf_histint(vel,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&z0)) z0=0.;
    if (!sf_histfloat(vel,"o2",&x0)) x0=0.;
    xmax = x0 + (nx-1)*dx;

    /* additional parameters */
    if(!sf_getbool("vel",&velocity)) velocity=true;
    if(!sf_getint("order",&order)) order=4;

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");

    if (!sf_getint("nr",&nr)) nr=1;
    if (!sf_getfloat("r0",&r0)) r0=x0;
    if (!sf_getfloat("dr",&dr)) dr=dx;
    rmax = r0 + (nr-1)*dr;

    /* sanity check */
    if (dr <= 0.) sf_error ("Need positive dr");
    if (r0 < x0 || r0 > xmax || rmax < x0 || rmax > xmax)
	sf_error ("Range of receivers [%f,%f] "
		  "is outside of the model range [%f,%f]",
		  r0,rmax,x0,xmax);
 
    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
	shots = sf_input("shotfile");
	if (!sf_histint(shots,"n1",&ndim) || 2 != ndim) 
	    sf_error("Must have n1=2 in shotfile");
	if (!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
  
	shot = sf_floatalloc2 (ndim,nshot);
	sf_read(shot[0],sizeof(float),ndim*nshot,shots);
	sf_fileclose (shots);
    } else {
	nshot = 1;
	ndim = 2;

	shot = sf_floatalloc2 (ndim,nshot);

	if (!sf_getfloat("zshot",shot[0]))   shot[0][0]=0.; 
	if (!sf_getfloat("yshot",shot[0]+1)) shot[0][1]=x0 + 0.5*(nx-1)*dx;
	
	sf_warning("Shooting from z=%f, x=%f",shot[0][0],shot[0][1]);
    }

    /* specify output dimensions */
    sf_putint (out,"n1",nr);
    sf_putfloat (out,"d1",dr);
    sf_putfloat (out,"o1",r0);
    sf_putint (out,"n2",nshot);
	    
    /* get velocity */
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

    angle = sf_floatalloc (nr);
    time = sf_floatalloc (nt);

    pi = acos(-1.);

    for( is = 0; is < nshot; is++) { /* loop over shots */
	/* initialize position */
	xs[0] = shot[is][0]; 
	xs[1] = shot[is][1];

	/* initialize directions */
	da = 0.01*pi;
	a0   = atan2f(r0  -xs[1],z0-xs[0]);
	amax = atan2f(rmax-xs[1],z0-xs[0]);

	xr = r0;
	for (a0 -= da; a0 > - pi; a0 -= da) {
	    r1 = shooting(a0);
	    if (r1 <= 0.) break;
	}
	r1 += xr;

	for (amax += da; amax < pi; amax += da) {
	    if (shooting(amax) >= 0.) break;
	}
	
	for (ir=0; ir < nr; ir++) {
	    time[ir] = angle[ir] = FLT_MAX;
	}

	da = (amax-a0)/nr;
	for (ia=0; ia < nr; ia++) {
	    a = a0+(ia+1)*da;
	    r2 = xr+shooting(a);
	    i1 = ceilf((r1-r0)/dr);
	    i2 = floorf((r2-r0)/dr);
	    for (ir=i1; ir <= i2; ir++) {
		if (ir >= nr) break;
		if (ir < 0) continue;
		xr = r0+ir*dr; /* target */
		a = fzero(shooting,a-da,a,r1-xr,r2-xr,0.1*dr,false);
		if (t < time[ir]) { /* select first arrivals */
		    angle[ir] = a;
		    time[ir] = t;
		}
	    }
	    r1 = r2;
      	}
	sf_write(angle,sizeof(float),nr,out);
    }
    
    exit (0);
}

static float shooting(float a)
{
    int it;
    float x[2], p[2];

    x[0]=xs[0];
    x[1]=xs[1];
    
    p[0] = cos(a);
    p[1] = sin(a);

    t = cell_trace (ct, x, p, &it, NULL);

    return (x[1]-xr);
}
