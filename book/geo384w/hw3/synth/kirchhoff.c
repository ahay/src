/* 2-D Post-stack Kirchhoff time migration. */
#include <rsf.h>

#include "aal.h" /* antialising routines */

int main(int argc, char* argv[])
{
    int nt, nx, nz, ix, iz, iy, i;
    float *trace, **out, **v;
    float x,z, dx, ti, tx, t0,dt, z0,dz, vi,aal;
    sf_file inp, mig, vel;

    sf_init (argc,argv);
    inp = sf_input("in");
    vel = sf_input("vel");
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2=");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2=");

    if (!sf_getint("nz",&nz)) nz=nt;
    if (!sf_getfloat("dz",&dz)) dz=dt;
    if (!sf_getfloat("z0",&z0)) z0=t0;

    if (!sf_getfloat("antialias",&aal)) aal=1.0;
    /* antialiasing */

    v = sf_floatalloc2(nz,nx);
    sf_floatread(v[0],nz*nx,vel);

    trace = sf_floatalloc(nt);
    out = sf_floatalloc2(nz,nx);

    for (i=0; i < nz*nx; i++) {
	out[0][i] = 0.;
    }

    /* loop over input traces */
    for (iy=0; iy < nx; iy++) { 
	sf_floatread (trace,nt,inp);
	aal_doubint(nt,trace);
        
	/* loop over output traces */
	for (ix=0; ix < nx; ix++) { 
	    x = (ix-iy)*dx;

	    /* loop over output time */
	    for (iz=0; iz < nz; iz++) {
		z = z0 + iz*dz;  
		vi = v[ix][iz];

		/* hypot(a,b) = sqrt(a*a+b*b) */
		ti = hypotf(z,2.0*x/vi);		

		/* tx = |dt/dx| */
		tx = 4.0*fabsf(x)/(vi*vi*(ti+dt));

		out[ix][iz] += 
		    aal_pick(ti,tx*dx*aal,trace,nt,dt,t0);
	    } 
	} 
    } 
    
    sf_floatwrite(out[0],nz*nx,mig);        

    exit(0);
}
