/* 2-D Poststack Kirchhoff time migration. */
#include <rsf.h>

static void doubint(int nt, float *trace)
/* causal and anticausal integratio */
{
    int it;
    float tt;

    tt = trace[0];
    for (it=1; it < nt; it++) {
	tt += trace[it];
	trace[it] = tt;
    }
    tt = trace[nt-1];
    for (it=nt-2; it >=0; it--) {
	tt += trace[it];
	trace[it] = tt;
    }
}
    
static float pick(float ti, float deltat, 
		  const float *trace,
		  int nt, float dt, float t0)
/* pick a traveltime sample from a trace */
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return 0.0;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return 0.0;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return 0.0;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    return imp*(
	2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	(1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	(1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
}

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

    for (iy=0; iy < nx; iy++) { 
	sf_floatread (trace,nt,inp);
	doubint(nt,trace);
        
	for (ix=0; ix < nx; ix++) { 
	    x = (ix-iy)*dx;

	    for (iz=0; iz < nz; iz++) {
		z = z0 + iz*dz;  
		vi = v[ix][iz];

		/* hypot(a,b) = sqrt(a*a+b*b) */
		ti = hypotf(z,2.0*x/vi);		

		/* tx = |dt/dx| */
		tx = 4.0*fabsf(x)/(vi*vi*(ti+dt));

		out[ix][iz] += 
		    pick(ti,tx*dx*aal,trace,nt,dt,t0);
	    } 
	} 
    } 
    
    sf_floatwrite(out[0],nz*nx,mig);        

    exit(0);
}
