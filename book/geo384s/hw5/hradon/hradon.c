/* Hyperbolic Radon transform */
#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_map4 map;
    int it,nt, ix,nx, ip,np, i3,n3;
    bool adj;
    float t0,dt,t, x0,dx,x, p0,dp,p;
    float **cmp, **rad, *trace;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1=");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1=");

    n3 = sf_leftsize(in,2);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	if (!sf_histint(in,"n2",&nx))   sf_error("No n2=");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2=");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");

	if (!sf_getint("np",&np)) sf_error("Need np=");
	if (!sf_getfloat("op",&p0)) sf_error("need p0=");
	if (!sf_getfloat("dp",&dp)) sf_error("need dp=");

	sf_putint(out,"n2",np);
	sf_putfloat(out,"o2",p0);
	sf_putfloat(out,"d2",dp);	
    } else {
	if (!sf_histint(in,"n2",&np))   sf_error("No n2=");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2=");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2=");

	if (!sf_getint("nx",&nx)) sf_error("Need nx=");
	if (!sf_getfloat("ox",&x0)) sf_error("need x0=");
	if (!sf_getfloat("dx",&dx)) sf_error("need dx=");

	sf_putint(out,"n2",nx);
	sf_putfloat(out,"o2",x0);
	sf_putfloat(out,"d2",dx);
    }

    trace = sf_floatalloc(nt);
    cmp = sf_floatalloc2(nt,nx);
    rad = sf_floatalloc2(nt,np);

    /* initialize half-order differentiation */
    sf_halfint_init (true,nt,1.0f-1.0f/nt);

    /* initialize spline interpolation */
    map = sf_stretch4_init (nt, t0, dt, nt, 0.01);
    
    for (i3=0; i3 < n3; i3++) { 
	if( adj) {
	    for (ix=0; ix < nx; ix++) { 
		sf_floatread(trace,nt,in);
		sf_halfint_lop(true,false,nt,nt,cmp[ix],trace);
	    }
	} else {	    
	    sf_floatread(rad[0],nt*np,in);	    
	}

	sf_adjnull(adj,false,nt*np,nt*nx,rad[0],cmp[0]);
	
	for (ip=0; ip < np; ip++) { 
	    p = p0 + ip*dp;
	    for (ix=0; ix < nx; ix++) { 
		x = x0 + ix*dx;
		
		for (it=0; it < nt; it++) {		
		    t = t0 + it*dt;
		    trace[it] = hypotf(t,p*x);
                    /* hypot(a,b)=sqrt(a*a+b*b) */
		}

		sf_stretch4_define(map,trace);
		
		if (adj) {
		    sf_stretch4_apply_adj(true,map,rad[ip],cmp[ix]);
		} else {
		    sf_stretch4_apply    (true,map,rad[ip],cmp[ix]);
		}
	    }
	}

	if( adj) {
	    sf_floatwrite(rad[0],nt*np,out);
	} else {
	    for (ix=0; ix < nx; ix++) {
		sf_halfint_lop(false,false,nt,nt,cmp[ix],trace);
		sf_floatwrite(trace,nt,out);
	    }
	}
    }

    exit(0);
}
