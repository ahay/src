/* Hyperbolic Radon transform */
#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_map4 map;
    int it,nt, ix,nx, iv,nv, i3,n3;
    bool adj;
    float t0,dt,t, x0,dx,x, v0,dv,v;
    float **cmp, **rad, *trace;
    sf_file in, out;

    /* initialize Madagascar */
    sf_init(argc,argv);

    /* input and output files */
    in = sf_input("in");
    out = sf_output("out");

    /* Time axis parameters from the input */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1=");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1=");

    /* number of CMPS */
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	if (!sf_histint(in,"n2",&nx))   sf_error("No n2=");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2=");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");

	if (!sf_getint("nv",&nv))   sf_error("Need nv=");
	if (!sf_getfloat("ov",&v0)) sf_error("need ov=");
	if (!sf_getfloat("dv",&dv)) sf_error("need dv=");

	sf_putint(out,"n2",nv);
	sf_putfloat(out,"o2",v0);
	sf_putfloat(out,"d2",dv);
	sf_putstring(out,"label2","Velocity");
	sf_putstring(out,"unit2","km/s");
    } else {
	if (!sf_histint(in,"n2",&nv))   sf_error("No n2=");
	if (!sf_histfloat(in,"o2",&v0)) sf_error("No o2=");
	if (!sf_histfloat(in,"d2",&dv)) sf_error("No d2=");

	if (!sf_getint("nx",&nx))   sf_error("Need nx=");
	if (!sf_getfloat("ox",&x0)) sf_error("need ox=");
	if (!sf_getfloat("dx",&dx)) sf_error("need dx=");

	sf_putint(out,"n2",nx);
	sf_putfloat(out,"o2",x0);
	sf_putfloat(out,"d2",dx);
	sf_putstring(out,"label2","Offset");
	sf_putstring(out,"unit2","km");
    }

    /* allocate storage */
    trace = sf_floatalloc(nt);
    rad = sf_floatalloc2(nt,nv);
    cmp = sf_floatalloc2(nt,nx);


    /* initialize half-order differentiation */
    sf_halfint_init (true,nt,1.0f-1.0f/nt);

    /* initialize spline interpolation */
    map = sf_stretch4_init (nt, t0, dt, nt, 0.01);
    
    for (i3=0; i3 < n3; i3++) { 
	if( adj) {
/*
	    for (ix=0; ix < nx; ix++) { 
		sf_floatread(trace,nt,in);
		sf_halfint_lop(XXXX,false,nt,nt,cmp[ix],trace);
	    }
*/
	    sf_floatread(cmp[0],nt*nx,in);
	} else {
	    sf_floatread(rad[0],nt*nv,in);
	}

	/* zero output */
	sf_adjnull(adj,false,nt*nv,nt*nx,rad[0],cmp[0]);
	
	for (iv=0; iv < nv; iv++) {	    
	    v = v0 + iv*dv;
	    for (ix=0; ix < nx; ix++) { 
		x = (x0 + ix*dx)/v;
		
		for (it=0; it < nt; it++) {		
		    t = t0 + it*dt;
		    trace[it] = hypotf(t,x);
                    /* hypot(a,b)=sqrt(a*a+b*b) */
		}

		/* define mapping */
		sf_stretch4_define(map,trace);
		
		if (adj) { /* rad <- cmp */
		    sf_stretch4_apply_adj(true,map,rad[iv],cmp[ix]);
		} else {   /* rad -> cmp */
		    sf_stretch4_apply    (true,map,rad[iv],cmp[ix]);
		}
	    }
	}

	if (adj) {
	    sf_floatwrite(rad[0],nt*nv,out);
	} else {
	    sf_floatwrite(cmp[0],nt*nx,out);
/*	    
	    for (ix=0; ix < nx; ix++) {
	        sf_halfint_lop(XXXX,false,nt,nt,cmp[ix],trace);
		sf_floatwrite(trace,nt,out);
	    } 
*/
	}
    }

    exit(0);
}
