#include <float.h>
#include <math.h>

#include <rsf.h>

#include "int2.h"
#include "interp.h"

int main (int argc, char* argv[])
{
    int id, nd, im, nm, nt, it, nx, ny, hdr[SF_NKEYS];
    int xkey, ykey, interp;
    float *mm, *count, *dd, **xy;
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, clip;
    sf_file in, out, head, fold;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nt)) sf_error("Need n2= in in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("xkey",&xkey)) xkey = sf_segykey("sx");
    if (!sf_getint("ykey",&ykey)) ykey = sf_segykey("sy");

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);
    head = sf_input("head");

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	sf_read (hdr,sizeof(int),SF_NKEYS,head);
	f = hdr[xkey]/1000.; 
	if (f < xmin) xmin=f;
	if (f > xmax) xmax=f;
	xy[id][0] = f;
	f = hdr[ykey]/1000.; 
	if (f < ymin) ymin=f;
	if (f > ymax) ymax=f;
	xy[id][1] = f;
    }

    sf_fileclose (head);

    /* create model */
    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    if (!sf_getint ("ny",&ny)) sf_error("Need ny=");

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nt);

    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(out,"o3",t0);
    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(out,"d3",dt);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    sf_getfloat ("xmax",&xmax);
    sf_getfloat ("ymin",&ymin);
    sf_getfloat ("ymax",&ymax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);
    if (ymax <= ymin) sf_error ("ymax=%f <= ymin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    if (!sf_getfloat("y0",&y0)) y0=ymin; 

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    if (!sf_getfloat("dx",&dx)) {
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }

    if (!sf_getfloat("dy",&dy)) {
	if (1 >= nx) {
	    dy = dx;
	} else {
	    dy = (ymax-ymin)/(ny-1);
	}
    }

    sf_putfloat (out,"d1",dx);
    sf_putfloat (out,"d2",dy);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=1;

    switch (interp) {
	case 1:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, bin_int, 1, nd);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    int2_init (xy, x0,y0,dx,dy,nx,ny, lin_int, 2, nd);
	    sf_warning("Using linear interpolation");
	    break;
	case 3:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);
    count  = sf_floatalloc(nm);

    /* compute fold */
    for (id=0; id<nd; id++) {
	dd[id]=1.;
    }

    int2_lop (true, false,nm,nd,count,dd);
 
    if (NULL != sf_getstring("fold")) {
	fold = sf_output("fold");
	sf_putint(fold,"n1",nx);
	sf_putint(fold,"n2",ny);
	sf_putint(fold,"n3",1);
	sf_putfloat(fold,"o1",x0);
	sf_putfloat(fold,"o2",y0);
	sf_putfloat(fold,"d1",dx);
	sf_putfloat(fold,"d2",dy);
	sf_write (count,sizeof(float),nm,fold);
	sf_fileclose (fold);
    }

    if (!sf_getfloat("clip",&clip)) clip = FLT_EPSILON;

    for (im=0; im<nm; im++) {
	if (clip < count[im]) count[im]=1./fabsf(count[im]);
	else                  count[im]=0.;
    }

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_read (dd,sizeof(float),nd,in);
	int2_lop (true,false,nm,nd,mm,dd);
	for (im=0; im<nm; im++) {
	    mm[im] *= count[im];
	}
	sf_write (mm,sizeof(float),nm,out);
    }

    exit(0);
}


