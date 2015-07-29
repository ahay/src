/* Data regulatization by inverse interpolation. */
#include <rsf.h>

static void lint (float x, int n, float* w) 
/*< linear interpolation>*/
{
    w[0] = 1.0f - x;
    w[1] = x;
}

static void regrid( int dim         /* dimensions */, 
		    const int* nold /* old size [dim] */, 
		    const int* nnew /* new size [dim] */, 
		    sf_filter aa    /* filter */) 
/*< change data size >*/
{
    int i, h0, h1, h, ii[SF_MAX_DIM];

    for (i=0; i < dim; i++) {
	ii[i] = nold[i]/2-1;
    }
  
    h0 = sf_cart2line( dim, nold, ii); 
    h1 = sf_cart2line( dim, nnew, ii); 
    for (i=0; i < aa->nh; i++) { 
	h = aa->lag[i] + h0;
	sf_line2cart( dim, nold, h, ii);
	aa->lag[i] = sf_cart2line( dim, nnew, ii) - h1;
    }
}

int main (int argc, char* argv[])
{
    int id, nd, nm, nx, ny, na, ia, niter, three, n[2], m[2];
    float *mm, *dd, **xy;
    float x0, y0, dx, dy, a0, eps;
    char *lagfile;
    bool prec;
    sf_filter aa;
    sf_file in, out, flt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* read data */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");
    if (!sf_histint(in,"n1",&three) || 3 != three) 
	sf_error("Need n1=3 in in");
    if (!sf_histint(in,"n2",&nd)) sf_error("Need n2=");

    xy = sf_floatalloc2(3,nd);
    sf_floatread(xy[0],3*nd,in);

    dd = sf_floatalloc(nd);
    for (id=0; id < nd; id++) dd[id] = xy[id][2];

    /* create model */

    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    if (!sf_getint ("ny",&ny)) sf_error("Need ny=");
    /* Number of bins */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);

    if (!sf_getfloat("x0",&x0)) sf_error("Need x0=");
    if (!sf_getfloat("y0",&y0)) sf_error("Need y0=");
    /* grid origin */

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    if (!sf_getfloat("dy",&dy)) sf_error("Need dy=");
    /* grid sampling */

    sf_putfloat (out,"d1",dx);
    sf_putfloat (out,"d2",dy);
 
    nm = nx*ny;
    mm = sf_floatalloc(nm);

    sf_int2_init (xy, x0,y0, dx,dy, nx,ny, lint, 2, nd);

    /* read filter */
    flt = sf_input("filt");

    if (NULL == (lagfile = sf_histstring(flt,"lag"))) 
	sf_error("Need lag= in filt");
    lag = sf_input(lagfile);

    n[0] = nx;
    n[1] = ny;
    if (!sf_histints (lag,"n",m,2)) {
	m[0] = nx;
	m[1] = ny;
    }

    if (!sf_histint(flt,"n1",&na)) sf_error("No n1= in filt");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(flt,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,flt);

    for( ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }
	
    sf_intread(aa->lag,na,lag);
    regrid (2, m, n, aa);

    if (!sf_getbool("prec",&prec)) prec=false;
    /* if use preconditioning */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.01; 
    /* regularization parameter */

    if (prec) {
	sf_polydiv_init (nm, aa);
	sf_solver_prec (sf_int2_lop, sf_cgstep, 
			sf_polydiv_lop,
			nm, nm, nd, 
			mm, dd, niter, eps, "end");
    } else {
	sf_igrad2_init (nx, ny);
	sf_solver_reg (sf_int2_lop, sf_cgstep, 
		       sf_igrad2_lop,
		       2*nm, nm, nd, 
		       mm, dd, niter, eps, "end");
    }
    
    sf_floatwrite (mm,nm,out);
    exit(0);
}

