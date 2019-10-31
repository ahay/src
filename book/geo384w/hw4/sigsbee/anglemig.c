/* 2-D angle-domain zero-offset migration. */
#include <rsf.h>

static float get_sample (float **dat, 
			 float t, float y,
			 float t0, float y0, 
			 float dt, float dy, 
			 int nt, int ny) 
/* extract data sample by linear interpolation */
{
    int it, iy;

    y = (y - y0)/dy; iy = floorf (y);
    y -= (float)iy;
    if (iy < 0 || iy >= (ny - 1)) return 0.0;
    t = (t - t0)/dt; it = floorf (t);
    t -= (float)it;
    if (it < 0 || it >= (nt - 1)) return 0.0;

    return (dat[iy][it]*(1.0 - y)*(1.0 - t) +
	    dat[iy][it + 1]*(1.0 - y)*t +
	    dat[iy + 1][it]*y*(1.0 - t) +
	    dat[iy + 1][it + 1]*y*t);
}

int main (int argc, char* argv[]) 
{
    int it, nt, ix, nx, ia, na;
    float dt, vel, da, a0, dx, z, t, x, y, a;
    float **dat, *img;
    sf_file data, imag;

    sf_init (argc, argv);

    data = sf_input ("in");
    imag = sf_output ("out");

    /* get dimensions */
    if (!sf_histint (data, "n1", &nt))   sf_error ("n1");
    if (!sf_histint (data, "n2", &nx))   sf_error ("n2");
    if (!sf_histfloat (data, "d1", &dt)) sf_error ("d1");
    if (!sf_histfloat (data, "d2", &dx)) sf_error ("d2");

    if (!sf_getint("na",&na)) sf_error("Need na="); 
    /* number of angles */
    if (!sf_getfloat("da",&da)) sf_error("Need da="); 
    /* angle increment */
    if (!sf_getfloat("a0",&a0)) sf_error("Need a0="); 
    /* initial angle */

    sf_shiftdim(data, imag, 1);

    sf_putint(imag,"n1",na);
    sf_putfloat(imag,"d1",da);
    sf_putfloat(imag,"o1",a0);
    sf_putstring(imag,"label1","Angle");

    /* degrees to radians */
    a0 *= SF_PI/180.;
    da *= SF_PI/180.;

    if (!sf_getfloat("vel",&vel)) vel=1.5; 
    /* constant velocity */

    dat = sf_floatalloc2(nt,nx);
    sf_floatread (dat[0],nt*nx, data);

    img = sf_floatalloc (na);
 
    /* loop over image location */
    for (ix = 0; ix < nx; ix++) {
	x = ix*dx;
        sf_warning ("CMP %d of %d;", ix, nx);
   
	/* loop over image time */
	for (it = 0; it < nt; it++) {
	    z = it*dt;

	    /* loop over angle */
            for (ia = 0; ia < na; ia++) { 
		a = a0+ia*da;
		
		t = z/cosf(a);           
                /* escape time */
		y = x+0.5*vel*t*sinf(a); 
                /* escape location */

		img[ia] = get_sample (dat,t,y,0.,0.,
				      dt,dx,nt,nx);
	    } /* ia */

            sf_floatwrite (img, na, imag);
        } /* it */
    } /* ix */
 
    exit(0);
}

