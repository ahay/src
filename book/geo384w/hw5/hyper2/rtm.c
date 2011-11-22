/* 2-D zero-offset reverse-time migration */
#include <stdio.h>

#include <rsf.h>

static int nz, nx;
static float c0, c11, c21, c12, c22;

static void laplacian(float **uin  /* [nx][nz] */, 
		      float **uout /* [nx][nz] */)
/* Laplacian operator, 4th-order finite-difference */
{
    int iz, ix;
    
    for (ix=2; ix < nx-2; ix++) {
	for (iz=2; iz < nz-2; iz++) {
	    uout[ix][iz] = 
		c11*(uin[ix][iz-1]+uin[ix][iz+1]) +
		c12*(uin[ix][iz-2]+uin[ix][iz+2]) +
		c21*(uin[ix-1][iz]+uin[ix+1][iz]) +
		c22*(uin[ix-2][iz]+uin[ix+2][iz]) +
		c0*uin[ix][iz];
	}
    }
}

int main(int argc, char* argv[])
{
    int it,ix,iz;        /* index variables */
    int nt,n0,n2, jt;
    float dt,dx,dz, dt2,d1,d2;

    float  **vv, **dd;          
    float  **u0,**u1,u2,**ud; /* tmp arrays */

    sf_file data, imag, modl, wave; /* I/O files */

    sf_init(argc,argv);

    /* setup I/O files */
    modl = sf_input ("in");   /* velocity model */
    imag = sf_output("out");  /* output image */
 
    data = sf_input ("data"); /* seismic data */
    wave = sf_output("wave"); /* wavefield */
    
    /* Dimensions */
    if (!sf_histint(modl,"n1",&nz)) sf_error("n1");
    if (!sf_histint(modl,"n2",&nx)) sf_error("n2");

    if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
    if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");

    if (!sf_histint  (data,"n1",&nt)) sf_error("n1");
    if (!sf_histfloat(data,"d1",&dt)) sf_error("d1");

    if (!sf_histint(data,"n2",&n2) || n2 != nx) 
	sf_error("Need n2=%d in data",nx);

    if (!sf_getint("n0",&n0)) n0=0;
    /* surface */
    if (!sf_getint("jt",&jt)) jt=1; 
    /* time interval */

    sf_putint(wave,"n3",1+(nt-1)/jt);
    sf_putfloat(wave,"d3",-jt*dt);
    sf_putfloat(wave,"o3",(nt-1)*dt);

    dt2 = dt*dt;

    /* set Laplacian coefficients */
    d1 = 1.0/(dz*dz);
    d2 = 1.0/(dx*dx);

    c11 = 4.0*d1/3.0;
    c12=  -d1/12.0;
    c21 = 4.0*d2/3.0;
    c22=  -d2/12.0;
    c0  = -2.0 * (c11+c12+c21+c22);

    /* read data and velocity */
    dd = sf_floatalloc2(nt,nx);
    sf_floatread(dd[0],nt*nx,data);

    vv = sf_floatalloc2(nz,nx);
    sf_floatread(vv[0],nz*nx,modl);

    /* allocate temporary arrays */
    u0=sf_floatalloc2(nz,nx);
    u1=sf_floatalloc2(nz,nx);
    ud=sf_floatalloc2(nz,nx);
    
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    u0[ix][iz]=0.0;
	    u1[ix][iz]=0.0;
	    ud[ix][iz]=0.0;
	    vv[ix][iz] *= vv[ix][iz]*dt2;
	}
    }

    /* Time loop */
    for (it=nt-1; it >= 0; it--) {
	sf_warning("%d;",it);

	laplacian(u1,ud);

	for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
		/* scale by velocity */
		ud[ix][iz] *= vv[ix][iz];

		/* time step */
		u2 = 
		    2*u1[ix][iz] 
		    - u0[ix][iz] 
		    + ud[ix][iz]; 
		
		u0[ix][iz] = u1[ix][iz];
		u1[ix][iz] = u2;
	    }

	    /* inject data */
	    u1[ix][n0] += dd[ix][it];
	}

	if (0 == it%jt) 
	    sf_floatwrite(u1[0],nx*nz,wave);       
    }
    sf_warning(".");	

    /* output image */
    sf_floatwrite(u1[0],nx*nz,imag);
    
    exit (0);
}
