/* 2-D zero-offset reverse-time migration */
#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int n1, n2;
static float c0, c11, c21, c12, c22;

static void laplacian(float **uin  /* [n2][n1] */, 
		      float **uout /* [n2][n1] */)
/* Laplacian operator, 4th-order finite-difference */
{
    int i1, i2;

#ifdef _OPENMP
#pragma omp parallel for	                  \
    private(i2,i1)		                  \
    shared(n2,n1,uout,uin,c11,c12,c21,c22,c0)
#endif	    
    for (i2=2; i2 < n2-2; i2++) {
	for (i1=2; i1 < n1-2; i1++) {
	    uout[i2][i1] = 
		c11*(uin[i2][i1-1]+uin[i2][i1+1]) +
		c12*(uin[i2][i1-2]+uin[i2][i1+2]) +
		c21*(uin[i2-1][i1]+uin[i2+1][i1]) +
		c22*(uin[i2-2][i1]+uin[i2+2][i1]) +
		c0*uin[i2][i1];
	}
    }
}


int main(int argc, char* argv[])
{
    int it,i1,i2;        /* index variables */
    int nt,n12,n0,nx, jt;
    float dt,dx,dz, dt2,d1,d2;

    float  **vv, **dd;          
    float  **u0,**u1,u2,**ud; /* tmp arrays */

    sf_file data, imag, modl, wave; /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

    /* initialize OpenMP support */ 
    omp_init();

    /* setup I/O files */
    modl = sf_input ("in");   /* velocity model */
    imag = sf_output("out");  /* output image */
 
    data = sf_input ("data"); /* seismic data */
    wave = sf_output("wave"); /* wavefield */
    
    /* Dimensions */
    if (!sf_histint(modl,"n1",&n1)) sf_error("n1");
    if (!sf_histint(modl,"n2",&n2)) sf_error("n2");

    if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
    if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");

    if (!sf_histint  (data,"n1",&nt)) sf_error("n1");
    if (!sf_histfloat(data,"d1",&dt)) sf_error("d1");

    if (!sf_histint(data,"n2",&nx) || nx != n2) 
	sf_error("Need n2=%d in data",n2);

    n12 = n1*n2;

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
    dd = sf_floatalloc2(nt,n2);
    sf_floatread(dd[0],nt*n2,data);

    vv = sf_floatalloc2(n1,n2);
    sf_floatread(vv[0],n12,modl);

    /* allocate temporary arrays */
    u0=sf_floatalloc2(n1,n2);
    u1=sf_floatalloc2(n1,n2);
    ud=sf_floatalloc2(n1,n2);
    
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    u0[i2][i1]=0.0;
	    u1[i2][i1]=0.0;
	    ud[i2][i1]=0.0;
	    vv[i2][i1] *= vv[i2][i1]*dt2;
	}
    }

    /* Time loop */
    for (it=nt-1; it >= 0; it--) {
	sf_warning("%d;",it);

	laplacian(u1,ud);

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(i2,i1,u2)		    \
    shared(ud,vv,it,u1,u0,dd,n0)
#endif
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1<n1; i1++) {
		/* scale by velocity */
		ud[i2][i1] *= vv[i2][i1];

		/* time step */
		u2 = 
		    2*u1[i2][i1] 
		    - u0[i2][i1] 
		    + ud[i2][i1]; 
		
		u0[i2][i1] = u1[i2][i1];
		u1[i2][i1] = u2;
	    }

	    /* inject data */
	    u1[i2][n0] += dd[i2][it];
	}

	if (0 == it%jt) 
	    sf_floatwrite(u1[0],n12,wave);       
    }
    sf_warning(".");	

    /* output image */
    sf_floatwrite(u1[0],n12,imag);
    
    exit (0);
}
