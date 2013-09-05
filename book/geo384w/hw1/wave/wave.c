/* 2-D finite-difference acoustic wave propagation */
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
    int nt,n12,ft,jt;
    float dt,d1,d2,dt2;

    float  *ww,**vv,**rr;       
    float **u0,**u1,**u2,**ud; 

    sf_file Fw,Fv,Fr,Fo;  /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

    /* initialize OpenMP support */ 
    omp_init();

    /* setup I/O files */
    Fr = sf_input ("in");    /* source position */
    Fo = sf_output("out");   /* output wavefield */

    Fw = sf_input ("wav");   /* source wavelet */
    Fv = sf_input ("v");     /* velocity */

    /* Read/Write axes */
    if (!sf_histint(Fr,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fr,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(Fr,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(Fr,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_histint(Fw,"n1",&nt))   sf_error("No n1= in wav");
    if (!sf_histfloat(Fw,"d1",&dt)) sf_error("No d1= in wav");

    n12 = n1*n2;

    if (!sf_getint("ft",&ft)) ft=0; 
    /* first recorded time */
    if (!sf_getint("jt",&jt)) jt=1; 
    /* time interval */

    sf_putint(Fo,"n3",(nt-ft)/jt);
    sf_putfloat(Fo,"d3",jt*dt);
    sf_putfloat(Fo,"o3",ft*dt);

    dt2 =    dt*dt;

    /* set Laplacian coefficients */
    d1 = 1.0/(d1*d1);
    d2 = 1.0/(d2*d2);

    c11 = 4.0*d1/3.0;
    c12=  -d1/12.0;
    c21 = 4.0*d2/3.0;
    c22=  -d2/12.0;
    c0  = -2.0 * (c11+c12+c21+c22);

    /* read wavelet, velocity & source position */
    ww=sf_floatalloc(nt);     sf_floatread(ww   ,nt ,Fw);
    vv=sf_floatalloc2(n1,n2); sf_floatread(vv[0],n12,Fv);
    rr=sf_floatalloc2(n1,n2); sf_floatread(rr[0],n12,Fr);

    /* allocate temporary arrays */
    u0=sf_floatalloc2(n1,n2);
    u1=sf_floatalloc2(n1,n2);
    u2=sf_floatalloc2(n1,n2);
    ud=sf_floatalloc2(n1,n2);
    
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    u0[i2][i1]=0.0;
	    u1[i2][i1]=0.0;
	    u2[i2][i1]=0.0;
	    ud[i2][i1]=0.0;
	    vv[i2][i1] *= vv[i2][i1]*dt2;
	}
    }

    /* Time loop */
    for (it=0; it<nt; it++) {
	laplacian(u1,ud);

#ifdef _OPENMP
#pragma omp parallel for	    \
    private(i2,i1)		    \
    shared(ud,vv,ww,it,rr,u2,u1,u0)
#endif	
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1<n1; i1++) {
		/* scale by velocity */
		ud[i2][i1] *= vv[i2][i1];
		/* inject wavelet */
		ud[i2][i1] += ww[it] * rr[i2][i1];
		/* time step */
		u2[i2][i1] = 
		    2*u1[i2][i1] 
		    - u0[i2][i1] 
		    + ud[i2][i1]; 
		u0[i2][i1] = u1[i2][i1];
		u1[i2][i1] = u2[i2][i1];
	    }
	}		
	
	/* write wavefield to output */
	if (it >= ft && 0 == (it-ft)%jt) {
	    sf_warning("%d;",it+1);
	    sf_floatwrite(u1[0],n12,Fo);
	}
    }
    sf_warning(".");
    
    exit (0);
}
