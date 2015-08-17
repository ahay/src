/* 2-D finite-difference born modeling */
#include <rsf.h>

static int n1, n2;
/* Laplacian coefficients */
static float c0, c11, c21, c12, c22;

static void laplacian(float **uin  /* [n2][n1] */, 
		      float **uout /* [n2][n1] */)
/* Laplacian operator, 4th-order finite-difference */
{
    int i1, i2;

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

static void deriv2(float *ww, int nt)
/* second time derivative */
{
	int it;
	float *temp;
	temp=sf_floatalloc(nt);
	temp[0]=ww[0];
	for(it=1; it<nt; it++)
		temp[it]=ww[it]-ww[it-1];
	for(it=0; it<nt-1; it++)
		ww[it]=temp[it+1]-temp[it];
	ww[nt-1]=temp[nt-1]-temp[nt-2];
}

int main(int argc, char* argv[])
/* main function */
{
	bool born; /* born modeling or forward modeling */

    int it,i1,i2;        /* index variables */
    int nt,n12,ft,jt;
    float dt,d1,d2,dt2;

	int nr, r0, rz, rx; /* receiver setup */

    float  *ww,**vv,**rr,**ref,**dd;  /* variables */
    float **u0,**u1,**u2,**ud,***wave,**vv2; /* temporary variables */

    sf_file Fw,Fv,Fr,Fo,Fd,Fdetv;  /* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

	/* born modeling or forward modeling */
	if(!sf_getbool("born", &born)) born=true;

    /* setup I/O files */
    Fr = sf_input ("in");    /* source position */
    Fd = sf_output("out");   /* output shot record */

	Fo = sf_output("snapshot"); /* output wavefield */
    Fw = sf_input ("wav");   /* source wavelet */
    Fv = sf_input ("v");     /* velocity */
	if (born) Fdetv = sf_input("detv");
	/* velocity perturbation/reflectivity */

    /* Read/Write axes */
    if (!sf_histint(Fr,"n1",&n1))   sf_error("No n1= in inp");
    if (!sf_histint(Fr,"n2",&n2))   sf_error("No n2= in inp");
    if (!sf_histfloat(Fr,"d1",&d1)) sf_error("No d1= in inp");
    if (!sf_histfloat(Fr,"d2",&d2)) sf_error("No d2= in inp");

    if (!sf_histint(Fw,"n1",&nt))   sf_error("No n1= in wav");
    if (!sf_histfloat(Fw,"d1",&dt)) sf_error("No d1= in wav");

    n12 = n1*n2;

	/* trace number of shot record */
	if (!sf_getint("nr", &nr)) nr=n2;
	/* starting position of shot record */
	if (!sf_getint("r0", &r0)) r0=0;
	/* depth of shot record */
	if (!sf_getint("rz", &rz)) rz=0;
    /* first recorded time */
    if (!sf_getint("ft",&ft)) ft=0; 
    /* time interval */
    if (!sf_getint("jt",&jt)) jt=1; 

	/* set the dimension of output data file */
	sf_putint(Fd, "n1", nt);
	sf_putfloat(Fd, "d1", dt);
	sf_putfloat(Fd, "o1", 0.);
	sf_putstring(Fd, "label1", "Time");
	sf_putstring(Fd, "unit1", "s");
	
	sf_putint(Fd, "n2", nr);
	sf_putfloat(Fd, "d2", d2);
	sf_putfloat(Fd, "o2", r0*d2);

	/* set the dimension of output wavefield file */
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

    /* read wavelet, velocity, source position & reflectivity*/
    ww=sf_floatalloc(nt);     sf_floatread(ww   ,nt ,Fw);
    vv=sf_floatalloc2(n1,n2); sf_floatread(vv[0],n12,Fv);
    rr=sf_floatalloc2(n1,n2); sf_floatread(rr[0],n12,Fr);
    if(born) {ref=sf_floatalloc2(n1,n2); sf_floatread(ref[0],n12,Fdetv);}

	/* allocate wavefield and data arrays */
	wave=sf_floatalloc3(n1, n2, nt);
	dd=sf_floatalloc2(nt, nr);
	vv2=sf_floatalloc2(n1, n2);

    /* allocate temporary arrays */
    u0=sf_floatalloc2(n1,n2);
    u1=sf_floatalloc2(n1,n2);
    u2=sf_floatalloc2(n1,n2);
    ud=sf_floatalloc2(n1,n2);
   
	/* initialization */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    u0[i2][i1]=0.0;
	    u1[i2][i1]=0.0;
	    u2[i2][i1]=0.0;
	    ud[i2][i1]=0.0;
		vv2[i2][i1]=vv[i2][i1]*vv[i2][i1]*dt2;
	}
    }

	/* second derivative */
	if(born) deriv2(ww, nt);

    /* Time loop 1 (U_0) */
    for (it=0; it<nt; it++) {
		sf_warning("Loop 1, it=%d;", it);
	
		/* wavefield storage */
		for(i2=0; i2<n2; i2++)
			for(i1=0; i1<n1; i1++)
				wave[it][i2][i1]=u1[i2][i1];
		
		/* call function */
		laplacian(u1,ud);
		
		for (i2=0; i2<n2; i2++) {
			for (i1=0; i1<n1; i1++) {
				
				/* scale by velocity */
				ud[i2][i1] *= vv2[i2][i1];
				
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
	} //end of it

	/* if forward modeling, solve wave equation once */ 
	if(!born){
		/* write wavefield to output */
		for(it=0; it<nt; it++){
			for(i2=0; i2<nr; i2++){
				rx=r0+i2;
				dd[i2][it]=wave[it][rx][rz];
			}
		
			if (it >= ft && 0 == (it-ft)%jt) {
				sf_floatwrite(wave[it][0],n12,Fo);
			}
		}

		sf_floatwrite(dd[0], nt*nr, Fd);
		
		/* if forward modeling, only solve wave equation once */
		exit(0);
	}
    
	/* second initialization */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    u0[i2][i1]=0.0;
	    u1[i2][i1]=0.0;
	    u2[i2][i1]=0.0;
	    ud[i2][i1]=0.0;
	}
    }

    /* Time loop 2 (det U) */
    for (it=0; it<nt; it++) {
		sf_warning("Loop 2, it=%d;", it);
	
		/* write shot record to output */
		for(i2=0; i2<nr; i2++){
			rx=r0+i2;
			dd[i2][it]=u1[rx][rz];
		}
		
		/* write wavefield to output */
		if (it >= ft && 0 == (it-ft)%jt) {
			sf_floatwrite(u1[0],n12,Fo);
		}

		/* call function */
		laplacian(u1,ud);
		
		for (i2=0; i2<n2; i2++) {
			for (i1=0; i1<n1; i1++) {
				
				/* scale by velocity */
				ud[i2][i1] *= vv2[i2][i1];
				
				/* inject source term */
				ud[i2][i1] += wave[it][i2][i1] 
					* ref[i2][i1] * 2./vv[i2][i1];
		
				/* time step */
				u2[i2][i1] = 
					2*u1[i2][i1] 
					- u0[i2][i1] 
					+ ud[i2][i1]; 
				u0[i2][i1] = u1[i2][i1];
				u1[i2][i1] = u2[i2][i1];
			}
		}
	} //end of it

	sf_floatwrite(dd[0], nt*nr, Fd);	

    exit (0);
}
