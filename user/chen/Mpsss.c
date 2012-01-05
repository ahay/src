/* EFD phase shift wave extrapolation. */

#include <rsf.h>
#include <fftw3.h>

#include "psss.h"


int main(int argc, char* argv[])
{
	int it,ix,iz,jt,njt;
	float dt,dx,dz;
	int nx,nz,nt,nz0,nw;
	float ox,**v0;
	double *v1,*v2;
	float **u2,**u3,**vel;  
	fftw_plan plan;

	sf_file modl,data,wave,imag;

	sf_init(argc,argv);

	modl = sf_input("in"); 		/* velocity model m(z,x) */
    imag = sf_output ("out"); 	/* image data: u(0,z,x) */
    data = sf_input("data");  	/* seismic data, zero offset u(t,0,x) */
    wave = sf_output("wave"); 	/* wavefield: u(t,z,x) */

    if (!sf_histint(modl,"n1",&nz0)) sf_error("n1");
    if (!sf_histint(modl,"n2",&nx)) sf_error("n2");
    if (!sf_histfloat(modl,"o2",&ox)) sf_error("o2");
    if (!sf_histint(data,"n1",&nt)) sf_error("n1");
    if (!sf_histfloat(data,"d1",&dt)) sf_error("d2");

    if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
    if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");

    if (!sf_getint("jt",&jt)) jt=40; 
    /* time step for wave data */
    if (!sf_getint("nz",&nz)) nz=nz0; 
    /* depth number */

    sf_putint(imag,"n1",nz);
    sf_putint(imag,"n2",nx);
    sf_putfloat(imag,"o2",ox);

	njt=1+(nt-1)/jt;
    sf_putint(wave,"n1",njt);
    sf_putfloat(wave,"d1",jt*dt);
    sf_putfloat(wave,"o1",0);
    sf_putint(wave,"n2",nx);
    sf_putfloat(wave,"d2",dx);
    sf_putfloat(wave,"o2",ox);
    sf_putint(wave,"n3",nz);
    sf_putfloat(wave,"d3",dz);
    sf_putfloat(wave,"o3",0);

	nw=nt/2+1;

    /* read data and velocity */
    vel = sf_floatalloc2(nz,nx);
    sf_floatread(vel[0],nz*nx,modl);

    v0 = sf_floatalloc2(nt,nx);	// U_z(t,x)
    v1 = (double*)malloc(nt*nx*sizeof(double));	// U_z(t,x)
    v2 = (double*)malloc(nw*2*nx*sizeof(double));	// U_z(t,x)
	
    u2 = sf_floatalloc2(njt,nx);	// u_z(x,t)
    u3 = sf_floatalloc2(nz,nx);		// u(x,z,0)

	plan = fftw_plan_dft_r2c_2d(nx, nt, v1, 
		(fftw_complex*)v2,  FFTW_ESTIMATE);

    sf_floatread(v0[0],nt*nx,data);

	for(ix=0;ix<nx;ix++)	
	{
		for(it=0;it<njt;it++)  { u2[ix][it]=v0[ix][it*jt]; v1[ix*nt+it]=v0[ix][it];}
        u3[ix][0]=u2[ix][0];  // imag iz =0
	}
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	sf_floatwrite(u2[0],njt*nx,wave); // wave slice iz=0


	plan = fftw_plan_dft_c2r_2d(nx,nt, (fftw_complex*)v2,
		v1,  FFTW_ESTIMATE);

	sf_psss_init(nw,nx,nz,
		1.0/(dt*nt),dx,dz,vel);


	for(iz=1;iz<nz;iz++)
	{
//		sf_psss_step(iz,v2);

		fftw_execute(plan);
		for(ix=0;ix<nx;ix++)	
		{
//			for(it=0;it<njt;it++)   u2[ix][it]=v1[ix][it*jt]/(nt*nx);
			u3[ix][iz]=u2[ix][0];  // imag iz =0
		}
		sf_floatwrite(u2[0],njt*nx,wave); 

		sf_warning("%d;",iz);
	}

	sf_floatwrite(u3[0],nz*nx,imag);
 
	fftw_destroy_plan(plan);
	sf_psss_exit();

	return 0;
}

