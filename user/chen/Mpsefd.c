/* EFD phase shift wave extrapolation. */

#include <rsf.h>
#include <fftw3.h>

#include "psefd.h"


int main(int argc, char* argv[])
{
	int iw,it,ix,iz,jt,njt;
	float dt,dx,dz;
	int nw,nx,nz,nt,nz0;
	float *v1; // fft buffer
	float ox;
	float **u1,**u2,**u3,**vel;  
	fftwf_plan plan;

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

    v1 = sf_floatalloc(nw*2);

    u1 = sf_floatalloc2(nx*2,nw);	// U_z(x,w)
    u2 = sf_floatalloc2(njt,nx);	// u_z(x,t)
    u3 = sf_floatalloc2(nz,nx);		// u(x,z,0)

	plan = fftwf_plan_dft_r2c_1d(nt, v1, 
		(fftwf_complex*)v1,  FFTW_ESTIMATE);

	for(ix=0;ix<nx;ix++)	
	{
    	sf_floatread(v1,nt,data);
		for(it=0;it<njt;it++)   u2[ix][it]=v1[it*jt];
        u3[ix][0]=v1[0];  // imag iz =0

		fftwf_execute(plan);
		for(iw=0;iw<nw;iw++) 
		{
			u1[iw][ix*2]=v1[iw*2]; 
			u1[iw][ix*2+1]=v1[iw*2+1]; 
		}
	}
	sf_floatwrite(u2[0],njt*nx,wave); // wave slice iz=0

	fftwf_destroy_plan(plan);

	plan = fftwf_plan_dft_c2r_1d(nt, (fftwf_complex*)v1,
		v1,  FFTW_ESTIMATE);

	sf_psefd_init(nx,nw,dx,dz,
		1.0/(nt*dt),vel);


	for(iz=1;iz<nz;iz++)
	{
		sf_psefd_step3(iz,u1);

		for(ix=0;ix<nx;ix++)
		{
			for(iw=0;iw<nw;iw++)
			{
				v1[iw*2]=u1[iw][ix*2];
				v1[iw*2+1]=u1[iw][ix*2+1];
			}
			fftwf_execute(plan);
			for(it=0;it<njt;it++)	u2[ix][it]=v1[it*jt]/nt;
			u3[ix][iz]=v1[0]/nt;
		}
		sf_floatwrite(u2[0],njt*nx,wave); 

		sf_warning("%d;",iz);
	}

	sf_floatwrite(u3[0],nz*nx,imag);
 
	fftwf_destroy_plan(plan);
	sf_psefd_exit();

	return 0;
}

