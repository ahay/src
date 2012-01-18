/* EFD phase shift wave extrapolation. */

#include <rsf.h>

#include "psss.h"


int main(int argc, char* argv[])
{
	int it,ix,iz,jt,njt;
	float dt,dx,dz;
	int nx,nz,nt,nz0,nw,nt1,nx1;
	float ox, **ptx, **ptx1,**pim,**pwvz,**vel;  
	kiss_fft_cpx **pfk;

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

	nw=fft2_init(0,1,nt,nx, &nt1, &nx1);

    /* read data and velocity */
    vel = sf_floatalloc2(nz,nx);
    sf_floatread(vel[0],nz*nx,modl);

    ptx  = sf_floatalloc2(nt,nx);	// U_z(t,x)
    ptx1 = sf_floatalloc2(nt1,nx1);	// U_z(t,x)
    pwvz = sf_floatalloc2(njt,nx);	// u_z(x,t)
    pim  = sf_floatalloc2(nz,nx);		// u(x,z,0)
	pfk  = sf_complexalloc2(nw,nx1);

    sf_floatread(ptx[0],nt*nx,data);

	for(ix=0;ix<nx;ix++)	
	{
		for(it=0;it<njt;it++)  	pwvz[ix][it]=ptx[ix][it*jt]; 
		for(it=0;it<nt;it++)  	ptx1[ix][it]=ptx[ix][it]; 
		for(it=nt;it<nt1;it++)  	ptx1[ix][it]=0.0; 
        pim[ix][0]=pwvz[ix][0];  // imag iz =0
	}
	for(ix=nx;ix<nx1;ix++)	memset(ptx1[ix],0,nt1*sizeof(float));

	fft2(ptx1[0],pfk[0]);

	sf_floatwrite(pwvz[0],njt*nx,wave); // wave slice iz=0

	sf_psss_init(nw,nx1,nz,
		1.0/(dt*nt1),dx,dz,vel);


	for(iz=1;iz<nz;iz++)
	{
		sf_psss_step(iz,pfk);

		ifft2(pfk[0],ptx1[0]);
		for(ix=0;ix<nx;ix++)	
		{
			for(it=0;it<njt;it++)   pwvz[ix][it]=ptx1[ix][it*jt];
			pim[ix][iz]=pwvz[ix][0];  // imag iz =0
		}
		sf_floatwrite(pwvz[0],njt*nx,wave); 

		sf_warning("%d;",iz);
	}

	sf_floatwrite(pim[0],nz*nx,imag);
 
	sf_psss_exit();

	return 0;
}

