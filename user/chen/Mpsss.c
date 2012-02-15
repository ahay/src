/* EFD phase shift wave extrapolation. */

#include <rsf.h>

#include "psss.h"
#include "rfft2.h"


int main(int argc, char* argv[])
{
	int ix, iz, jt, njt;
	float dt,dx,dz;
	int nx,nz,nt,nz0,nw,nt1,nx1;
	float ox, **ptx, **pim, **vel, *vv;  
	sf_complex **pfk;
	void *h;

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
	/* depth step for wave data */
	if (!sf_getint("nz",&nz)) nz=nz0; 
	/* depth number */

	sf_putint(imag,"n1",nz);
	sf_putint(imag,"n2",nx);
	sf_putfloat(imag,"o2",ox);

	njt=1+(nz-1)/jt;

	sf_putint(wave,"n1",nt);
	sf_putfloat(wave,"d1",dt);
	sf_putfloat(wave,"o1",0);
	sf_putint(wave,"n2",nx);
	sf_putfloat(wave,"d2",dx);
	sf_putfloat(wave,"o2",ox);
	sf_putint(wave,"n3",njt);
	sf_putfloat(wave,"d3",dz*jt);
	sf_putfloat(wave,"o3",0);

	nt1 = nt;
	nx1 = nx;
	h = sf_rfft2_init(&nt1, &nx1, &nw, 1);

	/* read data and velocity */
	vel = sf_floatalloc2(nz,nx);
	vv  = sf_floatalloc(nz);
	sf_floatread(vel[0],nz*nx,modl);
	for(iz=0; iz<nz; iz++)	
	{
		vv[iz] = 0.0;
		for(ix=0; ix<nx; ix++)	vv[iz] += vel[ix][iz];
		vv[iz] /= nx;
	}

	ptx  = sf_floatalloc2(nt,nx);	// U_z(t,x)
	pim  = sf_floatalloc2(nz,nx);		// u(x,z,0)
	pfk  = sf_complexalloc2(nw,nx1);

	sf_floatread(ptx[0],nt*nx,data);

	for(ix=0;ix<nx;ix++)	
		pim[ix][0]=ptx[ix][0];  // imag iz =0

	sf_rfft2(h, ptx, pfk);

	sf_floatwrite(ptx[0],nt*nx,wave); // wave slice iz=0

	sf_psss_init(nw,nx1,nz,
		1.0/(dt*nt1),dx,dz,vv);


	for(iz=1;iz<nz;iz++)
	{
		sf_psss_step(iz, pfk);

		sf_rifft2(h, pfk, ptx);
		for(ix=0;ix<nx;ix++)	
			pim[ix][iz]=ptx[ix][0];  // imag iz =0

		if(iz % jt == 0)
			sf_floatwrite(ptx[0],nt*nx,wave); 
		sf_warning("%d;",iz);
	}

	sf_floatwrite(pim[0],nz*nx,imag);
 
	sf_psss_release();
	sf_rfft2_release(h);

	return 0;
}

