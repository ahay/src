/* EFD phase shift wave extrapolation. */

#include <rsf.h>

#include "psefd.h"
#include "rfft1.h"


int main(int argc, char* argv[])
{
	int it,ix,iz,jt,njt;
	float dt,dx,dz;
	int nw,nx,nz,nt,nz0,nt2;

	float *v1; 
	sf_complex **u1; 

	float ox;
	float **u2,**u3,**vel;  
	void * h;

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

	nt2 = nt;
	h = sf_rfft1_init(&nt2, &nw);
    v1 = sf_floatalloc(nt);	

    /* read data and velocity */
    vel = sf_floatalloc2(nz,nx);
    sf_floatread(vel[0],nz*nx,modl);

    u1 = sf_complexalloc2(nw,nx);	// U_z(w,x)
    u2 = sf_floatalloc2(njt,nx);	// u_z(x,t)
    u3 = sf_floatalloc2(nz,nx);		// u(x,z,0)


	for(ix=0;ix<nx;ix++)	
	{
    	sf_floatread(v1,nt,data);
		for(it=0;it<njt;it++)   u2[ix][it]=v1[it*jt];
        u3[ix][0]=v1[0];  // imag iz =0

		sf_rfft1(h, v1, u1[ix]);
//		sf_rifft1(h,  u1[ix],v1);
	}
	sf_floatwrite(u2[0],njt*nx,wave); // wave slice iz=0

	sf_psefd_init(nx,nw,dx,dz,
		1.0/(nt2*dt),vel);


	for(iz=1;iz<nz;iz++)
	{
		sf_psefd_step3(iz,u1);

		for(ix=0;ix<nx;ix++)
		{
			sf_rifft1(h, u1[ix], v1);
			for(it=0;it<njt;it++)	u2[ix][it]=v1[it*jt];
			u3[ix][iz]=v1[0];
		}
		sf_floatwrite(u2[0],njt*nx,wave); 

		sf_warning("%d;",iz);
	}

	sf_floatwrite(u3[0],nz*nx,imag);
 
	sf_psefd_release();

	return 0;
}

