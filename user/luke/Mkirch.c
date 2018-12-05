
#include <rsf.h>

int main(int argc, char *argv[])
{
	int nz, nt, nx, ny;
	int iz, it, ix, iy;
	float z0, t0, x0, yo;
	float z, t, x, y, h, t1, rt;
	float dz, dt, dx, dy, vel;
	
	float **image=NULL;
	float **data=NULL;

	bool adj;
	
	sf_init (argc,argv);
	sf_file in=NULL, out=NULL;
	sf_file velfile=NULL;
	in = sf_input("in");
	out = sf_output("out");
	
	if (NULL == sf_getstring("vel") ) {sf_error("Need vel=");}
	velfile  = sf_input ("vel"); 
	
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if n, modeling, if y, migration */
	adj = !adj; // need to do it this way for conjgrad to work properly
		
	if (adj){
		
		if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");

		if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");

		if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&yo)) sf_error("No o1= in input");
		if (!sf_getint("nx",&nx)) nx=ny;
		/* image lateral spatial samples */
		if (!sf_getfloat("dx",&dx)) dx=dy;
		/* image lateral sampling */
		if (!sf_getfloat("x0",&x0)) x0=yo;
		/* image lateral start */
		if (!sf_getint("nt",&nt)) nt=nz;
		/* number of depth samplesw */
		if (!sf_getfloat("dt",&dt)) dt=dz;
		/* depth sampling */
		if (!sf_getfloat("t0",&t0)) t0=z0;
		/* depth start */
		

	
		data = sf_floatalloc2(nt,nx);
		image = sf_floatalloc2(nz,ny);
		sf_floatread (image[0],nz*ny,in);
		
		sf_putint(out,"n1",nt);
		sf_putfloat(out,"d1",dt);
		sf_putfloat(out,"o1",t0);
		
		sf_putint(out,"n2",nx);
		sf_putfloat(out,"d2",dx);
		sf_putfloat(out,"o2",x0);
	}
	
	else {
		if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

		if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

		if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&x0)) sf_error("No o1= in input");
		if (!sf_getint("ny",&ny)) ny=nx;
		/* image lateral spatial samples */
		if (!sf_getfloat("dy",&dy)) dy=dx;
		/* image lateral sampling */
		if (!sf_getfloat("y0",&yo)) yo=x0;
		/* image lateral start */
		if (!sf_getint("nz",&nz)) nz=nt;
		/* number of depth samplesw */
		if (!sf_getfloat("dz",&dz)) dz=dt;
		/* depth sampling */
		if (!sf_getfloat("z0",&z0)) z0=t0;
		/* depth start */

		data = sf_floatalloc2(nt,nx);
		image = sf_floatalloc2(nz,ny);
		sf_floatread (data[0],nt*nx,in);
		
		sf_putint(out,"n1",nz);
		sf_putfloat(out,"d1",dz);
		sf_putfloat(out,"o1",z0);
		
		sf_putint(out,"n2",ny);
		sf_putfloat(out,"d2",dy);
		sf_putfloat(out,"o2",yo);
	
	}
	
	float **velarray = sf_floatalloc2(nz,ny);
	sf_floatread (velarray[0],nz*ny,velfile);
	
	if (adj){
		for (it=0; it<nt; it++){
			for (ix=0;  ix<nx; ix++){
				data[ix][it] = 0.;
			}
		}
		
	}
	else{
		for (iz=0; iz<nz; iz++){
			for (iy=0;  iy<ny; iy++){
				image[iy][iz] = 0.;
			}
		}
	}
	for (ix=0; ix<nx; ix++){
		x = x0 + ix*dx;
		for (iy=0; iy<ny; iy++){
			y = yo + iy*dy;

			for (iz=0; iz<nz; iz++){
				vel = velarray[iy][iz];
				
				h = (x-y)/vel;
							
				z = z0 + iz*dz;
				t = sqrt(z*z+h*h);
				
				t1 = (t-t0)/dt;
				it = floor(t1);
				rt = t1-it;
				
				if (it>nt-2){continue;}
				if (adj){
					data[ix][it] += (1-rt)*image[iy][iz];
					data[ix][it+1] += rt*image[iy][iz];
				}else{
					image[iy][iz] += (1-rt)*data[ix][it] + rt*data[ix][it+1];
				}

			}
			
		}
		
	}
	if (adj){
		sf_floatwrite(data[0],nt*nx,out);
	}else{
		sf_floatwrite(image[0],nz*ny,out);
	}
	exit(0);
}