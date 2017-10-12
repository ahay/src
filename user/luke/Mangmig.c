
#include <rsf.h>



int main(int argc, char *argv[])
{
	int ia, na, htrack;
	float amax, amin, da, a, ar, dar, absar;
	int ig, ng;
	float g0, gmax, gmin, dg, g;
	int iz, nz;
	float z0, dz, z;
	int ixi, nxi;
	float xi0, dxi, xi;
	int ix, nx, nvx, ivx;
	float x0, dx, vx0, dvx, ivx1, vx;
	int ih, nh;
	float h0, dh;
	int it, nt, nvt, ivt;
	float t0, dt, vt0, dvt, ivt1, vt;
	float vel, cosa, sina, cosg, sing;
	float deno, tpre, hpre, xpre;
	float t1, x1, h1;
	float t2, x2, h2;
	float t, x, h, wt, inte;
	float dA, jac, tamppre, eps;
	bool sembool, velbool, adj;
	
	float **image=NULL;
	float ***data=NULL;
	float **velFile=NULL;
	float **sembl=NULL;
	float **stack=NULL;
	float **sqstack=NULL;

	sf_file in=NULL, out=NULL;
	sf_file	semb=NULL, vin=NULL;

	
	eps = 0.02;

	sf_init (argc,argv);
	in = sf_input("in");
	out = sf_output("out");
	
	
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if y modeling, if n, migration */
	adj = !adj; // needs to be this way for conjgrad to work
	
	if ( NULL != sf_getstring("semb") ) {
	/* output file containing Semblance */ 
			semb  = sf_output ("semb"); sembool = true; 

	} else { sembool=false; }

	if( adj){ sembool=false;}
	
	if ( adj){
		sf_warning("adj");
		if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&nxi)) sf_error("No n2= in input");
		if (!sf_histint(in,"n3",&na)) sf_error("No n3= in input");

		if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dxi)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"d3",&da)) sf_error("No d3= in input");

		if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&xi0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o3",&amin)) sf_error("No o3= in input");

		image = sf_floatalloc2(nz,nxi);
		
		amax = -amin;
	} else{
		sf_warning("fwd");
		if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
		if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input");

		if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"d3",&dh)) sf_error("No d3= in input");

		if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&x0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o3",&h0)) sf_error("No o3= in input");

		data = sf_floatalloc3(nt,nx,nh);
		
		
	} // end forward input

	if ( NULL != sf_getstring("vin") ) {
	/* input velocity file */ 
		vin  = sf_input ("vin"); 
		velbool = true;

		if (!sf_histint(vin,"n1",&nvt)) sf_error("No n1= in velocity");
		
		if (!sf_histint(vin,"n2",&nvx)) sf_error("No n2= in velocity");


		if (!sf_histfloat(vin,"d1",&dvt)) sf_error("No d1= in velocity");
		if (!sf_histfloat(vin,"d2",&dvx)) sf_error("No d2= in velocity");


		if (!sf_histfloat(vin,"o1",&vt0)) sf_error("No o1= in velocity");
		if (!sf_histfloat(vin,"o2",&vx0)) sf_error("No o1= in velocity");

		velFile = sf_floatalloc2(nvt,nvx);
		sf_floatread (velFile[0],nvt*nvx,vin);

	}
	else{ velbool=false;
			if (!sf_getfloat("vel",&vel)) {sf_error("Need vel=");};
			/*optional constant velocity*/
	}
	/* velocity */

	g0 = 0; //reflection angles centered at zero

	if (adj){ //modeling
		
		if (!sf_getint("nh",&nh)) sf_error("Need nh=");
		/* number of offsets */
		if (!sf_getfloat("h0",&h0)) sf_error("Need h0=");
		/* initial offset */
		if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
		/* offset increment*/
		
		if (!sf_getint("ng",&ng)) ng=na;
		/* number of reflection angles */
		if (!sf_getfloat("gmax",&gmax)) gmax=amax;
		/* maximum reflection angle*/

		if (!sf_getint("nx",&nx)) nx=nxi;
		/* data domain spatial samples */
		if (!sf_getfloat("x0",&x0)) x0=xi0;
		/* data domain spatial orgin */
		if (!sf_getfloat("dx",&dx)) dx=dxi;
		/* data domain spatial increment */
		if (!sf_getint("nt",&nt)) nt=nz;
		/* number time samples */
		if (!sf_getfloat("t0",&t0)) t0=z0;
		/* time orgin */
		if (!sf_getfloat("dt",&dt)) dt=dz;
		/* time increment */
		htrack = 0;
		if (nh==1){
			nh=2 ;
			htrack = 1;
		}

		data = sf_floatalloc3(nt,nx,nh);

		for (ih=0; ih<nh; ih++){
			for (ix=0; ix<nx; ix++){
				for (it=0; it<nt; it++){
					
					data[ih][ix][it] = 0.;
				}
				
			}
		}

		dar = da * SF_PI / 180.;

		dg = (gmax-g0) / (ng-1.);



		sf_putint(out,"n1",nt);
		sf_putfloat(out,"d1",dt);
		sf_putfloat(out,"o1",t0);
		sf_putint(out,"n2",nx);
		sf_putfloat(out,"d2",dx);
		sf_putfloat(out,"o2",x0);
		sf_putint(out,"n3",nh);
		sf_putfloat(out,"d3",dh);
		sf_putfloat(out,"o3",h0);
		
		
	}else{ //migration
		if (!sf_getint("na",&na)) sf_error("Need na=");
		/* number of dip angles */
		if (!sf_getfloat("amax",&amax)) sf_error("Need amax=");
		/* maximum dip angle */

		if (!sf_getint("ng",&ng)) ng=na;
		/* number of reflection angles */
		if (!sf_getfloat("gmax",&gmax)) gmax=amax;
		/* maximum reflection angle*/

		if (!sf_getint("nxi",&nxi)) nxi=nx;
		/* output samples */
		if (!sf_getfloat("xi0",&xi0)) xi0=x0;
		/* output orgin */
		if (!sf_getfloat("dxi",&dxi)) dxi=dx;
		/* output sampling */
		if (!sf_getint("nz",&nz)) nz=nt;
		/* output vertical samples */
		if (!sf_getfloat("z0",&z0)) z0=t0;
		/* output vertical orgin */
		if (!sf_getfloat("dz",&dz)) dz=dt;
		/* output vertical sampling */

		image = sf_floatalloc2(nz,nxi);


		sf_floatread (data[0][0],nt*nx*nh,in);

		amin = -1.*amax;
		da = (amax-amin)/(na-1);
		dar = da * SF_PI / 180.;

		dg = (gmax-g0) / (ng-1.);



		sf_putint(out,"n1",nz);
		sf_putfloat(out,"d1",dz);
		sf_putfloat(out,"o1",z0);
		sf_putint(out,"n2",nxi);
		sf_putfloat(out,"d2",dxi);
		sf_putfloat(out,"o2",xi0);
		sf_putint(out,"n3",na);
		sf_putfloat(out,"d3",da);
		sf_putfloat(out,"o3",amin);
	} // end 2nd round fwd setup

	if (sembool){
		sembl = sf_floatalloc2(nz,nxi);
		stack = sf_floatalloc2(nz,nxi);	
		sqstack = sf_floatalloc2(nz,nxi);
		sf_putint(semb,"n1",nz);
		sf_putfloat(semb,"d1",dz);
		sf_putfloat(semb,"o1",z0);
		sf_putint(semb,"n2",nxi);
		sf_putfloat(semb,"d2",dxi);
		sf_putfloat(semb,"o2",xi0);	
		sf_putint(semb,"n3",1);
		sf_putfloat(semb,"d3",da);
		sf_putfloat(semb,"o3",0.);
		

	}

	if (sembool){
	//initialize semblance
		for (ixi=0; ixi<nxi; ixi++){
			for (iz=0; iz<nz ; iz++){

				sembl[ixi][iz] = 0.;
				stack[ixi][iz] = 0.;
				sqstack[ixi][iz]=0.;
			}
		}
	}

	for (ia=0; ia<na; ia++){ // loop over dip angles
		a = amin + da * ia*1.;
		// to radians
		ar = a * SF_PI / 180.;

		sf_warning("angle %d of %d;",ia+1,na);
		// initialize array

		if (adj){
			sf_floatread (image[0],nxi*nz,in);
		}else{
			for (ixi=0; ixi<nxi; ixi++){
				for (iz=0; iz<nz ; iz++){

					image[ixi][iz] = 0.;
				}
			}
		}
		cosa = cosf(ar);
		sina = sinf(ar);

		absar = fabs(ar);
		

		gmax = fabs(SF_PI/2 - absar);

		
		gmin = 0.;
		//gmin = 0;
		dg = (gmax-gmin)/(ng-1.);


		for (ixi=0; ixi<nxi; ixi++){
			xi = xi0 + ixi*dxi;


			for (iz=0; iz<nz ; iz++){
				z = z0 + iz*dz;
				dA = dar*dg;
				
				
				// velocity
				if (velbool){
					ivt1 = (vt0-z)/dvt ;
					ivx1 = (vx0-xi)/dvx ;

					ivt = floor(ivt1);
					ivx = floor(ivx1);
					vt = ivt1-ivt;
					vx = ivx1-ivx;

					if (ivt < 0){ivt = 0; vt=0.;}
					if (ivt > nvt-2){ivt = nvt-2; vt=1. ;}

					if (ivx < 0){ivx = 0; vx = 0. ;}
					if (ivx > nvx-2){ivx = nvx-2; vx = 1.;}

					vel = (1.-vx)*(1.-vt)*velFile[ivx  ][ivt  ] +
					      (vx   )*(1.-vt)*velFile[ivx+1][ivt  ] +
					      (vx   )*(   vt)*velFile[ivx+1][ivt+1] +
					      (1.-vx)*(   vt)*velFile[ivx  ][ivt+1] ;
				}

				for (ig=0; ig<ng; ig++){ //loop over reflection angles

					g = gmin + ig*1. * dg ; 

					cosg = cosf(g);
					sing = sinf(g);

					deno = cosa*cosa - sing*sing ; 
			
					if (fabs(deno)<eps){continue;}
					tpre  = cosa*cosg/deno;
					hpre  = sing*cosg/deno;
					xpre = sina*cosa/deno;

					tamppre = 2./sqrt(2*SF_PI)*sina*(cosa*cosa+sing*sing)/sqrt(deno*deno*deno*deno*deno);
					
					jac = z/deno;

					 wt = dA * jac;




					t1 = 2*z*tpre;
					h1 = 0.5*z*vel*hpre;
					x1 = 0.5*z*vel*xpre+xi;


					t2 = (t1-t0)/dx;

					h2 = (h1-h0)/dh;
					x2 = (x1-x0)/dx;
					// convert to memory indicies
					it = floor(t2);

					if (it<0   ) {continue;}
					if (it>nt-2) {continue;}
//sf_warning("%g",h1);
					ih = floor(h2);
					if (ih<0   ) {continue;}
					if (ih>nh-2) {continue;}

					ix = floor(x2);
					if (ix<0   ) {continue;}
					if (ix>nx-2) {continue;}

					// modulos
					t = t2 - it;
					h = h2 - ih;
					x = x2 - ix;
					
					if(adj){
						inte = image[ixi][iz]*wt;
						
						data[ih  ][ix  ][it  ] += (1.-h) * (1.-x) * (1.-t) * inte;						
						data[ih+1][ix  ][it  ] += (   h) * (1.-x) * (1.-t) * inte;
						data[ih+1][ix+1][it  ] += (   h) * (   x) * (1.-t) * inte;
						data[ih+1][ix+1][it+1] += (   h) * (   x) * (   t) * inte;
						data[ih  ][ix+1][it+1] += (1.-h) * (   x) * (   t) * inte;				
						data[ih  ][ix  ][it+1] += (1.-h) * (1.-x) * (   t) * inte;
						data[ih+1][ix  ][it+1] += (   h) * (1.-x) * (   t) * inte;
						data[ih  ][ix+1][it  ] += (1.-h) * (   x) * (1.-t) * inte;
						
						
					}else{

						// interpolate data points
						inte =    (1.-h) * (1.-x) * (1.-t) * data[ih  ][ix  ][it  ]
								+ (   h) * (1.-x) * (1.-t) * data[ih+1][ix  ][it  ]
								+ (   h) * (   x) * (1.-t) * data[ih+1][ix+1][it  ]	
								+ (   h) * (   x) * (   t) * data[ih+1][ix+1][it+1]	
								+ (1.-h) * (   x) * (   t) * data[ih  ][ix+1][it+1]					
								+ (1.-h) * (1.-x) * (   t) * data[ih  ][ix  ][it+1]	
								+ (   h) * (1.-x) * (   t) * data[ih+1][ix  ][it+1]
								+ (1.-h) * (   x) * (1.-t) * data[ih  ][ix+1][it  ] ;
						inte *= wt;


						image[ixi][iz] += inte;
						
						if (sembool){
							stack[ixi][iz] += inte;
									
							sqstack[ixi][iz] += fabsf(inte);	
							
						}
					}				
					//image[ixi][iz] += wt;
//sf_warning("%g",image[ixi][iz]);
				}
			}


		} // end reflection angles
		if(adj){
		}else{
			sf_floatwrite(image[0],nz*nxi,out);
		}

	} // end dip angles loop
	if (sembool){ 
		for (ixi=0; ixi<nxi; ixi++){
			for (iz=0; iz<nz ; iz++){
	
				if(sqstack[ixi][iz] > 0){ 
					sembl[ixi][iz] = (stack[ixi][iz]/sqstack[ixi][iz])*(stack[ixi][iz]/sqstack[ixi][iz]);

				
				}
			}
		}	
		sf_floatwrite(sembl[0],nz*nxi,semb);
		
	}
	if (adj){
		
		sf_floatwrite(data[0][0],nh*nx*nt,out);
	}
    sf_fileclose (out);
    sf_fileclose (in);
	
    if(velbool){
		sf_fileclose (vin);
		free (velFile);
	}
	
	if(sembool){
		sf_fileclose(semb); 	
		free (sqstack);
		free (sembl);
		free (stack);
	}
	
    free (image);
    free (data);


	exit(0);
}
