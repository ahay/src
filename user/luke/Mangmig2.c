/* Angle-gather constant-velocity time migration. */
/*
  Copyright (C) 2017 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// zero out partial image
void zeroimage(float **image, int ntau, int nxi){
	int tau, xi;
	for( xi=0; xi < nxi; xi++){
		for (tau = 0; tau<ntau; tau++){
			
			image[xi][tau] = 0.;
		}				
	}	
	return;	
}
// dero out data
void zerodata(float ***data, int nt, int nx, int nh){
	int t, x, h;
	for( h=0; h < nh; h++){
		for (x = 0; x<nx; x++){
			for ( t = 0; t<nt; t++){			
				data[h][x][t] = 0.;
			}
		}				
	}	
	return ;	
}


void moduloso(float x, float dx, float x0, int *ix, float *rx){
	float xpre;
	
	xpre = (x-x0)/dx;
	*ix = floor(xpre);
	*rx = xpre - floor(xpre);
	return;
}

void interpolate3(float ***data, int it, float t, int ix, float x, int ih, float h, float *inter){
	*inter = 
		          (1.-h) * (1.-x) * (1.-t) * data[ih  ][ix  ][it  ]
				+ (   h) * (1.-x) * (1.-t) * data[ih+1][ix  ][it  ]
				+ (   h) * (   x) * (1.-t) * data[ih+1][ix+1][it  ]	
				+ (   h) * (   x) * (   t) * data[ih+1][ix+1][it+1]	
				+ (1.-h) * (   x) * (   t) * data[ih  ][ix+1][it+1]					
				+ (1.-h) * (1.-x) * (   t) * data[ih  ][ix  ][it+1]	
				+ (   h) * (1.-x) * (   t) * data[ih+1][ix  ][it+1]
				+ (1.-h) * (   x) * (1.-t) * data[ih  ][ix+1][it  ] ;
	return;
}
void interpolate2p5(float ***data, int it, float t, int ix, float x, int ih, float *inter){
	*inter = 
		          (1.-x) * (1.-t) * data[ih  ][ix  ][it  ]
				+ (   x) * (   t) * data[ih  ][ix+1][it+1]					
				+ (1.-x) * (   t) * data[ih  ][ix  ][it+1]	
				+ (   x) * (1.-t) * data[ih  ][ix+1][it  ] ;
	return;
}

void interpolate2(float **data, int it, float t, int ix, float x,float *inter){
	*inter = 
		          (1.-x) * (1.-t) * data[ix  ][it  ]
				+ (   x) * (1.-t) * data[ix+1][it  ]	
				+ (   x) * (   t) * data[ix+1][it+1]						
				+ (1.-x) * (   t) * data[ix  ][it+1];
	return;
}

void datawriting( float ***data, int it, float t, int ix, float x, int ih, float h, float imgpt){
	
     data[ih  ][ix  ][it  ] += (1.-h) * (1.-x) * (1.-t) * imgpt;
	 data[ih+1][ix  ][it  ]	+= (   h) * (1.-x) * (1.-t) * imgpt;
	 data[ih+1][ix+1][it  ] += (   h) * (   x) * (1.-t) * imgpt;
	 data[ih+1][ix+1][it+1] += (   h) * (   x) * (   t) * imgpt;
	 data[ih  ][ix+1][it+1]	+= (1.-h) * (   x) * (   t) * imgpt;			
	 data[ih  ][ix  ][it+1]	+= (1.-h) * (1.-x) * (   t) * imgpt;
	 data[ih+1][ix  ][it+1] += (   h) * (1.-x) * (   t) * imgpt;
	 data[ih  ][ix+1][it  ] += (1.-h) * (   x) * (1.-t) * imgpt;
	 
	return;
}
void datawriting2p5( float ***data, int it, float t, int ix, float x, int ih, float imgpt){
	
     data[ih  ][ix  ][it  ]  += (1.-x) * (1.-t) * imgpt;
	 data[ih  ][ix+1][it+1]	+= (   x) * (   t) * imgpt;			
	 data[ih  ][ix  ][it+1]	+= (1.-x) * (   t) * imgpt;
	 data[ih  ][ix+1][it  ]  += (   x) * (1.-t) * imgpt;
	 
	return;
}

void angle_engine(bool adj,
	float ***data,int nt,int nx,int nh,float dt,float dx,float dh,float t0,float x0,float h0,
	float **image,int ntau,int nxi,float dtau,float dxi,float tau0,float xi0,
	float **velFile,int nvtau,int nvxi,float dvtau,float dvxi,float vtau0,float vxi0,
	float tpre,float hpre,float xpre, float wtpre){
		
		int itau, ixi;
		float tau, xi;
		float z, t1, h1, x1;
		int it, ix, ih, ivtau, ivxi;
		float rt, rx, rh, rvtau, rvxi;
		float wt, inter, vel;
		
		vel = 1.;
		
		for (ixi=0; ixi<nxi; ixi++){
			xi = xi0 + dxi*ixi;	
			
			ivxi=0; rvxi=0;
			moduloso(xi,dvxi,vxi0,&ivxi,&rvxi);
			if (ivxi<0){continue;}
			if (ivxi>nvxi-2){continue;}		
			
			for (itau=0; itau<ntau; itau++){
				
				tau = tau0 + itau*dtau;
				
				ivtau=0; rvtau=0;
				moduloso(tau,dvtau,vtau0,&ivtau,&rvtau);
				if (ivtau<0){continue;}
				if (ivtau>nvtau-2){continue;}
				
				interpolate2(velFile,ivtau, rvtau, ivxi, rvxi, &vel);
				
				z = tau*vel/2;				
				t1 = tau*tpre;
				h1 = z*hpre;
				x1 = z*xpre+xi;
				wt = wtpre*z*z;
				
				//sf_warning("%g",wt);
				
				it=0; rt=0;
				moduloso(t1,dt,t0,&it,&rt);
				if (it<0){continue;}
				if (it>nt-2){continue;}
								
				ix=0; rx=0;
				moduloso(x1,dx,x0,&ix,&rx);
				if (ix<0){continue;}
				if (ix>nx-2){continue;}	
				
				ih=0; rh=0;
				moduloso(h1,dh,h0,&ih,&rh);
				if (ih<0){continue;}
				if (ih == nh-1){
					
						if (adj){

							inter = image[ixi][itau]*wt;
							datawriting2p5(data,it, rt, ix, rx, ih,inter);
							
						}else{											
							inter = 0.;
							interpolate2p5(data,it, rt, ix, rx, ih,&inter);

							inter *= wt;				
							image[ixi][itau] += inter;
						}
									
				}else{
					if (ih>nh-1){continue;}	
				
					if(adj){
						inter = image[ixi][itau]*wt;					
						datawriting(data,it, rt, ix, rx, ih, rh, inter);
					}else{
						inter = 0.;
						interpolate3(data,it, rt, ix, rx, ih, rh, &inter);
						inter *= wt;				
						image[ixi][itau] += inter;
					}
					
				}
			}
		}
			
		return;
}

void stack(float **image,int ntau,int nxi,float **stk, float **sqstk){

	float img;
	int ixi, itau;
	for (ixi=0; ixi<nxi; ixi++){
		for (itau=0; itau<ntau; itau++){
			img = image[ixi][itau];
			stk[ixi][itau] += img;
			sqstk[ixi][itau] += img*img;
			
		}
		
	} 
	
	return;
}

void semblance(float **stk,float **sqstk,int ntau,int nxi){
	
	int ixi, itau;
	float sqst;
	float eps = 0.00001;
	
	for (ixi=0; ixi<nxi; ixi++){
		for (itau=0; itau<ntau; itau++){
			sqst = sqstk[ixi][itau];

			sqstk[ixi][itau] = stk[ixi][itau]*sqst/(eps+sqst*sqst);
		}
	}
	
	
	return;
}

int main(int argc, char *argv[])
{
	int ia, na, ig, ng, nx, nxi, ntau, nt, nh;
	float a, amax, amin, da, g, dg, gmax;
	float dtau, dt, dh, dxi, dx, tau0, t0, h0, xi0, x0;
	float cosa, sina, sing, cosg;
	float tpre, hpre, xpre, deno, wtpre;
	int nvtau, nvxi;
	float vtau0, vxi0, dvtau, dvxi;

	float **image=NULL;
	float ***data=NULL;
	float **sqstk=NULL;
	float **stk=NULL;
	float **velFile=NULL;
	
	bool adj, sembool;
	
	sf_file in=NULL, out=NULL;
	sf_file	semb=NULL;
	sf_file	vin=NULL;
		
	sf_init (argc,argv);
	in = sf_input("in");
	out = sf_output("out");

		
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if y modeling, if n, migration */
	if ( NULL != sf_getstring("semb") ) {
	/* output file containing Semblance */ 
			semb  = sf_output ("semb"); sembool = true; 

	} else { sembool=false; }
	
	if( adj){ sembool=false;} // no semblance if modelig
	
	if (!adj){ // migration
		if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
		if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input");

		if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"d3",&dh)) sf_error("No d3= in input");

		if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&x0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o3",&h0)) sf_error("No o3= in input");
		
		if (!sf_getint("na",&na)) sf_error("Need na=");
		/* number of dip angles */
		if (!sf_getfloat("amax",&amax)) sf_error("Need amax=");
		/* maximum dip angle */
		
		amin = -1.*amax;
		da = (amax-amin)/(na-1);

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
		if (!sf_getint("ntau",&ntau)) ntau=nt;
		/* output vertical samples */
		if (!sf_getfloat("tau0",&tau0)) tau0=t0;
		/* output vertical orgin */
		if (!sf_getfloat("dtau",&dtau)) dtau=dt;
		/* output vertical sampling */
		
		// output parameters
		sf_putint(out,"n1",ntau);
		sf_putfloat(out,"d1",dtau);
		sf_putfloat(out,"o1",tau0);
		sf_putint(out,"n2",nxi);
		sf_putfloat(out,"d2",dxi);
		sf_putfloat(out,"o2",xi0);
		sf_putint(out,"n3",na);
		sf_putfloat(out,"d3",da);
		sf_putfloat(out,"o3",amin);
				
	} else { // modeling
		
		if (!sf_histint(in,"n1",&ntau)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&nxi)) sf_error("No n2= in input");
		if (!sf_histint(in,"n3",&na)) sf_error("No n3= in input");

		if (!sf_histfloat(in,"d1",&dtau)) sf_error("No d1= in input");
		if (!sf_histfloat(in,"d2",&dxi)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"d3",&da)) sf_error("No d3= in input");

		if (!sf_histfloat(in,"o1",&tau0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o2",&xi0)) sf_error("No o1= in input");
		if (!sf_histfloat(in,"o3",&amin)) sf_error("No o3= in input");
		
		amax = -1.*amin;
		
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
		if (!sf_getint("nt",&nt)) nt=ntau;
		/* number time samples */
		if (!sf_getfloat("t0",&t0)) t0=tau0;
		/* time orgin */
		if (!sf_getfloat("dt",&dt)) dt=dtau;
		/* time increment */
		
		// output parameters
		sf_putint(out,"n1",nt);
		sf_putfloat(out,"d1",dt);
		sf_putfloat(out,"o1",t0);
		sf_putint(out,"n2",nx);
		sf_putfloat(out,"d2",dx);
		sf_putfloat(out,"o2",x0);
		sf_putint(out,"n3",nh);
		sf_putfloat(out,"d3",dh);
		sf_putfloat(out,"o3",h0);
				
	}
	if ( NULL != sf_getstring("vin") ) {
	/* input velocity file */ 
		vin  = sf_input ("vin"); 

		if (!sf_histint(vin,"n1",&nvtau)) sf_error("No n1= in velocity");	
		if (!sf_histint(vin,"n2",&nvxi)) sf_error("No n2= in velocity");


		if (!sf_histfloat(vin,"d1",&dvtau)) sf_error("No d1= in velocity");
		if (!sf_histfloat(vin,"d2",&dvxi)) sf_error("No d2= in velocity");


		if (!sf_histfloat(vin,"o1",&vtau0)) sf_error("No o1= in velocity");
		if (!sf_histfloat(vin,"o2",&vxi0)) sf_error("No o1= in velocity");

		velFile = sf_floatalloc2(nvtau,nvxi);
		sf_floatread (velFile[0],nvtau*nvxi,vin);

	}
	else{
		sf_error("need vel=");
	}
			
	data = sf_floatalloc3(nt,nx,nh);
	image = sf_floatalloc2(ntau,nxi);
	
	
	if (sembool){
		stk = sf_floatalloc2(ntau,nxi);
		sqstk = sf_floatalloc2(ntau,nxi);
		zeroimage(stk,ntau,nxi);
		zeroimage(sqstk,ntau,nxi);
		
		sf_putint(semb,"n1",ntau);
		sf_putfloat(semb,"d1",dtau);
		sf_putfloat(semb,"o1",tau0);
		sf_putint(semb,"n2",nxi);
		sf_putfloat(semb,"d2",dxi);
		sf_putfloat(semb,"o2",xi0);	
		sf_putint(semb,"n3",1);
		sf_putfloat(semb,"d3",da);
		sf_putfloat(semb,"o3",0.);
		
	}
	if (adj){zerodata(data,nt,nx,nh);}
	else{sf_floatread(data[0][0],nt*nx*nh,in);}

	for (ia=0; ia<na; ia++){
		sf_warning("angle %d of %d;",ia+1,na);
				
		if (!adj){zeroimage(image,ntau,nxi);}
		else{sf_floatread (image[0],nxi*ntau,in);} // read partial image		
		
		a = amin + da*ia;
		
		a = a * SF_PI/180; // radians

		cosa = cosf(a);
		sina = sinf(a);
		
		gmax = SF_PI/2 - fabs(a);
		
		if (gmax<=0){continue;}
		
		dg = gmax / (ng-1);
#ifdef _OPENMP 
#pragma omp parallel for
#endif	
		for (ig=0; ig<ng; ig++){
			
			g = ig*dg;
			
			cosg = cosf(g);
			sing = sinf(g);
			
			deno = (cosa*cosa-sing*sing);
			
			if (deno <= 0){continue;}
			
			tpre = cosa*cosg/deno;
			hpre = sing*cosg/deno;
			xpre = sina*cosa/deno;
			
			wtpre = 1/(deno*deno)*dg; // jacobian weighting
			
			angle_engine(adj,data,nt,nx,nh,dt,dx,dh,t0,x0,h0,image,ntau,nxi,dtau,dxi,tau0,xi0,velFile,nvtau,nvxi,dvtau,dvxi,vtau0,vxi0,tpre,hpre,xpre,wtpre);			

		}

		if(!adj){
			sf_floatwrite(image[0],ntau*nxi,out); //write gather				
		}
		if (sembool){
			stack(image,ntau,nxi,stk,sqstk);
			
		}
	}
	if (adj){
		sf_floatwrite(data[0][0],nh*nx*nt,out); //write data
	}
	if (sembool){
		semblance(stk,sqstk,ntau,nxi);
		sf_floatwrite(sqstk[0],nxi*ntau,semb); 
	}
	exit(0);
}


