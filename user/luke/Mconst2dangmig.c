#include <rsf.h>

// dero out partial image
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

void angle_engine(bool adj,float vel,
	float ***data,int nt,int nx,int nh,float dt,float dx,float dh,float t0,float x0,float h0,
	float **image,int ntau,int nxi,float dtau,float dxi,float tau0,float xi0,
	float tpre,float hpre,float xpre){
		
		int itau, ixi;
		float tau, xi;
		float z, t1, h1, x1;
		int it, ix, ih;
		float rt, rx, rh;
		float wt, inter;
		
		wt = 1.;
		
		for (ixi=0; ixi<nxi; ixi++){
			xi = xi0 + dxi*ixi;			
			for (itau=0; itau<ntau; itau++){
				

				tau = tau0 + itau*dtau;
				
				z = tau*vel/2;				
				t1 = tau*tpre;
				h1 = z*hpre;
				x1 = z*xpre+xi;
				
				it=0; rt=0;
				moduloso(t1,dt,t0,&it,&rt);
				if (it<0){continue;}
				if (it>nt+2){continue;}
				
				ih=0; rh=0;
				moduloso(h1,dh,h0,&ih,&rh);
				if (ih<0){continue;}
				if (ih>nh-2){continue;}
				
				ix=0; rx=0;
				moduloso(x1,dx,x0,&ix,&rx);
				if (ix<0){continue;}
				if (ix>nx-2){continue;}				
				
				if(adj){
					inter = image[ixi][itau]*wt;
					
					datawriting(data,it, rt, ix, rx, ih, rh, inter);
				}
				else{
					inter = 0.;
					interpolate3(data,it, rt, ix, rx, ih, rh, &inter);
					inter *= wt;
				
					image[ixi][itau] += inter;
					
				}
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
	float tpre, hpre, xpre, deno;
	float vel;

	float **image=NULL;
	float ***data=NULL;
	
	bool adj;
	
	sf_file in=NULL, out=NULL;
	
	sf_init (argc,argv);
	in = sf_input("in");
	out = sf_output("out");
	
	
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if y modeling, if n, migration */
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
	
	if (!sf_getfloat("vel",&vel)) {sf_error("Need vel=");};
	/*optional constant velocity*/
			
	data = sf_floatalloc3(nt,nx,nh);
	image = sf_floatalloc2(ntau,nxi);
	
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
		
		for (ig=0; ig<ng; ig++){
			
			g = ig*dg;
			
			cosg = cosf(g);
			sing = sinf(g);
			
			deno = (cosa*cosa-sing*sing);
			
			if (deno <= 0){continue;}
			
			tpre = cosa*cosg/deno;
			hpre = sing*cosg/deno;
			xpre = sina*cosa/deno;
			
			angle_engine(adj,vel,data,nt,nx,nh,dt,dx,dh,t0,x0,h0,image,ntau,nxi,dtau,dxi,tau0,xi0,tpre,hpre,xpre);			

		}

		if(!adj){
			sf_floatwrite(image[0],ntau*nxi,out); //write gather				
		}
	}
	if (adj){
		sf_floatwrite(data[0][0],nh*nx*nt,out); //write data
	}
	exit(0);
}


