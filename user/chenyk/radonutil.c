/* Time domain high-resolution Hyperbolic Radon transform utilities. 
m(tau,p) = \sum_{ih=1}^{nh} d(tau, \sqrt{tau^2+h[ih]^2/p^2)}
*/
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
  
  Reference:
  M.D. Sacchi's lecture in China University of Petroleum-Beijing, 2012.5.
*/
#include <rsf.h>

#include "radonutil.h"
#include "vecoper.h"

static int nt,nh,nv;
static float dt,*v,*h;

void hradon_init(int Nt, int Nh, int Nv, float Dt, float *V, float *H)
/*< Initialization for linear radon operator >*/
{
	nt=Nt;
	nh=Nh;
	nv=Nv;
	dt=Dt;
	v=V;
	h=H;
}

void hradon1(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< hyperbolic radon operator >*/
{   
	int itau,ih,iv,it;
	float tau,h_v,t;
	
	sf_adjnull(adj,add,nx,ny,x,y);
	
	for(itau=0;itau<nt;itau++)
		for(ih=0;ih<nh;ih++)
			for (iv=0;iv<nv;iv++)
			{
				tau=itau*dt;
				h_v=h[ih]/v[iv];
				t=sqrtf(tau*tau+h_v*h_v);
				it=floorf(t/dt)+1;
				if(it<=nt)
				{
					if (adj)
					x[iv*nt+itau]+=y[ih*nt+it];
					else
					y[ih*nt+it]+=x[iv*nt+itau]; 
				}	
			}			

}


void hradon(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< hyperbolic radon operator >*/
{   
  int itau,ih,iv,it,i;
	float tau,h_v,t;
	
	nv=nx/nt;  nh=ny/nt;
	sf_adjnull(adj,add,nx,ny,x,y);

	
	for(i=0;i<nx;i++)
		x[i]=x[i]+0.0001;
			
	for(i=0;i<ny;i++)
		y[i]=y[i]+0.0001;

	for(itau=0;itau<nt;itau++)
		for(ih=0;ih<nh;ih++)
			for (iv=0;iv<nv;iv++)
			{
				tau=itau*dt;
				h_v=h[ih]/v[iv];
				t=sqrtf(tau*tau+h_v*h_v);
				it=floorf(t/dt)+1;
				if(it<=nt)
				{
					if (adj)
					x[iv*nt+itau]+=y[ih*nt+it];
					else
					y[ih*nt+it]+=x[iv*nt+itau]; 
				}	
			}		
}

void hradon_pcg(sf_operator oper, 
			float *d, 
			float *x, 
			int itmax_internal,
			int itmax_external, 
			int verb,
			float *misfit)
/*< Preconditioned CG for sparsity promotion hyperbolic radon transform >*/
{
	int l,k,kc,nx,ny;
	float *z,*P,*s,*g,*xtemp,*di,*ss,*y,*r;
	float gammam,gamma,beta,den,alpha;
	nx=nt*nv;
	ny=nt*nh;
	
	/* Allocate memory */
	xtemp=sf_floatalloc(nx);
	y=sf_floatalloc(nx);														
	z=sf_floatalloc(nx);			
	P=sf_floatalloc(nx);
	di=sf_floatalloc(ny);
	r=sf_floatalloc(ny);
	g=sf_floatalloc(nx);
	s=sf_floatalloc(nx);
	ss=sf_floatalloc(ny);

	
	consvec(0,nx,z);				// z=zeros(size(m0));				
	consvec(1,nx,P);				// P=ones(size(z));	
			
	kc=1;						// kc=1;
	misfit=sf_floatalloc(itmax_internal*itmax_external); // Misfit=[];
	scale(z,x,nx,1);						// x=z;
	
	for(l=1;l<=itmax_external;l++)
	{
	vecmul(P,z,xtemp,nx);
	oper(false,false,nx,ny,xtemp,di);	// feval(operator,P.*z,Param,1);
	scalesum(d, di, r, ny, 1, -1);	// r=d-di;
	oper(true,false,nx,ny,g,r);		// g=feval(operator,r,Param,-1);
	vecmul(g,P,g,nx);				// g=g.*P;
	scale(g,s,nx,1);				// s=g;
	gammam=cblas_sdot(nx,g,1,g,1);	// gammam=cgdot(g);
	k=1;							// k=1;
	while (k<=itmax_internal)
	{
		vecmul(P,s,xtemp,nx);
		oper(false,false,nx,ny,xtemp,ss);	// ss=feval(operator,P.*s,Param,1);
		den = cblas_sdot(ny,ss,1,ss,1);	// den = cgdot(ss);	
		alpha = gammam/(den+SF_EPS);		// alpha=gammam/(den+1.e-8);
		scalesum(z,s,z,nx,1,alpha);			// z=z+alpha*s;
		scalesum(r,ss,r,ny,1,-alpha);			// r=r-alpha*ss;
		misfit[kc-1]=cblas_sdot(ny,r,1,r,1);	// misfit(kc)=cgdot(r);
		oper(true,false,nx,ny,g,r);			// g=feval(operator,r,Param,-1);
		vecmul(g,P,g,nx);					// g=g.*P;
		gamma = cblas_sdot(nx,g,1,g,1);		// gamma=cgdot(g);
		beta=gamma/(gammam+SF_EPS);			// beta=gamma/(gammam+1.e-7);	
		gammam=gamma;						// gammam=gamma;
		
		scalesum(g,s,s,nx,1,beta);			// s=g+beta*s;
		if(verb == 1)
			sf_warning("Iteration = %d Misfit=%0.5g",k,misfit[kc-1]);
		k=k+1;							// k=k+1;
		kc=kc+1;}							// kc=kc+1;
	vecmul(P,z,x,nx);						// x=P.*z;
	vecabs(x,xtemp,nx);						// y=x/max(abs(x(:)));
	scale(x,y,nx,1.0/vecmax(xtemp,nx));		// y=x/max(abs(x(:)));
	vecabs(y,P,nx);						// P=abs(y)+0.001;
	scalesumreal(P,0.001,P,nx,1);				// P=abs(y)+0.001;
	}	
}
























