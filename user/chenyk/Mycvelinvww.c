/* Inverse velocity spectrum with interpolation by modeling from inversion result (C version) */
/*
  Copyright (C) 2020 University of Texas at Austin
   
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

void velxf( bool adj, 	/*adj flag*/
			bool add, 	/*add flag*/
			int nt, 	/*nt*/
			float dt, 	/*dt*/
			float ot, 	/*ot*/
			float *x2, 	/*squared x axis*/
			int nx, 	/*number of space samples*/
			float *data, /*input data, e.g., CMP gather*/
			float *mask, /*mask*/
			float *s,	/*velocity axis*/
			int ns, 	/*number of velocities*/
			float *modl,/*velocity spectrum domain*/
			float *z2 	/*squared t axis*/)
{
	int is, ix, iz, it;
	float x2s, t, ft, gt;
	int endt;

	sf_adjnull(adj,add,nt*ns,nt*nx,modl,data);
	
	if(!adj)
	{
	for(is=0;is<ns;is++)
	{
		for(ix=0;ix<nx;ix++)
		{
			x2s=x2[ix]*s[is];
			if(x2s>z2[nt-1]) break;		//?
			if (mask[ix]==0) continue; 	//?
			endt=sqrtf(z2[nt-1]-x2s)/dt+1.0;
			for(iz=0;iz<endt;iz++)
			{
				t=sqrtf(z2[iz]+x2s);
				it=1.0+(t-ot)/dt;
				ft=(it*dt-t)/dt;
				gt=1.0-ft;
				data[it+ix*nt]=data[it+ix*nt]+ft*modl[iz+is*nt];
				data[it+1+ix*nt]=data[it+1+ix*nt]+gt*modl[iz+is*nt];
			}
		}
	}
	}else{
	for(is=0;is<ns;is++)
	{	for(ix=0;ix<nx;ix++)
		{
			x2s=x2[ix]*s[is];
			if(x2s>z2[nt-1])break;		//?
			if(mask[ix]==0)continue;	//?
			endt=sqrtf(z2[nt-1]-x2s)/dt+1.0;
			for(iz=0;iz<endt;iz++)
			{
				t=sqrtf(z2[iz]+x2s);
				it=1.0+(t-ot)/dt;
				ft=(it*dt-t)/dt;
				gt=1.0-ft;
				modl[iz+nt*is]=modl[iz+nt*is]+ft*data[it+ix*nt]+gt*data[it+1+ix*nt];
			}
		}
	}
	}
}


float mmin(int n, float *a)
{
int i;
float min=a[0];
for(i=0;i<n;i++)
	if(min>a[i]) min=a[i];
return min;
}

float mmax(int n, float *a)
{
int i;
float max=a[0];
for(i=0;i<n;i++)
	if(max<a[i]) max=a[i];
return max;
}


void copy(int n, float *x, float *y)
{
int i;
for(i=0;i<n;i++)
y[i]=x[i];
}

void zero(int n, float *x)
{
int i;
for(i=0;i<n;i++)
x[i]=0.0;
}

float my_ddot(int n, float *x, float *y)
{
int i;
float dot=0;
for(i=0;i<n;i++)
dot=dot+x[i]*y[i];
return dot;
}

void powwt(int iter, float *rr, /*input vector*/
		  int n, float *wt, 	/*weighting vector*/
		  float *wrr, 	/*weighting vector*/
		  float pow,	/*weighting power*/
		  float eps,
		  int huber,
		  float srate)
{

  int i; 
  float small, maxval;
  float *arr, *awt; 
  arr=sf_floatalloc(n);
  awt=sf_floatalloc(n);

  if(pow==0.0 || iter==0)
  {		
  copy(n,rr,wrr);
  return;	
  }
  
  for(i=0;i<n;i++)
  	awt[i]=fabsf(wt[i]);
  
  if(pow<0.0)
  {
  	maxval=mmax(n,awt);
  	small=srate*maxval;
  	for(i=0;i<n;i++)
  		if(awt[i]<small) awt[i]=small;
  }
  
  if(huber==1)
  {
  	for(i=0;i<n;i++)
  		arr[i]=fabsf(rr[i]);
  	maxval=mmax(n,arr);
  	small=eps*maxval;
  	for(i=0;i<n;i++)
  		if(arr[i]>small)
  			wrr[i]=rr[i]*powf(awt[i],pow);
  		else
  			wrr[i]=rr[i];
  }else{
  	for(i=0;i<n;i++)
  		wrr[i]=rr[i]*powf(awt[i],pow);
  }
}

void cgstep( int iter, int n, float *x, float *g, float *s, int m, float *rr, float *gg, float *ss)
/*! solution, residual
! gradient, conjugate gradient
! step, conjugate step*/
{
  int i;
  float dot, den,num, sds, gdg, gds, determ, gdr, sdr, alfa, beta;

  if(iter==0)
  {
  	for(i=0;i<n;i++)
  		s[i]=0.0;
   	for(i=0;i<m;i++)
  		ss[i]=0.0; 
  	dot=my_ddot(m,gg,gg);
  	if(dot==0)
  	{	sf_warning("cgstep: grad vanishes identically");
  		return;
  	}		
  	den=my_ddot(m,gg,rr);
  	num=my_ddot(m,gg,gg);
  	alfa=den/num;
  	beta=0.0;
  }else{
  	/*search plane by solving 2-by-2*/
  	gdg=my_ddot(m,gg,gg);	/*G . (R - G*alfa - S*beta) = 0*/
  	sds=my_ddot(m,ss,ss);	/*S . (R - G*alfa - S*beta) = 0*/
  	gds=my_ddot(m,gg,ss);
  	determ=gdg * sds - gds * gds + 1.e-15;
  	gdr=my_ddot(m,gg,rr);
  	sdr=my_ddot(m,ss,rr);
  	alfa=( sds * gdr - gds * sdr ) / determ;
  	beta=(-gds * gdr + gdg * sdr ) / determ;
  }

  for(i=0;i<n;i++) /* s = model step */
  	s[i]=alfa*g[i]+beta*s[i];

  for(i=0;i<m;i++) /* ss = conjugate */
  	ss[i]=alfa*gg[i]+beta*ss[i];
  	
  for(i=0;i<n;i++) /* update solution */
  	x[i]=x[i]+s[i];
  	
  for(i=0;i<m;i++) /* update residual */
  	rr[i]=rr[i]-ss[i];
}


void velinvww( int nt, float dt, float t0, 
			   float *x2, int nx, /*squared space axis and number of samples*/
			   float *data, 	  /*CMP gather*/
			   float *s, int ns,  /*velocity axis*/
			   float *mm, 			/*velocity spectrum*/
			   float *z2,			/*squared z axis*/
			   float *mask, 		/*mask*/
			   float rwt, float mwt, int niter, /*inversion parameter*/
			   int huber, float eps, int irls, int nstep, float srate)
{
  int iter, iiter, ndata, nmodl, i;
  float *dm, *sm, *wdm, *tmm, *rr, *dr, *sr, *wrr, *trr;
  
  dm=sf_floatalloc(nt*ns);
  sm=sf_floatalloc(nt*ns);
  wdm=sf_floatalloc(nt*ns);
  tmm=sf_floatalloc(nt*ns);
  
  rr=sf_floatalloc(nt*nx);
  dr=sf_floatalloc(nt*nx);
  sr=sf_floatalloc(nt*nx);
  wrr=sf_floatalloc(nt*nx);
  trr=sf_floatalloc(nt*nx);
   
  rwt = 0.5 * (rwt - 2.);   //! residual weight exponent (p-2)/2 for L-p norm

  ndata = nt*nx;
  nmodl = nt*ns;

  zero( nmodl, mm);
  
  copy( ndata, data,  rr);
  copy( ndata, rr,   trr);

  if(irls==0.0)			/*CGG inversion*/
  {
  	for(iter=0;iter<=niter;iter++)
  	{
  		powwt(iter,rr,ndata,rr,wrr,rwt,eps,huber,srate);
  		velxf(1,0,nt,dt,t0,x2,nx,wrr,mask,s,ns,dm,z2);
  		if(mwt!=0)
  			powwt(iter,dm,nmodl,mm,dm,mwt,eps,0,srate);
  		velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2);
  		cgstep(iter,nmodl,mm,dm,sm,ndata,rr,dr,sr);
  	}
  }else{				/*IRLS inversion*/
  		
  		if(mwt!=0)
  			mwt=0.5*(2.0-mwt);	/*model weight exponent (2-p)/2 for L-p norm*/
  		for(iter=0;iter<=niter;iter++)
  		{
  			velxf( 1, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, dm,z2);
  			velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, dm,z2);
  			cgstep( iter, nmodl, mm,dm,sm, ndata, rr,dr,sr);
  		}
  		for(iiter=1;iiter<=nstep;iiter++)
  		{
  			copy( nmodl, mm, tmm);
  			copy( ndata, rr, trr);
  			powwt( iiter, mm,nmodl, tmm,wdm, mwt,eps,0,srate);
  			velxf( 0, 0, nt,dt,t0, x2,nx, rr,mask, s,ns, wdm,z2);
  			for(i=0;i<ndata;i++)
  				rr[i]=data[i]-rr[i];
  			powwt( iiter, rr,ndata, trr,rr, rwt,eps,huber,srate);
  			for(iter=1;iter<=niter;iter++)
  			{
              powwt( iter, rr,ndata, trr,wrr, rwt,eps,huber,srate);       
              velxf( 1, 0, nt,dt,t0, x2,nx, wrr,mask, s,ns, dm,z2);
              powwt( iter, dm,nmodl, tmm, dm, mwt,eps,0,srate);
              powwt( iter, dm,nmodl, tmm,wdm, mwt,eps,0,srate);
              velxf( 0, 0, nt,dt,t0, x2,nx, dr,mask, s,ns, wdm,z2);
              powwt( iter, dr,ndata, trr, dr, rwt,eps,huber,srate);
              cgstep(iter, nmodl, mm,dm,sm, ndata, rr,dr,sr);
  			}
  			powwt( iter, mm,nmodl, tmm,mm, mwt,eps,0,srate);
  		}
  }

}



int main(int argc, char* argv[])
{
    int nt, ns, nh, it, is, ih;
    float dt, ds, dh, ot, os, oh;
    float *mask, *z2,z, *s, *oh2, *cmp, *intcmp, *vel;
    sf_file in, out, vout;
    
    sf_init(argc,argv);
    
    in = sf_input("in");
    out = sf_output("out");

	if(!sf_histint(in,"n1",&nt)) 	sf_error("Need nt");
	if(!sf_histfloat(in,"d1",&dt))	sf_error("Need dt");
	if(!sf_histfloat(in,"o1",&ot)) 	sf_error("Need ot");
	if(!sf_histint(in,"n2",&nh)) 	sf_error("Need nh");
	if(!sf_histfloat(in,"d2",&dh)) 	sf_error("Need dh");	
	if(!sf_histfloat(in,"o2",&oh)) 	sf_error("Need oh");
    if (!sf_getint("ns",&ns))   ns=nh;
    if (!sf_getfloat("ds",&ds)) sf_error("Need ds");
    if (!sf_getfloat("os",&os)) sf_error("Need os");
 
	int huber, irls, nstep, niter, savevel;
	float srate, rwt, mwt, eps, h;
	
	if (!sf_getint("huber",&huber))   huber=0;
	if (!sf_getint("irls",&irls))   irls=0;
	if (!sf_getint("nstep",&nstep))   nstep=1;
	if (!sf_getint("niter",&niter))   niter=20;
	if (!sf_getint("savevel",&savevel))   savevel=0;
	
    if (!sf_getfloat("rwt",&rwt)) rwt=0.0;	
    if (!sf_getfloat("mwt",&mwt)) mwt=0.0;	
    if (!sf_getfloat("epw",&eps)) eps=0.01;	
    if (!sf_getfloat("srate",&srate)) srate=0.01;	
	
	cmp=sf_floatalloc(nt*nh);		
	intcmp=sf_floatalloc(nt*nh);	
	vel=sf_floatalloc(nt*ns);			
	mask=sf_floatalloc(nh);			
	oh2=sf_floatalloc(nh);
	z2=sf_floatalloc(nt);			
	s=sf_floatalloc(ns);
  	
	for(ih=0;ih<nh;ih++)
	{	h=oh+dh*ih; oh2[ih]=h*h;} 
  
	for(it=0;it<nt;it++)
	{	z=ot+dt*it; z2[it]=z*z;}

	for(is=0;is<ns;is++)
		s[is]=os+ds*is;
	
	sf_floatread(cmp,nt*nh,in);
	
	for(ih=0;ih<nh;ih++)
	{	mask[ih]=0.0;
		for(it=0;it<nt;it++)
			mask[ih]=mask[ih]+cmp[it+ih*nt]*cmp[it+ih*nt];
	}

	sf_warning("niter=%d, ns=%d, ds=%g, os=%g",niter,ns,ds,os);
	sf_warning("rwt=%g, mwt=%g, savevel=%d",rwt,mwt,savevel);

  	velinvww( nt,dt,ot, oh2,nh,cmp, s, ns,vel, z2,mask,rwt,mwt,niter,huber,eps,irls,nstep,srate);

	if(savevel==1)
	{	
		vout=sf_output("velout");
		sf_putint(vout,"n1",nt);
		sf_putint(vout,"n2",ns);	
		sf_putfloat(vout,"d1",dt);
		sf_putfloat(vout,"d2",ds);
		sf_putfloat(vout,"o1",ot);
		sf_putfloat(vout,"o2",os);
		sf_floatwrite(vel,nt*ns,vout);
	}

	for(ih=0;ih<nh;ih++)
			mask[ih]=1.0;

  	velxf( 0,0, nt,dt,ot, oh2,nh,intcmp, mask, s,ns, vel, z2);

	for(ih=0;ih<nh;ih++)
	{	mask[ih]=0.0;
		for(it=0;it<nt;it++)
			mask[ih]=mask[ih]+cmp[it+ih*nt]*cmp[it+ih*nt];
	}
	
	
	for(ih=0;ih<nh;ih++)/*meaningless, right?*/
	{	
		if (mask[ih] == 0)
		  for(it=0;it<nt;it++)
			cmp[it+ih*nh]=intcmp[it+ih*nh];
	}
		
		
	sf_floatwrite(intcmp,nt*nh,out);

    exit(0);
}
