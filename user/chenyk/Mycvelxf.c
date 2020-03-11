/* Velocity transform for generating velocity spectra and its corresponding hyperbolic modeling (C version) */
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


int main(int argc, char* argv[])
{
    int nt, nx, ns, it, ix, is;
    float dt, dx, ds, ot, ox, os;
    float *data, *modl, *mask, *x2, *z2, x,z, *s;
    bool adj;
    sf_file in, out;
    
    sf_init(argc,argv);
    
    in = sf_input("in");
    out = sf_output("out");

	if(!sf_getbool("adj",&adj)) adj=false;
	/* if implement the adjoint transform instead of the inverse transform */
	/*adj=n: 	data to velocity spectrum; adj=y: velocity spectrum to data; */

	if(!adj) 
	{
	if(!sf_histint(in,"n1",&nt)) 	sf_error("Need nt");
	if(!sf_histfloat(in,"d1",&dt))	sf_error("Need dt");
	if(!sf_histfloat(in,"o1",&ot)) 	sf_error("Need ot");
	if(!sf_histint(in,"n2",&ns)) 	sf_error("Need ns");
	if(!sf_histfloat(in,"d2",&ds)) 	sf_error("Need ds");	
	if(!sf_histfloat(in,"o2",&os)) 	sf_error("Need os");
    if (!sf_getint("nx",&nx)) 		sf_error("Need nx");
    if (!sf_getfloat("dx",&dx)) 	sf_error("Need dx");
    if (!sf_getfloat("ox",&ox)) 	sf_error("Need ox");
	sf_putint(out,"n1",nt);
	sf_putint(out,"n2",nx);	
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o1",ot);
	sf_putfloat(out,"o2",ox);		
	}else{
	if(!sf_histint(in,"n1",&nt))   sf_error("Need nt");
	if(!sf_histfloat(in,"d1",&dt)) sf_error("Need dt");
	if(!sf_histfloat(in,"o1",&ot)) sf_error("Need ot");
	if(!sf_histint(in,"n2",&nx))   sf_error("Need nx");
	if(!sf_histfloat(in,"d2",&dx)) sf_error("Need dx");
	if(!sf_histfloat(in,"o2",&ox)) sf_error("Need ox");
    if (!sf_getint("ns",&ns))   sf_error("Need ns");
    if (!sf_getfloat("ds",&ds)) sf_error("Need ds");
    if (!sf_getfloat("os",&os)) sf_error("Need os");
	sf_putint(out,"n1",nt);
	sf_putint(out,"n2",ns);	
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"d2",ds);
	sf_putfloat(out,"o1",ot);
	sf_putfloat(out,"o2",os);	
	}
	data=sf_floatalloc(nt*nx);
	modl=sf_floatalloc(nt*ns);
	s=sf_floatalloc(ns);
	x2=sf_floatalloc(nx);
	mask=sf_floatalloc(nx);
	z2=sf_floatalloc(nt);	
	if(!adj)
		sf_floatread(modl,nt*ns,in);
	else
		sf_floatread(data,nt*nx,in);

	for(ix=0;ix<nx;ix++)
		mask[ix]=1.0;
	for(ix=0;ix<nx;ix++)
	{	x=ox+dx*ix; x2[ix]=x*x;}
	for(it=0;it<nt;it++)
	{	z=ot+dt*it; z2[it]=z*z;}
	for(is=0;is<ns;is++)
		s[is]=os+ds*is;
		
	velxf( adj, 0,nt,dt,ot, x2,nx, data,mask, s,ns, modl,z2);

	if(!adj)
		sf_floatwrite(data,nt*nx,out);
	else
		sf_floatwrite(modl,nt*ns,out);
		

    exit(0);
}
