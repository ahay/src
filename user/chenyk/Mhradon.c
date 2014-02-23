/* Time domain high-resolution hyperbolic Radon transform. 
m(tau,p) = \sum_{ih=1}^{nh} d(tau=\sqrt{tau^2+h[ih]^2/p^2),h}
inv=true do inverse
adj=true do adjoint
inv=false && adj=false do forward
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

#include <math.h>

#include <rsf.h>
#include "radonutil.h"

int main (int argc, char *argv[])
{
	int i3,i;						/* index */
	int n2,n3;						/* second and third axis size */
	bool inv,adj;						/* inverse and adjoint flags */
	bool solver;						/* if use Madagascar bigsolver */
	int verb;						/* if output verbosity */
	int N_internal,N_external; 				/* internal and external iterations in PCG  */
	int nh,nt,nv;						/* number of points in t,h,v axises */
	float *data=NULL;					/* t-x domain */
	float *model=NULL;					/* tau-p domain */
	float *misfit=NULL;					/* esitimated and real data domain misfit */	
	float *v=NULL;						/* velocity vector */
	float *h=NULL;						/* offset vector */
	float t0,dt,h0,dh,v0,dv;				/* origin and increments for t,h,v axises. */
	sf_file in,out;						/* standard input and output files */
	sf_file vel,offset;				/* Input velocity and offset vector files*/

/***************************************************/
/*	Initialization				   */
/***************************************************/
	sf_init(argc,argv);
	in= sf_input("in");
	out=sf_output("out");

/***************************************************/
/* Getting and putting dimensions 	 		*/
/***************************************************/	
	if(!sf_histint(in,"n1",&nt)) sf_error("No n1 in input");
	if(!sf_histfloat(in,"d1",&dt)) dt=0.004;
	if(!sf_histfloat(in,"o1",&t0)) t0=0;		
	if(!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");
	if(!sf_histint(in,"n3",&n3)) n3=1;
			
	if(!sf_getbool("inv",&inv)) inv=true;
	/* if implement the inverse transform */
	
	if(!sf_getbool("adj",&adj)) adj=false;
	/* if implement the adjoint transform instead of the inverse transform */

	if(!sf_getbool("solver",&solver)) solver=false;
	/* if use Madagascar bigsolver, default is not */

/***************************************************/
/* Setting parameters for inv and for respectively */
/***************************************************/	
	if(inv || adj ) /*inverse or adjoint*/
	{
		nh=n2;	
		if(NULL!=sf_getstring("vel"))
		{
			vel=sf_input("vel");
			if(!sf_histint(vel,"n1",&nv)) sf_error("No n1 in velocity file");
			v=sf_floatalloc(nv);
			sf_floatread(v,nv,vel);	
		}
		else
		{
			if(!sf_getfloat("v0",&v0)) sf_error("No v0 specified");
			if(!sf_getfloat("dv",&dv)) sf_error("No dv specified");		
			if(!sf_getint("nv",&nv)) sf_error("No nv specified");	
			v=sf_floatalloc(nv);
			for(i=0;i<nv;i++)
				v[i]=v0+i*dv;	
			sf_putfloat(out,"d2",dv);
			sf_putfloat(out,"o2",v0);		
		}
		
		if(NULL!=sf_getstring("offset"))
		{
			offset=sf_input("offset");
			if(!sf_histint(offset,"n1",&nh)) sf_error("No n1 in offset file");
			h=sf_floatalloc(nh);
			sf_floatread(h,nh,offset);	
		}
		else
		{
			if(!sf_histfloat(in,"d2",&dh)) sf_error("No d2 in Input");
			if(!sf_histfloat(in,"o2",&h0)) sf_error("No o2 in Input");		

			h=sf_floatalloc(nh);
			for(i=0;i<nh;i++)
				h[i]=h0+i*dh;			
		}	
				
		sf_putint(out,"n2",nv);
		if(!sf_getint("N1",&N_internal)) N_internal=10;
		/* CG Iterations (Internal loop) */
		if(!sf_getint("N2",&N_external)) N_external=3; 
		/* Update of weights for the sparse solution, N1 = 1 LS , N2 > 3 for High Res (Sparse) solution */
		if(!sf_getint("verb",&verb)) verb=0;
		/* If output the debugging process */
		
		model=sf_floatalloc(nt*nv);	
		data=sf_floatalloc(nt*nh);		
	}
	else /*forward*/
	{	
		nv=n2;
		if(NULL!=sf_getstring("offset"))
		{
			offset=sf_input("offset");
			if(!sf_histint(offset,"n1",&nh)) sf_error("No n1 in offset file");
			h=sf_floatalloc(nh);
			sf_floatread(h,nh,offset);	
		}
		else
		{
			if(!sf_getfloat("h0",&h0)) sf_error("No h0 specified");
			if(!sf_getfloat("dh",&dh)) sf_error("No dh specified");		
			if(!sf_getint("nh",&nh)) sf_error("No nh specified");	
			h=sf_floatalloc(nh);
			for(i=0;i<nh;i++)
				h[i]=h0+i*dh;	
			sf_putfloat(out,"d2",dh);
			sf_putfloat(out,"o2",h0);
		}	
		
			if(!sf_histfloat(in,"d2",&dv)) sf_error("No d2 in Input");
			if(!sf_histfloat(in,"o2",&v0)) sf_error("No o2 in Input");		

			v=sf_floatalloc(nv);
			for(i=0;i<nv;i++)
				v[i]=v0+i*dv;			

		
		sf_putint(out,"n2",nh);
		model=sf_floatalloc(nt*nv);
		data=sf_floatalloc(nt*nh);
	}

/***************************************************/
/* Initialization for hyperbolic Radon transform 	 */
/***************************************************/
	hradon_init(nt, nh, nv, dt, v, h);	
	
/***************************************************/
/* Forward or Inverse hyperbolic Radon transform 	 */
/***************************************************/	
	for(i3=0;i3<n3;i3++)
	{
		if(inv) 	/* do high resolution inverse hyperbolic radon (m ~= F^{-1} d) */
		{
			sf_floatread(data,nt*nh,in);
			if(!solver) /* use my own solver */
			{
			hradon_pcg(hradon, data, model, N_internal,N_external, verb, misfit);
			}else{ /* use Madagascar big solver */
			int ix,iter;
    			float *w, *p, eps=0.01;
			float *error; error=sf_floatalloc(N_internal);
	
    			sf_cdstep_init();

    			w = sf_floatalloc(nt*nv);
    			p = sf_floatalloc(nt*nv);

    			for (ix=0; ix < nt*nv; ix++) {w[ix] = 1.f;model[ix]=0.0;p[ix]=0.0;}

    			for (iter = 0; iter < N_external; iter++) {
        			sf_solver_prec (hradon, sf_cgstep, sf_copy_lop, nt*nv, nt*nv, nt*nh, model, data, 
                        N_internal, eps, "err",error , "verb", verb, "mwt", w, "xp", p, "end");
        			for (ix=0; ix < nt*nv; ix++) w[ix] = fabsf(p[ix]);
    			}
    			sf_cdstep_close();
			}

		}
		else if (!adj) /*  do forward hyperbolic radon transform (d=Fm) */
			{
			sf_floatread(model,nt*nv,in);	
			hradon(false,false,nt*nv,nt*nh,model,data);	
			}
			else		  /* do adjoint hyperbolic radon tranform (m=F^T d) */
			{
			sf_floatread(data,nt*nh,in);
			hradon(true,false,nt*nv,nt*nh,model,data);	
			}	
	 	
		if(inv || adj)
		{ sf_floatwrite(model,nt*nv,out);} //model must not be zero
		else 	{sf_floatwrite(data,nt*nh,out);}
	}
	exit(0);
}

