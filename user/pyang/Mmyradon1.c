/* Hyperbolic radon transform 
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int nt,nh,nv;
static float dt,*v,*h;

void hradon_init(int nt_, int nh_, int nv_, float dt_, float *v_, float *h_)
/*< Initialization of hyperbolic radon operator, borrowed from chenyk >*/
{
	nt=nt_;
	nh=nh_;
	nv=nv_;
	dt=dt_;
	v=v_;
	h=h_;
}

void hradon_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< hyperbolic radon operator, borrowed from chenyk
adjoint: m(tau,p) = sum_{ih=1}^{nh} d(tau, sqrt{tau^2+h[ih]^2/p^2)}
>*/
{   
	int itau, ih, iv, it;
	float tau, hov, t;
	
	sf_adjnull(adj, add, nx, ny, x, y);/* x=m; y=d; nx=nt*nv; ny=nt*nh; */ 
	
	for(ih=0; ih<nh; ih++)
	for(iv=0; iv<nv; iv++)
	{
		hov=h[ih]/v[iv];// h over v
		for(itau=0; itau<nt; itau++)
		{
			tau=itau*dt;
			t=sqrtf(tau*tau+hov*hov);
			it=(int)(t/dt);
			if(it<nt)
			{
				if (adj) x[iv*nt+itau]+=y[ih*nt+it];
				else	 y[ih*nt+it]+=x[iv*nt+itau]; 
			}	
		}
	}
}


int main(int argc, char* argv[])
{
	bool adj, inv, verb;
	int ix, iv, nv, nt, nx, niter;
	float ox, dx, dv, v0, dt, t0;
	float *vv, *xx, *mm, *dd;
	sf_file in, out, offset=NULL, vel=NULL;

    	sf_init(argc,argv);
	in = sf_input("in");	/* input data or radon  */
	out =sf_output("out");	/* output radon or data */

    	if (!sf_getbool("adj",&adj)) adj=false;
	/* if y, perform adjoint operation */
    	if (!sf_getbool("inv",&inv)) inv=true; 
	/* if y, perform inverse operation */
    	if (!sf_getbool ("verb",&verb)) verb=false;
	/* verbosity flag */

    	/* read input file parameters */
    	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    	if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    	if (adj||inv) { 
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

		/* specify slope axis */
		if (!sf_getint  ("nv",&nv)) sf_error("Need nv=");
		/* number of p values (if adj=y) */
		if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
		/* p sampling (if adj=y) */
		if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
		/* p origin (if adj=y) */
		if (inv && !sf_getint("niter",&niter)) niter=100;
		/* number of CGLS iterations */

		sf_putint(  out,"n2",nv);
		sf_putfloat(out,"d2",dv);
		sf_putfloat(out,"o2",v0);
    	} else { /* modeling */
		if (!sf_histint  (in,"n2",&nv)) sf_error("No n2= in input");
		/* number of ray parameter if input in radon domain */
		if (!sf_histfloat(in,"d2",&dv)) sf_error("No d2= in input");
		/* p sampling interval if input in radon domain */
		if (!sf_histfloat(in,"o2",&v0)) sf_error("No o2= in input");
		/* p origin if input in radon domain */
		if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
		/* number of offsets (if adj=n) */
	
		sf_putint(out,"n2",nx);
    	}

	vv=sf_floatalloc(nv);
	xx=sf_floatalloc(nx);
	dd=sf_floatalloc(nt*nx);
	mm=sf_floatalloc(nt*nv);

    	if (NULL != vel) {// velocity axis
		sf_floatread(vv,nv,vel);
		sf_fileclose(vel);
    	} else {
		for(iv=0; iv<nv; iv++) vv[iv]=v0+iv*dv;	
	}
	if (adj||inv) {// m(tau,p)=sum_{i=0}^{nx} d(t=sqrt(tau^2+(x_i/v_j)^2),x_i)
		sf_floatread(dd, nt*nx, in);

	    	if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
		/* data origin in x */
	    	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
		/* sampling interval in x */
	} else {// d(t,h)=sum_{i=0}^{np} m(tau=sqrt(t^2-(x_i/v_j)^2),p_i)
		sf_floatread(mm, nt*nv, in);
	    	if (!sf_getfloat("ox",&ox)) sf_error("Need ox=");
		/* x origin */
	    	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
		/* sampling interval in x */

	    	sf_putfloat(out,"o2",ox);
	    	sf_putfloat(out,"d2",dx);
	}
    	if (NULL != offset) {// offset axis
		sf_floatread(xx,nx,offset);
		sf_fileclose(offset);
    	} else {
		for(ix=0; ix<nx; ix++) xx[ix]=ox+ix*dx;
	}

	hradon_init(nt, nx, nv, dt, vv, xx);
	hradon_lop(adj, false, nt*nv, nt*nh, mm, dd);
	/* least-squares inversion for hyperbolic radon transform */
	if(inv) sf_solver(hradon_lop, sf_cgstep, nt*nv, nt*nx, mm, dd, niter,"x0", mm, "verb", verb, "end");

	free(vv);
	free(xx);
	free(mm);
	free(dd);

    	exit(0);
}

