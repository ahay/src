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
	float tau, h_v, t;
	
	sf_adjnull(adj, add, nx, ny, x, y);/* x=m; y=d; nx=nt*nv; ny=nt*nh; */ 
	
	for(ih=0; ih<nh; ih++)
	for(iv=0; iv<nv; iv++)
	{
		h_v=h[ih]/v[iv];
		for(itau=0; itau<nt; itau++)
		{
			tau=itau*dt;
			t=sqrtf(tau*tau+h_v*h_v);
			it=(int)(t/dt)+1;
			if(it<=nt)
			{
				if (adj) x[iv*nt+itau]+=y[ih*nt+it];
				else	 y[ih*nt+it]+=x[iv*nt+itau]; 
			}	
		}
	}
}

void inv_hradon(float *mod, float *dat, float *v, float *h, int nt, 
	int nh, int nv, float dt, int niter, bool verb)
/* least-squares inversion for hyperbolic radon transform */
{
	hradon_init(nt, nh, nv, dt, v, h);
   	sf_solver(hradon_lop, sf_cgstep, nt*nv, nt*nh, mod, dat, niter, "verb", verb, "end");
}


int main(int argc, char* argv[])
{
	bool adj, inv, verb;
	int i, j, k, np, nt, nx;
	float dp, p0, dt, t0;
	sf_file in, out;

    	sf_init(argc,argv);
	in = sf_input("in");	/* veloctiy model */
	out =sf_output("out");	/* shot records */

    	if (!sf_getbool("adj",&adj)) adj=true;
	/* if y, perform adjoint operation */
    	if (!sf_getbool("inv",&inv)) inv=adj; 
	/* if y, perform inverse operation */
    	if (!sf_getbool ("verb",&verb)) verb=false;
	/* verbosity flag */

    	/* read input file parameters */
    	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    	if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    	if (adj) { 
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

		/* specify slope axis */
		if (!sf_getint  ("np",&np)) sf_error("Need np=");
		/* number of p values (if adj=y) */
		if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
		/* p sampling (if adj=y) */
		if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
		/* p origin (if adj=y) */

		sf_putint(  out,"n2",np);
		sf_putfloat(out,"d2",dp);
		sf_putfloat(out,"o2",p0);
    	} else { /* modeling */
		if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
		if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
		if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");
		if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
		/* number of offsets (if adj=n) */
	
		sf_putint(out,"n2",nx);
    	}

    	exit(0);
}

