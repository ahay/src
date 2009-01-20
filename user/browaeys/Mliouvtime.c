/* Linear operator (centered FD stencil) for time evolution in phase space of 1D harmonic oscillator. */
/*
   Equation is v.dt/dx - c.x.dt/dv = 1 with c > 0 
   Analytical solution is t = 1/sqrt(c).arctan[sqrt(c).x/v] 
   for initial condition t=0 at x=0 and any nonzero v
*/
/*
  Copyright (C) 2009 University of Texas at Austin

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

#include <math.h>
#include <float.h>

#include <rsf.h>


int main(int argc, char* argv[])
{
    int np,nin;                    /* linear system matrix dimension */
    int i;                         /* matrix row counters */

    bool adj;                      /* adjoint */

    int nx,ix;                     /* number of grid points in x, counter */
    int nv,iv;                     /* number of grid points in v, counter */

    float dx,ox;                   /* increment in x, starting position */
    float dv,ov;                   /* increment in v, starting position */ 

    float signadj;                 /* sign for adjoint implementation */
    float fx,fv;                   /* factors in the system matrix */

    float c;                       /* harmonic oscillator stiffness constant */

    float *b, *t;

    sf_file in,out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, perform adjoint operation */

    if (!sf_getfloat("c",&c)) sf_error("Need c=");
    if (c <= 0.0) sf_error("Need strictly positive c");
    /* potential positive constant */

    if (!sf_getint("nx",&nx)) sf_error("Need nx=");
    /* number of x values. */
    if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    /* x sampling. */
    if (!sf_getfloat("ox",&ox)) ox = 0.0;
    /* x origin. */

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* number of v values. */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* v sampling. */
    if (!sf_getfloat("ov",&ov)) ov = dv;
    /* v origin. */

    np = nv*nx;

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nin)) sf_error("No n1= in input");
    if (nin != np) sf_error("Incorrect size of input vector");

    /* adjoint is a^t = -a */
    if (adj) {
	signadj = -1.;
    } else { 
	signadj = 1.;
    }

    /* memory allocations */
    b = sf_floatalloc(np);
    t = sf_floatalloc(np);
    sf_floatread(b,np,in);

 
    /* ------------------------------------ */
    /* SQUARE LINEAR SYSTEM MATRIX b = A.t  */
    /* ------------------------------------ */

    /* Linear system depends on stencils and indexing of t(ix,iv) vector */
    /* t^T = [t(1,1), ..., t(nx,1), ...., t(1,nv), ..., t(nx,nv)] */
    /* FD centered stencil */

    fv = signadj*0.5*c/dv;
    fx = signadj*0.5*ov/dx;

    b[0] = fx*t[1] - fv*ox*t[nx];

    for (i = 1; i < nx; i++) b[i] = fx*(t[i+1] - t[i-1]) - fv*(ox + ix*dx)*t[i+nx];

    for (iv = 1; iv < (nv-1); iv++) {
	fx = signadj*0.5*(ov + iv*dv)/dx;
	for (ix = 0; ix < nx; ix++){
	    i = iv*nx + ix;
	    b[i] = fv*(ox + ix*dx)*(t[i-nx]-t[i+nx]) + fx*(t[i+1]-t[i-1]);
	}
    }

    fx = signadj*0.5*(ov+(nv-1)*dv)/dx;
    for (ix = 0; ix < nx; ix++){
	i = (nv-1)*nx + ix;
	b[i] = fv*(ox + ix*dx)*t[i-nx] + fx*(t[i+1]-t[i-1]);
    }

    b[np-1] = fv*(ox + (nx-1)*dx)*t[np-1-nx] - fx*t[np-2];

    /* output */
    sf_floatwrite(b,np,out);

    exit(0);
}
