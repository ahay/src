/* Solving 1-D heat equation using explicit finite difference 
 \partial(u)/\partial(t)=a^2\partial^2(u)/\partial(x^2), 0<x<l & t>0,
   u(0,t)=u(l,t)=0, t>0, 
   u(x,0)=f(x), 0<=x<=l.  */

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
*/

#include<stdio.h>
#include<rsf.h>

float f (float x);
float f_true(float x, float t);

int main(int argc, char* argv[])
{
   int j, nx, nt, n;
   float dt, dx, lambda; /* dt corresponds to k, dx corresponds to h */
   float *uo, *up, *u_true, *dif, *t; /*uo->u(n), up->u(n+1)*/
   sf_file out;
   sf_init(argc,argv);
   out=sf_output("out");
   
  if(!sf_getint("nt",&nt)) sf_error("No nt in command line!");
  /* number of temporal points */
  if(!sf_getfloat("dt",&dt)) sf_error("No dt in command line!");
  /* temporal sampling */

  if(!sf_getint("nx",&nx)) sf_error("No nx in command line!");
  /* number of spatial points */
  if(!sf_getfloat("dx",&dx)) sf_error("No dx in command line!");
  /* spatial sampling */
   
  sf_putint(out,"n1",nx);
  sf_putfloat(out,"d1",dx);
  sf_putfloat(out,"o1",0);
 
   lambda=dt/(dx*dx); /* correpsonding to a^2k/h^2 */

   uo=sf_floatalloc(nx);
   up=sf_floatalloc(nx);
   u_true=sf_floatalloc(nx);
   dif=sf_floatalloc(nx);

/* initialize it=0 */
   for(j=0;j<nx;j++)
	uo[j]=f(dx*j);

/* loop over it until finish computing up[nt-1] */
   for(n=0;n<nt-1;n++)  /* (nt-1) iterations */
	{
	up[0]=0;up[nx-1]=0;
	for(j=1;j<nx-1;j++)
		up[j]=lambda*uo[j-1]+(1-2*lambda)*uo[j]+lambda*uo[j+1];
	t=uo;
	uo=up;
	up=t;
	}
   for(j=0;j<nx;j++)
	{u_true[j]=f_true(j*dx,(nt-1)*dt); dif[j]=fabs(u_true[j]-uo[j]);}

   for(j=0;j<nx;j++)
   {
	sf_warning("%.1f   %.8f   %.8f   %.3e",j*dx, u_true[j], uo[j], dif[j]);
   }

   /* Here uo refers to the final result. */
   sf_floatwrite(uo, nx, out);
   exit(0);
}

float f(float x)
{
   float y;
   y=sinf(SF_PI*x);
   return y;
}

float f_true(float x, float t)
{
   float y;
   y=expf(-SF_PI*SF_PI*t)*sin(SF_PI*x);
   return y;
}
