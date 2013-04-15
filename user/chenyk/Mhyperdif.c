/* Solving 1-D transportation equation using finite difference algorithm 
 \partial(u)/\partial(t)+a(x,t)\partial(u)/\partial(x), 0<t<=T, x->unlimited,
   u(x,0)=f(x).  */

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

#include<rsf.h>

float f (float x);
float f_true(float x, float t);
float axt(float x, float t);

int main(int argc, char* argv[])
{
   int j, nx, nt, n;
   float dt, dx; /* dt corresponds to k, dx corresponds to h */
   float *uo, *up, *u_true, *t, **a, k, h, ox; /*uo->u(n), up->u(n+1)*/
   bool iftrue, ifinitial;
   char *type;
   sf_file out, dtrue, dinit;
   sf_init(argc,argv);
   out=sf_output("out");
   
  if(!sf_getbool("wanttrue",&iftrue)) iftrue=false;
  /* [y/n] if want true solution. y: want, n: don't want. */

  if(!sf_getbool("wantinit",&ifinitial)) ifinitial=false;
  /* [y/n] if want initial value. y: want, n: don't want. */

  if(NULL==(type=sf_getstring("type"))) type="upwind";
  /* [upwind, friedrichs, wendroff] get the type for solving hyperbola partial differential equation, the default is upwind */
  
  if(!sf_getint("nt",&nt)) sf_error("No nt in command line!");
  /* number of temporal points */
  if(!sf_getfloat("dt",&dt)) sf_error("No dt in command line!");
  /* temporal sampling */

   if(!sf_getfloat("ox",&ox)) sf_error("No ox in command line!");
   /* spatial starting point */
   if(!sf_getint("nx",&nx)) sf_error("No nx in command line!");
   /* number of spatial points */
   if(!sf_getfloat("dx",&dx)) sf_error("No dx in command line!");
   /* spatial sampling */
  
   sf_putint(out,"n1",nx);
   sf_putfloat(out,"d1",dx);
   sf_putfloat(out,"o1",ox);


   uo=sf_floatalloc(nx);
   up=sf_floatalloc(nx);
   u_true=sf_floatalloc(nx);
   a=sf_floatalloc2(nt,nx);

   k=dt;h=dx;
   
   for(j=0;j<nx;j++)
	for(n=0;n<nt;n++)
		a[j][n]=axt(ox+j*dx,n*dt);

/* initialize it=0 */
   for(j=0;j<nx;j++)
	uo[j]=f(ox+dx*j);

   if(ifinitial==true) 
   {  
        dinit=sf_output("dinit");
	sf_putint(dinit,"n1",nx);
	sf_putfloat(dinit,"d1",dx);
	sf_putfloat(dinit,"o1",ox);
	sf_floatwrite(uo,nx,dinit);
   }

/* loop over it until finish computing up[nt-1] */
   for(n=0;n<nt-1;n++)  /* (nt-1) iterations */
	{
	up[0]=0;up[nx-1]=0;	
	for(j=1;j<nx-1;j++)
	{
		switch(type[0])
		{
		case 'u': if(a[j][n]>=0) {up[j]=-a[j][n]*k/h*(uo[j]-uo[j-1])+uo[j];} else {up[j]=-a[j][n]*k/h*(uo[j+1]-uo[j])+uo[j];}  break;
		case 'f': up[j]=-a[j][n]*k/2/h*(uo[j+1]-uo[j-1])+0.5*(uo[j-1]+uo[j+1]);   break;
		case 'w': up[j]=(a[j][n]*a[j][n]*k/2/h/h*(uo[j+1]-2*uo[j]+uo[j-1])-a[j][n]/2/h*(uo[j+1]-uo[j-1]))*k+uo[j];
		}
	}  
	t=uo;
	uo=up;
	up=t;
	}
   for(j=0;j<nx;j++)
	{
	u_true[j]=f_true(ox+j*dx,(nt-1)*dt);
	}

   if(iftrue==true) 
   {  
        dtrue=sf_output("dtrue");
	sf_putint(dtrue,"n1",nx);
	sf_putfloat(dtrue,"d1",dx);
	sf_putfloat(dtrue,"o1",ox);
	sf_floatwrite(u_true,nx,dtrue);
   }

   /* Here uo refers to the final result. */
   sf_floatwrite(uo, nx, out);

   exit(0);
}

float f(float x)
{
   float y;
   y=expf(-0.8*x*x);
   return y;
}

float f_true(float x, float t)
{
   float y;
   y=expf(-0.8*(x-t)*(x-t));
   return y;
}

float axt(float x, float t)
{
   float y;
   y=1;
   return y;
}
