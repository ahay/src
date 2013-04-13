/* Solving 1-D heat equation using implicit finite difference 
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

static int i, n;
static float *a, *b, *c;

float f (float x);
float f_true(float x, float t);
void trid_init(int N/* matrix size */, float *A/*input matrix*/);
void trid_solve( float *d /* in - right-hand side */, float *x /* out - solution */);
void trid_close(void);

int main(int argc, char* argv[])
{
   int j, nx, nt, n;
   float dt, dx, lambda; /* dt corresponds to k, dx corresponds to h */
   float *uo, *up, *u_true, *dif, *t, *A; /*uo->u(n), up->u(n+1)*/
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
   A=sf_floatalloc(3*(nx-2));

   for(j=0;j<nx-1;j++)
	{A[j*3]=-lambda;A[j*3+1]=1+2*lambda;A[j*3+2]=-lambda;}

   /*initialize tridiagonal matrix*/
   trid_init(nx-2,A);

   /* initialize it=0 */
   for(j=0;j<nx;j++)
	uo[j]=f(dx*j);

   /* loop over it until finish computing up[nt-1] */
   for(n=0;n<nt-1;n++)  /* (nt-1) iterations */
	{
	up[0]=0;up[nx-1]=0;
        trid_solve(uo+1,up+1);
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

   /* release allocated memory */
   free(uo);
   free(up);
   free(u_true);
   free(dif);
   free(A);
   trid_close();
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

void trid_init(int N/* matrix size */,
		       float *A/*input matrix*/)
/*< initialize >*/
/*|ac   |    */		/*alpha		    *//* 1 gamma 	*/
/*|bac  |    */		/* beta alpha	    *//*     1 gamma   	*/
/*| bac |x=d */  /*->*/ /*       beta alpha *//*        1 gamma	*/
/*|  ...|    */ 	/*	       beta *//*	    1	*/
/*|   ba|    */		/*		    *//*	  	*/
{
    n=N;
    a=sf_floatalloc(n);
    b=sf_floatalloc(n);
    c=sf_floatalloc(n);
    for(i=0;i<n;i++)
	{a[i]=A[i*3+1];}

    for(i=1;i<n;i++)
	{b[i]=A[i*3];}
    for(i=0;i<n-1;i++)
	{c[i]=A[i*3+2];}	
}

void trid_solve( float *d /* in - right-hand side */, 
			   float *x /* out - solution */)
/*< invert the matrix >*/
{
    float *alpha, *beta, *gamma, *y;
    alpha=sf_floatalloc(n);
    beta=sf_floatalloc(n);
    gamma=sf_floatalloc(n);
    y=sf_floatalloc(n);

    /* LU decomposition */
    gamma[0]=c[0]/(a[0]+0.000000000000001);
    alpha[0]=a[0];
    for(i=1;i<n-1;i++)
	{alpha[i]=a[i]-b[i]*gamma[i-1];
    	gamma[i]=c[i]/(alpha[i]+0.000000000000001);}
    alpha[n-1]=a[n-1]-b[n-1]*gamma[n-2];
    for(i=1;i<n;i++)beta[i]=b[i];

    /* Solving Ly=d */    
    y[0]=d[0]/(alpha[0]+0.000000000000001);
    for(i=1;i<n;i++)
	y[i]=(d[i]-b[i]*y[i-1])/(alpha[i]+0.000000000000001);

    /* Solving Ux=y */ 
    x[n-1]=y[n-1];
    for(i=n-2;i>=0;i--)
	x[i]=y[i]-gamma[i]*x[i+1];

}

void trid_close(void)
/*< free allocated memory >*/
{
    free(a);
    free(b);
    free(c);
}


