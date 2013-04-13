/* A forward interpolation using Lagrange method. 
Specify ox= dx= nx=
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
*/

#include<rsf.h>
static float *dd, *xd;
static int N;

void lagrange_init(int n);
void lagrange_apply( float *xx, float *yy, int nx);
float L( float xx);
void lagrange_close(void);

int main(int argc, char* argv[])
{
   int i;
   float d1,o1,dx,ox;
   int n1,nx;
   float *x,*data;
   sf_file in, out;
   sf_init(argc,argv);
   
   in=sf_input("in");
   out=sf_output("out");

   if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input!");
   if(!sf_histfloat(in,"d1",&d1)) sf_error("No d1 in input!");
   if(!sf_histfloat(in,"o1",&o1)) sf_error("No o1 in input!");
   if(!sf_getfloat("ox",&ox)) ox=o1;
   if(!sf_getfloat("dx",&dx)) dx=d1;
   if(!sf_getint("nx",&nx)) nx=n1;

   lagrange_init(n1);

   for(i=0;i<n1;i++)
	xd[i]=o1+i*d1;
   sf_floatread(dd,n1,in);

   x=sf_floatalloc(nx);
   data=sf_floatalloc(nx);
   for(i=0;i<nx;i++)
	x[i]=ox+dx*i;
 
   lagrange_apply(x,data,nx);
   sf_putint(out,"n1",nx);
   sf_putfloat(out,"d1",dx);
   sf_putfloat(out,"o1",ox);
   sf_floatwrite(data,nx,out);
   
   lagrange_close();
   free(x);
   free(data);
   exit(0);
}

void lagrange_init(int n)
{
dd=sf_floatalloc(n);
xd=sf_floatalloc(n);
N=n;
}


void lagrange_apply( float *xx, float *yy, int nx)
{
   int i;
   for(i=0;i<nx;i++)
	yy[i]=L(xx[i]);
}

float L( float xx)
{  
   int i,j;
   float yy, t; 
   yy=0;
   for(i=0;i<N;i++)
   { 
        t=1;
	for(j=0;j<N;j++)
	    {if(j!=i) t*=(xx-xd[j])/(xd[i]-xd[j]);}
	yy+=t*dd[i];
   }
   return yy;
}

void lagrange_close(void)
{
free(dd);
free(xd);
}

