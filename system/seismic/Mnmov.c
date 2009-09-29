/*Least-squares fitting of t^2-t_0^2 surfaces for isotropic V_{nmo}.*/

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

#include <rsf.h>

int main(int argc, char* argv[])
{
  /*Declare variables*/
  int nx, ny, nt, N, ix, iy, it, j, nw;
  float dx, dy, dt, x0, y0, t0, x, y, h;
  float g1;
  float M11, M11i;
  float *dT=NULL, *w=NULL;

  /*Declare and initialize Madagascar files*/
  sf_file inp=NULL, out=NULL;

  sf_init(argc, argv);
  inp = sf_input("in");
  /*Axes: 1-->x, 2-->y, 3-->t.*/
  out = sf_output("out");

  /*Read in input file dimensions*/
  /* input dT vector (t^2-t0^2), output w vector (Wx, Wy, Wxy)' */
  if (!sf_histint(inp,"n1",&nx)) sf_error("No n1=");
  if (!sf_histint(inp,"n2",&ny)) sf_error("No n2=");
  if (!sf_histint(inp,"n3",&nt)) sf_error("No n3=");
  if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1=");
  if (!sf_histfloat(inp,"d2",&dy)) sf_error("No d2=");
  if (!sf_histfloat(inp,"d3",&dt)) sf_error("No d3=");
  if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1=");
  if (!sf_histfloat(inp,"o2",&y0)) sf_error("No o2=");
  if (!sf_histfloat(inp,"o3",&t0)) sf_error("No o3=");

  /*Set size of output file (w is model vector at each time-slice)*/
  nw=1;
  sf_putint(out,"n1",nw);
  sf_putint(out,"n2",nt);
  sf_putint(out,"n3",1);

  /*Allocate memory for model and data vectors*/
  N=nx*ny;
  w = sf_floatalloc(nw);
  dT = sf_floatalloc(N);

  /*M=GTG (In this case size=1x1)*/
  M11=0;

  for (iy=0; iy < ny; iy++){
    y=(y0+iy*dy);
    for (ix=0; ix < nx; ix++){
      x=(x0+ix*dx);
      h=sqrt(x*x+y*y);
      M11=M11+h*h*h*h;
    }
  }

  /*Minv: Inverse of 1x1 is just reciprocal.*/
  M11i=1.0/M11;

   /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {
    sf_floatread(dT,N,inp);
    g1=0;
    j=0;
    for (iy=0; iy < ny; iy++){
      y=(y0+iy*dy);
      for (ix=0; ix < nx; ix++){
	x=(x0+ix*dx);
	h=sqrt(x*x+y*y);
	/*g=GTd*/
	g1=g1+h*h*dT[j];
	j++;
      }
    }

    /*Computed value is actually slowness; take reciprocal to get velocity.*/
    w[0]=1.0/(M11i*g1);

    /*Output velocity at each time value.*/
    sf_floatwrite(w,nw,out);
  }

  exit(0);
}
