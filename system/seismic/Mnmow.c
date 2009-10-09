/*Least-squares fitting of t^2-t_0^2 surfaces for elliptical slowness matrix, W.*/

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
  bool half;
  int nx, ny, nt, N, ix, iy, it, j, nw, nhx, nhy, nhdim;
  float dx, dy, dt, x0, y0, t0, x, y, detM, halfscale;
  float g1, g2, g3;
  float M11, M12, M13, M21, M22, M23, M31, M32, M33; 
  float M11i, M12i, M13i, M21i, M22i, M23i, M31i, M32i, M33i; 
  float *dT, *w, *offx, *offy;

  /*Declare and initialize Madagascar files*/
  sf_file inp, out, offset;

  sf_init(argc, argv);
  inp = sf_input("in");
  /*Axes: 1-->x, 2-->y, 3-->t.*/
  out = sf_output("out");
 
  /* input dT vector (t^2-t0^2), output w vector (Wx, Wy, Wxy)' */
  if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");
  if (!sf_histint(inp,"n2",&ny)) sf_error("No n2= in input");
  if (!sf_histint(inp,"n3",&nt)) sf_error("No n3= in input");
  if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1= in input");
  if (!sf_histfloat(inp,"d2",&dy)) sf_error("No d2= in input");
  if (!sf_histfloat(inp,"d3",&dt)) sf_error("No d3= in input");
  if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1= in input");
  if (!sf_histfloat(inp,"o2",&y0)) sf_error("No o2= in input");
  if (!sf_histfloat(inp,"o3",&t0)) sf_error("No o3= in input");

  halfscale=1.0;
  if (!sf_getbool("half",&half)) half=true;
  /* if y, the second axis is half-offset instead of full offset */
  if (half) halfscale=2.0;

  /*Setup offset vector*/
  if (NULL != sf_getstring("offset")) {
    /*If offset file is provided, it must be of the form:*/
    /*Axes: 1-->x, 2-->y, 3-->dim.*/
    /*All points in offset file with indices (i,j,0) are x-offset values*/
    /*All points in offset file with indices (i,j,1) are y-offset values*/
    /*The offset file should therefore be of size nx*ny*2*/
    offset = sf_input("offset");
    if (!sf_histint(offset,"n1",&nhx)) sf_error("No n1= in offset");
    if (!sf_histint(offset,"n2",&nhy)) sf_error("No n2= in offset");
    if (!sf_histint(offset,"n3",&nhdim)) sf_error("No n3= in offset");
    N = nhx*nhy;
    if (N != nx*ny) sf_error("Wrong dimensions in offset");
    if (nhdim != 2) sf_error("Wrong dimensions in offset");
    
    offx = sf_floatalloc(N);
    offy = sf_floatalloc(N);
    sf_floatread (offx,N,offset);
    sf_floatread (offy,N,offset);
    sf_fileclose(offset);
  } else {
     
    N = nx*ny;
    offx = sf_floatalloc(N);
    offy = sf_floatalloc(N);
    j=0;
    for (iy=0; iy < ny; iy++){
      y=2*(y0+iy*dy);
      for (ix=0; ix < nx; ix++){
	offx[j] = x0 + ix*dx; 
	offy[j] = y0 + iy*dy;
	j=j+1;
      }
    
      offset = NULL;
    }
  }

  /*Set size of output file (w is model vector at each time-slice)*/
  if (!sf_getint("nw",&nw)) nw=3;
  sf_putint(out,"n1",nw);
  sf_putint(out,"n2",nt);
  sf_putint(out,"n3",1);

    /*Allocate memory for model and data vectors*/
   w = sf_floatalloc(nw);
  dT = sf_floatalloc(N);

  /*M=GTG (In W case, size(GTG)=3x3)*/
  M11=0;
  M12=0;
  M13=0;
  M21=0;
  M22=0;
  M23=0;
  M31=0;
  M32=0;
  M33=0;

  for (j=0; j < N; j++){
    x=offx[j]*halfscale;
    y=offy[j]*halfscale;
    M11=M11+x*x*x*x;
    M12=M12+x*x*y*y;
    M13=M13+x*x*x*y;
    M22=M22+y*y*y*y;
    M23=M23+x*y*y*y;
  }


  M13=2*M13;
  M21=M12;
  M23=2*M23;
  M31=M13;
  M32=M23;
  M33=4*M12;

  /*Minv: Inverse of M using cascaded determinants*/
  detM=1.0/((M11*M22*M33)-(M11*M23*M32)-(M12*M21*M33)+(M12*M23*M31)+(M13*M21*M32)-(M13*M22*M31));
  M11i=detM*(M33*M22-M23*M32);
  M12i=detM*(M32*M13-M33*M12);
  M13i=detM*(M23*M12-M22*M13);
  M21i=detM*(M31*M23-M33*M21);
  M22i=detM*(M33*M11-M31*M13);
  M23i=detM*(M21*M13-M23*M11);
  M31i=detM*(M32*M21-M31*M22);
  M32i=detM*(M31*M12-M32*M11);
  M33i=detM*(M22*M11-M21*M12);


  /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {
    sf_floatread(dT,N,inp);
    g1=0;
    g2=0;
    g3=0;
    j=0;
    for (j=0; j < N; j++){
      x=offx[j]*halfscale;
      y=offy[j]*halfscale;
      /*g=GTd*/
      g1=g1+dT[j]*x*x;
      g2=g2+dT[j]*y*y;
      g3=g3+dT[j]*2*x*y;
    }


    /*Compute W model vector at each t_0 (slowness-squared values)*/
    w[0]=M11i*g1+M12i*g2+M13i*g3;
    w[1]=M21i*g1+M22i*g2+M23i*g3;
    w[2]=M31i*g1+M32i*g2+M33i*g3;

    /*Write 3 W values to output at every time coordinate*/
    sf_floatwrite(w,nw,out);
  }

  exit(0);
}
