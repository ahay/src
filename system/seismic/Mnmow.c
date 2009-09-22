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
  int nx, ny, nt, N, ix, iy, it, j, nw;
  float dx, dy, dt, x0, y0, t0, x, y, detM;
  float g1, g2, g3;
  float M11, M12, M13, M21, M22, M23, M31, M32, M33; 
  float M11i, M12i, M13i, M21i, M22i, M23i, M31i, M32i, M33i; 
  float *dT=NULL, *Mat=NULL, *w=NULL;

  /*Declare and initialize Madagascar files*/
  sf_file inp=NULL, out=NULL, GTG=NULL;

  sf_init(argc, argv);
  inp = sf_input("in");
  /*Axes: 1-->x, 2-->y, 3-->t.*/
  out = sf_output("out");
  GTG = sf_output("M");

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
  if (!sf_getint("nw",&nw)) nw=3;
  sf_putint(out,"n1",nw);
  sf_putint(out,"n2",nt);
  sf_putint(out,"n3",1);

  /*Output GTG matrix for troubleshooting.*/
  sf_putint(GTG,"n1",10);
  sf_putint(GTG,"n2",1);
  sf_putint(GTG,"n3",1);

  /*Allocate memory for model and data vectors*/
  N=nx*ny;
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
  for (iy=0; iy < ny; iy++){
    y=2*(y0+iy*dy);
    for (ix=0; ix < nx; ix++){
      x=2*(x0+ix*dx);
      M11=M11+x*x*x*x;
      M12=M12+x*x*y*y;
      M13=M13+x*x*x*y;
      M22=M22+y*y*y*y;
      M23=M23+x*y*y*y;
    }
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

  /*Output for troubleshooting*/
  Mat = sf_floatalloc(10);
  Mat[0]=M11i;
  Mat[1]=M12i;
  Mat[2]=M13i;
  Mat[3]=M21i;
  Mat[4]=M22i;
  Mat[5]=M23i;
  Mat[6]=M31i;
  Mat[7]=M32i;
  Mat[8]=M33i;
  Mat[9]=detM;
  sf_floatwrite(Mat,10,GTG);

  /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {
    sf_floatread(dT,N,inp);
    g1=0;
    g2=0;
    g3=0;
    j=0;
    for (iy=0; iy < ny; iy++){
      y=2*(y0+iy*dy);
      for (ix=0; ix < nx; ix++){
	x=2*(x0+ix*dx);
	/*g=GTd*/
	g1=g1+dT[j]*x*x;
	g2=g2+dT[j]*y*y;
	g3=g3+dT[j]*2*x*y;
	j++;
      }
    }

    /*Compute W model vector at each t_0 (slowness-squared values)*/
    w[0]=M11i*g1+M12i*g2+M13i*g3;
    w[1]=M21i*g1+M22i*g2+M23i*g3;
    w[2]=M31i*g1+M32i*g2+M33i*g3;

    /*Write 3 W values to output at every time coordinate*/
    sf_floatwrite(w,nw,out);
  }
  sf_close();
  exit(0);
}
