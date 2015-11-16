/* Fwd-Adj of 3D NMO GMA for iterative LS coefficient solve */
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
  bool adj;
  int nx, ny, nt, N, ix, iy, it, i, j, nw, iw;
  float dx, dy, dt, x0, y0, t0, x, y, t0f;
  float w1,w2,w3,A1,A2,A3,A4,A5,B1,B2,B3,C1,C2,C3,C4,C5;
  float *dT, *X, *w, *coeff, **t0sq;

  sf_file inp, out, gather, inicoef, t0sqf;
  sf_init(argc, argv);
  inp = sf_input("in");
  out = sf_output("out");
  gather = sf_input("gather");
  inicoef = sf_input("mod");
  t0sqf = sf_input("t0sq");

  /* Adjoint flag */
  if (!sf_getbool("adj",&adj)) adj=true;

  if (adj) {
    /* input dT vector delta(t^2-t0^2), output dw vector delta(Wx, Wxy,Wy,A_i,B_i,C_i)' */
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&ny)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&nt)) sf_error("No n3=");
    if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dy)) sf_error("No d2=");
    if (!sf_histfloat(inp,"d3",&dt)) sf_error("No d3=");
    if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"o2",&y0)) sf_error("No o2=");
    if (!sf_histfloat(inp,"o3",&t0)) sf_error("No o3=");

    if (!sf_getint("nw",&nw)) nw=16;
    /* 16 parameters of 3D GMA*/

    sf_putint(out,"n1",nw);
    sf_putint(out,"n2",nt);
    sf_putint(out,"n3",1);

  } else {
    /* input dw vector delta(Wx, Wxy,Wy,A_i,B_i,C_i), output dT vector delta(t^2-t0^2) */
    if (!sf_histint(inp,"n1",&nw)) sf_error("No n1=");

    if (!sf_histint(gather,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(gather,"n2",&ny)) sf_error("No n2=");
    if (!sf_histint(gather,"n3",&nt)) sf_error("No n3=");
    if (!sf_histfloat(gather,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(gather,"d2",&dy)) sf_error("No d2=");
    if (!sf_histfloat(gather,"d3",&dt)) sf_error("No d3=");
    if (!sf_histfloat(gather,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(gather,"o2",&y0)) sf_error("No o2=");
    if (!sf_histfloat(gather,"o3",&t0)) sf_error("No o3=");

  }

  N=nx*ny;
 
  w = sf_floatalloc(nw);
  X = sf_floatalloc(nw);
  dT = sf_floatalloc(N);
  coeff = sf_floatalloc(16);
  t0sq = sf_floatalloc2(nx,ny);
  
  sf_floatread(coeff,nw,inicoef);
  sf_floatread(t0sq[0],nx*ny,t0sqf);
    
    /* Initial coefficients*/
    w1 = coeff[0];w2 = coeff[1];w3 = coeff[2];
    A1 = coeff[3];A2 = coeff[4];A3 = coeff[5];A4 = coeff[6];A5 = coeff[7];
    B1 = coeff[8];B2 = coeff[9];B3 = coeff[10];
    C1 = coeff[11];C2 = coeff[12];C3 = coeff[13];C4 = coeff[14];C5 = coeff[15];

  /* Initial Filter */
  if (adj) {
    for(i=0; i < nw; i++) w[i]=0.0;
  } else {
    sf_floatread(w,nw,inp);
  }

  
  /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {

    /* Initial Data */
    if (adj) {
        sf_floatread(dT,N,inp);
        for(i=0; i < nw; i++) {
            w[i]=0.0;
        }
    } else {
          for(i=0; i < N; i++) dT[i]=0.0;
    }
    
    j=0; /*CHANGED FROM JUST ABOVE it LOOP!!!!*/
    /* Loop through x-y coordinates*/
    for (iy=0; iy < ny; iy++){    
        y = y0+iy*dy;
        for (ix=0; ix < nx; ix++){
	        x = x0+ix*dx;
            
            t0f = sqrt(t0sq[iy][ix]);
            
            /* Compute Derivatives */
	        X[0]=x*x; // Wx
	        X[1]=x*y; // Wxy
	        X[2]=y*y; // Wy
	        
	        X[3] = pow(x,4)/(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))); // A1

            X[4] = (pow(x,3)*y)/(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))); // A2

            X[5] = (pow(x,2)*pow(y,2))/(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))); // A3

            X[6] = (x*pow(y,3))/(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))); // A4

            X[7] = pow(y,4)/(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))); // A5

            X[8] = -(((A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4))*(pow(x,2) + (pow(t0f,2)*pow(x,2))/sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))))/pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // B1

            X[9] = -(((A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4))*(x*y + (pow(t0f,2)*x*y)/sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))))/pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // B2

            X[10] = -(((A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4))*(pow(y,2) + (pow(t0f,2)*pow(y,2))/sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))))/pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // B3

            X[11] = -(pow(x,4)*(A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4)))/(2.*sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))*pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // C1

            X[12] = -(pow(x,3)*y*(A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4)))/(2.*sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))*pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // C2

            X[13] = -(pow(x,2)*pow(y,2)*(A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4)))/(2.*sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))*pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // C3

             X[14] = -(x*pow(y,3)*(A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4)))/(2.*sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))*pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // C4

             X[15] = -(pow(y,4)*(A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4)))/(2.*sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2)))*pow(pow(t0f,2) + B1*pow(x,2) + B2*x*y + B3*pow(y,2) + sqrt(pow(t0f,4) + C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4) + 2*pow(t0f,2)*(B1*pow(x,2) + B2*x*y + B3*pow(y,2))),2)); // C5

	/* Loop through slowness parameters*/
	        for (iw=0; iw < nw; iw++){

	/* Solve */
	              if (adj) {
	                w[iw] += X[iw]*dT[j];
	              } else {
	                dT[j] += X[iw]*w[iw];
	            }
	        }
	/*x-y loops corresponds to loop 0:N-1*/
	    j++;
      }
    }
    j=0;
    if (adj) sf_floatwrite(w,nw,out);
    if (!adj) sf_floatwrite(dT,N,out);
  }

  exit(0);
}
