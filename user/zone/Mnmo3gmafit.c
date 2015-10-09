/* 3D NMO GMA  linearized operator preparation for lsfit*/
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
  int nx, ny, nt, npara, nc;
  float x, y, dx, dy, dt, x0, y0, tini, t0;
  float **t0sq,**time, *coeff, ***dtime; 
  float w1,w2,w3,A1,A2,A3,A4,A5,B1,B2,B3,C1,C2,C3,C4,C5,A,B,C;
  int i,j,k,test;
  int count = 0;
  bool verb;

  sf_file inp, out, fit ,inicoef;

  sf_init(argc, argv);
    inp = sf_input("in");
    inicoef = sf_input("coef");
    out = sf_output("out");
    fit = sf_output("fit");

    /* input vector t0^2' */
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&ny)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&nt)) sf_error("No n3=");
    if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dy)) sf_error("No d2=");
    if (!sf_histfloat(inp,"d3",&dt)) sf_error("No d3=");
    if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"o2",&y0)) sf_error("No o2=");
    if (!sf_histfloat(inp,"o3",&tini)) sf_error("No o2=");

    /*Number of fitting parameters*/
    npara=16;
    
    /*Verbal flag*/
    if (!sf_getbool("verb",&verb)) verb=false;
    
    /*Memory allocation*/
    coeff = sf_floatalloc(npara);
    t0sq = sf_floatalloc2(nx,ny);
    time = sf_floatalloc2(nx,ny);
    dtime = sf_floatalloc3(nx,ny,npara);
    
    /*Output dimension*/
    if (!sf_histint(inicoef,"n1",&nc) || nc != npara) 
	sf_error("Need n1=%d in inicoef",npara);
    
    /*Shift the third dimension to 4th to insert coefficients*/
     sf_shiftdim(inp, fit, 3);
     sf_putint(fit,"n3",npara);
     sf_putint(fit,"d3",1);
     sf_putint(fit,"o3",0);
    
    /* Loop over time slices */
      for(k=0; k < nt; k++) {
        /*Read intial parameters*/
        sf_floatread(coeff,npara,inicoef);
        w1 = coeff[0];w2 = coeff[1];w3 = coeff[2];
        A1 = coeff[3];A2 = coeff[4];A3 = coeff[5];A4 = coeff[6];A5 = coeff[7];
        B1 = coeff[8];B2 = coeff[9];B3 = coeff[10];
        C1 = coeff[11];C2 = coeff[12];C3 = coeff[13];C4 = coeff[14];C5 = coeff[15];
        sf_floatread(t0sq[0],nx*ny,inp);
        

        /* Loops for each x and y*/
            for(j=0;j<ny;j++){
                y = y0 + j*dy;
                
                for(i=0;i<nx;i++){
                    x = x0 + i*dx;
                    t0 = sqrt(t0sq[j][i]);
                    
                    /* Avoid dividing by zero*/
                    if (x!=0 || y!=0 || t0!=0) {
                    
                    A = A1*pow(x,4) + A2*pow(x,3)*y + A3*pow(x,2)*pow(y,2) + A4*x*pow(y,3) + A5*pow(y,4);
                    B = B1*pow(x,2) + B2*x*y + B3*pow(y,2);
                    C = C1*pow(x,4) + C2*pow(x,3)*y + C3*pow(x,2)*pow(y,2) + C4*x*pow(y,3) + C5*pow(y,4);
                    
                    /* Compute traveltime (t^2)*/
                    time[j][i] = pow(t0,2) + w1*pow(x,2) + w2*x*y + w3*pow(y,2) + A/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B));
                    /* Compute Derivatives */
                    dtime[0][j][i] = pow(x,2); // Wx
                    dtime[1][j][i] = x*y; // Wxy
                    dtime[2][j][i] = pow(y,2); //Wy

                    dtime[3][j][i] = pow(x,4)/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)); // A1

                    dtime[4][j][i] = (pow(x,3)*y)/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)); // A2

                    dtime[5][j][i] = (pow(x,2)*pow(y,2))/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)); // A3

                    dtime[6][j][i] = (x*pow(y,3))/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)); // A4

                    dtime[7][j][i] = pow(y,4)/(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)); // A5

                    dtime[8][j][i] = (-A*x*x*(1+pow(t0,2)/(sqrt(pow(t0,4) + C + 2*pow(t0,2)*B))))/pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2); // B1

                    dtime[9][j][i] = (-A*x*y*(1+pow(t0,2)/(sqrt(pow(t0,4) + C + 2*pow(t0,2)*B))))/pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2); // B2

                    dtime[10][j][i] = (-A*y*y*(1+pow(t0,2)/(sqrt(pow(t0,4) + C + 2*pow(t0,2)*B))))/pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2); // B3
                    
                    dtime[11][j][i] = (-pow(x,4)*A)/(2*sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)*pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2)); // C1

                    dtime[12][j][i] = (-pow(x,3)*y*A)/(2*sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)*pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2)); // C2

                    dtime[13][j][i] = (-pow(x,2)*pow(y,2)*A)/(2*sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)*pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2)); // C3

                    dtime[14][j][i] = (-x*pow(y,3)*A)/(2*sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)*pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2)); // C4

                    dtime[15][j][i] = (-pow(y,4)*A)/(2*sqrt(pow(t0,4) + C + 2*pow(t0,2)*B)*pow(pow(t0,2) + B + sqrt(pow(t0,4) + C + 2*pow(t0,2)*B),2)); // C5
                    
                    for (test=0;test<16;test++) {
                    
                            if(verb && isnan(dtime[test][j][i]) != 0) {
                                sf_warning("The sqrt is NaN at dtime %d k %d j %d i %d",test,k+1,j,i);
                                sf_warning("x %f y %f B %f %f %f C %f %f %f %f %f",x,y,B1,B2,B3,C1,C2,C3,C4,C5);
                                sf_warning("pow(t0,4) + C + 2*pow(t0,2)*B : %f \n",pow(t0,4) + C + 2*pow(t0,2)*B);
                            }
                    }
                    
                    } else {
                    
                    time[j][i] = 0;
                    dtime[0][j][i] = 0;
                    dtime[1][j][i] = 0;
                    dtime[2][j][i] = 0;
                    dtime[3][j][i] = 0;
                    dtime[4][j][i] = 0;
                    dtime[5][j][i] = 0;
                    dtime[6][j][i] = 0;
                    dtime[7][j][i] = 0;
                    dtime[8][j][i] = 0;
                    dtime[9][j][i] = 0;
                    dtime[10][j][i] = 0;
                    dtime[11][j][i] = 0;
                    dtime[12][j][i] = 0;
                    dtime[13][j][i] = 0;
                    dtime[14][j][i] = 0;
                    dtime[15][j][i] = 0;
                    
                    }
                    
                    count++;
    /*                sf_warning("%d of %d surface locations ;",count,nx*ny);*/
            }
        }
        
        sf_warning(" Time step: %d; of %d;",k+1,nt);
        sf_floatwrite(time[0],nx*ny,out);
        sf_floatwrite(dtime[0][0],nx*ny*npara,fit);
    
    }

  exit(0);
}












