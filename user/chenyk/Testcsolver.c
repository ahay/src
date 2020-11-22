/* Test for solving complex linear equation. 
Sometimes the solver is not stable.*/
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

#include <rsf.h>
#include "csolver.h"
/*
A=
3+1i 1+1i 1+1i
1+1i 3+2i 1+1i
1+1i 1+1i 3+3i
*/
/*
x=
1+1i
2+2i
3+3i
*/
/*
b=
2+14i
2+18i
0+24i
*/
/*The target is to solve Ax=b*/
int main(void)
{
    int j,niter=15,N=3;
    kiss_fft_cpx **a, *b, *x1, *x2, *x3, *x4, *x5;

    a=(kiss_fft_cpx**) sf_complexalloc2(N,N);
    b=(kiss_fft_cpx*) sf_complexalloc(N);
    x1=(kiss_fft_cpx*) sf_complexalloc(N);
    x2=(kiss_fft_cpx*) sf_complexalloc(N);
    x3=(kiss_fft_cpx*) sf_complexalloc(N);
    x4=(kiss_fft_cpx*) sf_complexalloc(N);
    x5=(kiss_fft_cpx*) sf_complexalloc(N);

    a[0][0].r=3; a[0][0].i=1;    a[0][1].r=1; a[0][1].i=1;    a[0][2].r=1; a[0][2].i=1;
    a[1][0].r=1; a[1][0].i=1;    a[1][1].r=3; a[1][1].i=2;    a[1][2].r=1; a[1][2].i=1;
    a[2][0].r=1; a[2][0].i=1;    a[2][1].r=1; a[2][1].i=1;    a[2][2].r=3; a[2][2].i=3;
    b[0].r=2; b[0].i=14;    b[1].r=2; b[1].i=18;    b[2].r=0; b[2].i=24;  

    for(j=0;j<3;j++)
      sf_warning("a[%d][0]=%.3f+i%.3f, a[%d][1]=%.3f+i%.3f, a[%d][2]=%.3f+i%.3f", j, a[j][0].r, a[j][0].i, j, a[j][1].r, a[j][1].i, j, a[j][2].r, a[j][2].i);

    sf_warning("b[0]=%.3f+i%.3f, b[1]=%.3f+i%.3f, b[2]=%.3f+i%.3f",b[0].r,b[0].i,b[1].r,b[1].i,b[2].r,b[2].i);

    csolve_init(a,b,N);
    
    gaussel_csolve(x1,niter);

    gs_csolve(x2,niter);

    jacobi_csolve(x3,niter);

    sor_csolve(x4,niter,0.5);
 
    sd_csolve(x5,niter);
     
    sf_warning("Gauss ellimination:");
    sf_warning("x[0]=%.3f+i%.3f, x[1]=%.3f+i%.3f, x[2]=%.3f+i%.3f",x1[0].r,x1[0].i,x1[1].r,x1[1].i,x1[2].r,x1[2].i);
    sf_warning("Gauss Seidel:");
    sf_warning("x[0]=%.3f+i%.3f, x[1]=%.3f+i%.3f, x[2]=%.3f+i%.3f",x2[0].r,x2[0].i,x2[1].r,x2[1].i,x2[2].r,x2[2].i);
    sf_warning("Jacobi:");
    sf_warning("x[0]=%.3f+i%.3f, x[1]=%.3f+i%.3f, x[2]=%.3f+i%.3f",x3[0].r,x3[0].i,x3[1].r,x3[1].i,x3[2].r,x3[2].i);
    sf_warning("Succesive over relaxation:");
    sf_warning("x[0]=%.3f+i%.3f, x[1]=%.3f+i%.3f, x[2]=%.3f+i%.3f",x4[0].r,x4[0].i,x4[1].r,x4[1].i,x4[2].r,x4[2].i);
    sf_warning("Steepest descent:");
    sf_warning("x[0]=%.3f+i%.3f, x[1]=%.3f+i%.3f, x[2]=%.3f+i%.3f",x5[0].r,x5[0].i,x5[1].r,x5[1].i,x5[2].r,x5[2].i);


     exit(0);
}

