/* Test for conjugate gradient, steepest descent, jacobi iteration, gauss-seidel iteration, successive over relaxation (SOR) iteration */
/* In this test, A is symmetric,  because of inappropriate initial value, Jacobi iteration fails, all other four algorithm converge.*/

#include <stdio.h>
#include <rsf.h>
#include "cgsolver.h"
#include "sdsolver.h"
#include "classicsolver.h"

int main (void) {
  float** a;
  float x1[4], x2[4], x3[4], x4[4], x5[4], y[4];
  int i,j, iter;
  a = sf_floatalloc2(4,4);
  //symmetric positve definite
  a[0][0] = 4.; a[0][1] = 2.;	a[0][2] = 1.; a[0][3] = 0.;
  a[1][0] = 2.; a[1][1] = 4.;	a[1][2] = 2.; a[1][3] = 1.;
  a[2][0] = 1.; a[2][1] = 2.;	a[2][2] = 4.; a[2][3] = 2.;
  a[3][0] = 0.; a[3][1] = 1.;	a[3][2] = 2.; a[3][3] = 4.;
  
  y[0]=11.; y[1]=20.; y[2]=25.; y[3]=24.;
  // reference x: 1,2,3,4
  
  //not symmetric 
  //    a[0][0] = 4.; a[0][1] = 3.;	a[0][2] = 1.; a[0][3] = 0.;
  //    a[1][0] = 2.; a[1][1] = 4.;	a[1][2] = 3.; a[1][3] = 1.;
  //    a[2][0] = 1.; a[2][1] = 2.;	a[2][2] = 4.; a[2][3] = 3.;
  //    a[3][0] = 0.; a[3][1] = 1.;	a[3][2] = 2.; a[3][3] = 4.;
  //
  //    y[0]=13.; y[1]=23.; y[2]=29.; y[3]=24.;
  // reference x: 1,2,3,4
  
  printf ("y = \n");
  for (i=0; i < 4; i ++) {
    printf (" %10.2f",y[i]);
  }
  printf ("\n");
  printf ("a = \n");
  for (j=0; j < 4; j ++) {
    for (i=0; i < 4; i ++) {
      printf (" %10.2f",a[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  
  printf ("conjugate gradient\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<4;i++){x1[i]=0;}
    cg_solve(a, x1, y, 4, iter);
    
    printf ("x = ");
    for (i=0; i < 4; i ++) {
      printf (" %12.8f",x1[i]);
    }
    printf ("\n");  
  }
  
  printf ("steepest descent\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<4;i++){x2[i]=0;}
    sd_solve(a, x2, y, 4, iter);
    
    printf ("x = ");
    for (i=0; i < 4; i ++) {
      printf (" %12.8f",x2[i]);
    }
    printf ("\n");
    
  }
  
  printf ("jacobi iteration\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<4;i++){x3[i]=0;}
    
    jacobi_solve(a, x3, y, 4, iter);
    
    printf ("x = ");
    for (i=0; i < 4; i ++) {
      printf (" %12.8f",x3[i]);
    }
    printf ("\n");
    
  }
  
  printf ("gauss-seidel iteration\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<4;i++){x4[i]=0;}
    gs_solve(a, x4, y, 4, iter);
    
    printf ("x = ");
    for (i=0; i < 4; i ++) {
      printf (" %12.8f",x4[i]);
    }
    printf ("\n");
    
  }

  printf ("successive over relaxation\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<4;i++){x5[i]=0;}
    sor_solve(a, x5, y, 4, iter, 0.5);
    
    printf ("x = ");
    for (i=0; i < 4; i ++) {
      printf (" %12.8f",x5[i]);
    }
    printf ("\n");
    
  }
  
    free(a[0]);
    free(a);
    exit(0);
}
