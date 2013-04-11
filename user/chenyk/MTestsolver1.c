/* Test for conjugate gradient, steepest descent, jacobi iteration, gauss-seidel iteration, successive over relaxation (SOR) iteration */
/* In this test, A is not symmetric, conjugate gradient doesn't converge, all the other four algorithm converge.*/

#include <stdio.h>
#include <rsf.h>
#include "cgsolver.h"
#include "sdsolver.h"
#include "classicsolver.h"

int main (void) {
  float** a;
  float x1[3], x2[3], x3[3], x4[3], x5[3], y[3];
  int i,j, iter;
  
  a = sf_floatalloc2(3,3);
  //symmetric positve definite
  a[0][0] = 8.; a[0][1] = -3.;a[0][2] = 2.; 
  a[1][0] = 4.; a[1][1] = 11.;a[1][2] = -1.; 
  a[2][0] = 6.; a[2][1] = 3. ;a[2][2] = 12.; 
  
  y[0]=20.; y[1]=33.; y[2]=36.;
  // reference x:3,2,1
  
  printf ("y = \n");
  for (i=0; i < 3; i ++) {
    printf (" %10.2f",y[i]);
  }
  printf ("\n");
  printf ("a = \n");
  for (j=0; j < 3; j ++) {
    for (i=0; i < 3; i ++) {
      printf (" %10.2f",a[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  
  printf ("conjugate gradient\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<3;i++){x1[i]=0;}
    cg_solve(a, x1, y, 3, iter);
    
    printf ("x = ");
    for (i=0; i < 3; i ++) {
      printf (" %12.8f",x1[i]);
    }
    printf ("\n");
    
  }
  
  printf ("steepest descent\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<3;i++){x2[i]=0;}
    sd_solve(a, x2, y, 3, iter);
    
    printf ("x = ");
    for (i=0; i < 3; i ++) {
      printf (" %12.8f",x2[i]);
    }
    printf ("\n");
    
  }
  
  printf ("jacobi iteration\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<3;i++){x3[i]=0;}
    
    jacobi_solve(a, x3, y, 3, iter);
    
    printf ("x = ");
    for (i=0; i < 3; i ++) {
      printf (" %12.8f",x3[i]);
    }
    printf ("\n");
    
  }
  
  printf ("gauss-seidel iteration\n------\n");
  for (iter =0; iter < 20; iter++) {
    for(i=0;i<3;i++){x4[i]=0;}
    gs_solve(a, x4, y, 3, iter);
   
    printf ("x = ");
	for (i=0; i < 3; i ++) {
	  printf (" %12.8f",x4[i]);
	}
	printf ("\n");
  }
  
  printf ("successive over relaxation\n------\n");
    for (iter =0; iter < 20; iter++) {
      for(i=0;i<3;i++){x5[i]=0;}
      sor_solve(a, x5, y, 3, iter, 0.5);
      
      printf ("x = ");
      for (i=0; i < 3; i ++) {
	printf (" %12.8f",x5[i]);
      }
      printf ("\n");
      
    }
    
    free(a[0]);
    free(a);   
    exit(0);
}
