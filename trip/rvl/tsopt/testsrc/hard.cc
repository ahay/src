#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>

int main() {

  const int nt=9;
  float dt=0.01;
  
  float u=0.5;
  float du=1.0;

  int i;

  float stateHist[nt];

  stateHist[0] = 0.5;

  printf( "Starting Forward Evolution\n");

  printf("step = 0 u = 5.000000e-01 du = 1.000000e+00\n");
  for (i=0;i<nt;i++) {
    du = du*(1.0-2.0*dt*u);
    u = u + dt*(1.0-u*u);

    stateHist[i+1] = u;
    printf("step = %d u = %e du = %e\n",i+1,u,du);
  }

  printf( "Starting Backward Evolution\n");

  float au = 0.5;
  
  printf("step = 9 u[9] = %e au = 5.000000e-01\n", stateHist[9]); 
  for (i=nt-1; i>=0; --i) {
    
    au = au*(1.0-2.0*dt*stateHist[i]);
    printf("step = %d u[%d] = %e, au = %e\n", i, i, stateHist[i], au);
    
  }

}

