#ifndef DYN1_SOLN_HH
#define DYN1_SOLN_HH

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
#include <vector>

float Dyn1_Fwd_Soln(const int tLevel) {
  
  const int nt=tLevel;
  float dt=0.01;  
  float u=0.5;
  float du=1.0;
  int i;
 
  for (i=0;i<nt;i++) {
    du = du*(1.0-2.0*dt*u);
    u = u + dt*(1.0-u*u);

   
  }
  
  
  return u;
}


float Dyn1_Lin_Soln(const int tLevel) {
  const int nt=tLevel;
  float dt=0.01;
  
  float u=0.5;
  float du=1.0;

  int i;

   
  for (i=0;i<nt;i++) {
    du = du*(1.0-2.0*dt*u);
    u = u + dt*(1.0-u*u);
  
  }

  return du;
}


float Dyn1_Adj_Soln(const int tLevel) {

  const int nt=9;
  float dt=0.01;
  
  float u=0.5;
 
  int i;

  float stateHist[nt];
 
  stateHist[0] = 0.5;

  
  for (i=0;i<nt-1;i++) {
    u = u + dt*(1.0-u*u);
    stateHist[i+1] = u;

  }

 
  float au = 0.5;
  //adjHist[nt] = 0.5; 
 
  for (i=nt-1; i>=tLevel; --i) {
    //cout << "i = " << i << endl;
    au = au*(1.0-2.0*dt*stateHist[i]);

  }
  
  return au;


}

#endif
