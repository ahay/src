#ifndef __RVL_OPTTESTS_
#define __RVL_OPTTESTS_

#include "utility.hh"

using namespace RVL;

// Rosenbrock function
template<typename T>
int rosie(int m, int n, T const * x, T * y) {
  if ((m!=2) || (n!=2)) return 1;
  y[0]=10*(x[1]-x[0]*x[0]);
  y[1]=ScalarFieldTraits<T>::One()-x[0];
  return 0;
}

// Rosenbrock jacobian
template<typename T>
int rosiejac(int m, int n, T const * x, T * y) {
  if ((m!=2) || (n!=2)) return 1;
  // df[0]/dx[0]
  y[0]=-20*x[0];
  // df[1]/dx[0]
  y[1]=-ScalarFieldTraits<T>::One();
  // df[0]/dx[1]
  y[2]=10*ScalarFieldTraits<T>::One();
  // df[1]/dx[1]
  y[3]=ScalarFieldTraits<T>::Zero();
  return 0;
}

// diagonal quartic function
template<typename T>
int diagquart(int m, int n, T const * x, T * y) {
  if ((m != n) || (n < 1)) return 1; 
  for (int i=0;i<n;i++) {
    T z = (1.0+(T)i)*x[i];
    y[i] = z - 0.1*z*z+0.05*z*z*z*z;
  }
  return 0;
}

// diagonal quartic jacobian
template<typename T> 
int diagquartjac(int m, int n, T const * x, T * y) {
  if ((m != n) || (n < 1)) return 1; 
  for (int i=0;i<n*n;i++) y[i]=ScalarFieldTraits<T>::Zero();
  for (int i=0;i<n;i++) {
    T z =  (1.0+(T)i)*x[i];
    y[i+i*n]=(1.0-0.2*z+0.2*z*z*z);
  }
  return 0;
}

#endif
