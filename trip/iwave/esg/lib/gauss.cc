#include "gauss.h"

/*
float * getgauss(int * iw, float dt, float fpeak, float rref, float cref) {
  int i;
  float st;
  float pi=3.1415927;
  float * f = NULL;
  float fac=2*cref*cref*rref/(pi*fpeak*fpeak);
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      st = pi*fpeak*dt*i;
      f[i+*iw]= fac*exp(-st*st);
    }
  }
  return f;
}

float * getdgauss(int * iw, float dt, float fpeak, float rref, float cref) {
  int i;
  float st;
  float pi=3.1415927;
  float * f = NULL;
  float fac=2*cref*cref*rref/(fpeak);
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      st = pi*fpeak*dt*i;
      f[i+*iw]= -2 * st * fac*exp(-st*st);
    }
  }
  return f;
}

float * getrick(int * iw, float dt, float fpeak, float rref, float cref) {
  int i;
  float st;
  float pi=3.1415927;
  float fac=2*cref*cref*rref/(pi*fpeak*fpeak);
  float * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      st = pi*fpeak*dt*i;
      f[i+*iw]= -2*pi*pi*fpeak*fpeak*fac*exp(-st*st)*(1.0-2.0*st*st);
    }
  }
  return f;
}
*/

/* these functions return Gaussian, derivative of Gaussian,
   and Ricker pulse arrays truncated at approximately single
   precision round-off error (i.e. O(10^7)). The Gaussian is
   normalized so that the Ricker is its second derivative and 
   has unit amplitude at t=0.*/

ireal * igetgauss(int * iw, ireal dt, ireal fpeak) {
  int i;
  ireal * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (ireal *)malloc((2*(*iw)+1)*sizeof(ireal));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      f[i+*iw]= compgauss(i*dt,fpeak);
    }
  }
  return f;
}



float * getgauss(int * iw, float dt, float fpeak) {
  int i;
  float * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      f[i+*iw]= compgauss(i*dt,fpeak);
    }
  }
  return f;
}

ireal * igetdgauss(int * iw, ireal dt, ireal fpeak) {
  int i;
  ireal * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (ireal *)malloc((2*(*iw)+1)*sizeof(ireal));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      f[i+*iw]=compdgauss(i*dt,fpeak);
    }
  }
  return f;
}

float * getdgauss(int * iw, float dt, float fpeak) {
  int i;
  float * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      f[i+*iw]=compdgauss(i*dt,fpeak);
    }
  }
  return f;
}


float * getrick(int * iw, float dt, float fpeak) {
  int i;
  float * f = NULL;
  *iw = 1+(int)(1.2/(fpeak*dt)+0.1);
  f = (float *)malloc((2*(*iw)+1)*sizeof(float));
  if (f) {
    for (i=-*iw;i<*iw+1;i++) {
      f[i+*iw]=comprick(i*dt,fpeak);
    }
  }
  return f;
}

/* in the following functions, the Gaussian is normalized to produce a 
   Ricker of unit amplitude at t=0. */

float compgauss(float t, float fpeak) {
  float pi=3.1415927;
  float st=pi*pi*fpeak*fpeak;
  return -exp(-st*t*t)/(2*st);
}

float compdgauss(float t, float fpeak) {
  float pi=3.1415927;
  float st=pi*pi*fpeak*fpeak;
  return t*exp(-st*t*t);
}

float comprick(float t, float fpeak) {
  float pi=3.1415927;
  float st=pi*pi*fpeak*fpeak;
  float t2=t*t;
  return (1-2*st*t2)*exp(-st*t2);
}
