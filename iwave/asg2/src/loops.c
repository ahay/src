#include "sgn.h"

void ploop(int n0,
	   ireal c1lam0, 
	   ireal c1lam1,
	   ireal c2lam0,
	   ireal c2lam1,
	   ireal tmp_ep1p,
	   ireal tmp_ep1pp,
	   ireal * restrict _p0,
	   ireal * restrict _p1,
	   ireal * restrict _mp,
	   ireal * restrict _v0,
	   ireal * restrict _v1p0,
	   ireal * restrict _v1p1,
	   ireal * restrict _v1m1,
	   ireal * restrict _v1m2,
	   ireal * restrict ep0p,
	   ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
  }
}

void v0loop(int n0,
	    ireal c1lam0,
	    ireal c2lam0,
	    ireal * restrict _p0,
	    ireal * restrict _mv0,
	    ireal * restrict _v0,
	    ireal * restrict ev0p,
	    ireal * restrict ev0pp) {

  int i0;

  for (i0=0; i0<n0; i0++) {   
    _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
      _mv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0])+ 
		c2lam0*(_p0[i0+2]-_p0[i0-1]));
  }
}

void v1loop(int n0,
	    ireal c1lam1,
	    ireal c2lam1,
	    ireal tmp_ev1p,
	    ireal tmp_ev1pp,
	    ireal * restrict _p1p0,
	    ireal * restrict _p1p1,
	    ireal * restrict _p1p2,
	    ireal * restrict _p1m1,
	    ireal * restrict _mv1,
	    ireal * restrict _v1) {
  int i0;

  for (i0=0; i0<n0; i0++) {
    _v1[i0] = tmp_ev1pp*_v1[i0] + tmp_ev1p*      
      _mv1[i0]*(c1lam1*(_p1p1[i0]-_p1p0[i0])+
		c2lam1*(_p1p2[i0]-_p1m1[i0]));
  }
}

