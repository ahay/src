#include "sgn.h"

void dploop(int n0,
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
	    ireal * restrict _rmp,
	    ireal * restrict _rv0,
	    ireal * restrict _rv1p0,
	    ireal * restrict _rv1p1,
	    ireal * restrict _rv1m1,
	    ireal * restrict _rv1m2,
	    ireal * restrict ep0p,
	    ireal * restrict ep0pp) {

  int i0;
  ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _rmp[i0]*
      (c1lam0*(_v0[i0]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0])) +
      _mp[i0]*
      (c1lam0*(_rv0[i0]-_rv0[i0-1])+c2lam0*(_rv0[i0+1]-_rv0[i0-2]) +
       c1lam1*(_rv1p0[i0]-_rv1m1[i0])+c2lam1*(_rv1p1[i0]-_rv1m2[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
  }
}

void dv0loop(int n0,
	     ireal c1lam0,
	     ireal c2lam0,
	     ireal * restrict _p0,
	     ireal * restrict _mv0,
	     ireal * restrict _rp0,
	     ireal * restrict _rmv0,
	     ireal * restrict _v0,
	     ireal * restrict ev0p,
	     ireal * restrict ev0pp) {

  int i0;
  for (i0=0; i0<n0; i0++) {
    _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
      (_rmv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0]) + 
		  c2lam0*(_p0[i0+2]-_p0[i0-1])) +
       _mv0[ i0]*(c1lam0*(_rp0[i0+1]-_rp0[i0]) + 
		  c2lam0*(_rp0[i0+2]-_rp0[i0-1])));
  }
}

void dv1loop(int n0,
	     ireal c1lam1,
	     ireal c2lam1,
	     ireal tmp_ev1p,
	     ireal tmp_ev1pp,
	     ireal * restrict _p1p0,
	     ireal * restrict _p1p1,
	     ireal * restrict _p1p2,
	     ireal * restrict _p1m1,
	     ireal * restrict _mv1,
	     ireal * restrict _rp1p0,
	     ireal * restrict _rp1p1,
	     ireal * restrict _rp1p2,
	     ireal * restrict _rp1m1,
	     ireal * restrict _rmv1,
	     ireal * restrict _v1) {
	     
  int i0;

  for (i0=0; i0<n0; i0++) {
    _v1[i0] = tmp_ev1pp*_v1[i0] + tmp_ev1p*      
      (_rmv1[i0]*(c1lam1*(_p1p1[i0]-_p1p0[i0]) +
		  c2lam1*(_p1p2[i0]-_p1m1[i0])) +
       _mv1[ i0]*(c1lam1*(_rp1p1[i0]-_rp1p0[i0])+
		  c2lam1*(_rp1p2[i0]-_rp1m1[i0])));
  }
}

void aploop(int n0,
	    ireal c1lam0,
	    ireal c1lam1,
	    ireal c2lam0,
	    ireal c2lam1,
	    ireal tmp_ep1p,
	    ireal tmp_ep1pp,
	    ireal * restrict _p0,
	    ireal * restrict _p1,
	    ireal * restrict _mp,
	    ireal * restrict _v0p0,
	    ireal * restrict _v0p1,
	    ireal * restrict _v0m1,
	    ireal * restrict _v0m2,
	    ireal * restrict _v1p0,
	    ireal * restrict _v1p1,
	    ireal * restrict _v1m1,
	    ireal * restrict _v1m2,
	    ireal * restrict _rmp,
	    ireal * restrict _rv0,
	    ireal * restrict _rv1p0,
	    ireal * restrict _rv1p1,
	    ireal * restrict _rv1m1,
	    ireal * restrict _rv1m2,
	    ireal * restrict ep0p,
	    ireal * restrict ep0pp) {

  int i0;
  ireal sdiv;
  ireal tmp __attribute__ ((__aligned__(16))); 
  ireal tmpq __attribute__ ((__aligned__(16))); 

  for (i0=0;i0<n0;i0++) {      
    
    tmp        = ep0pp[i0]*_rmp[i0]*_p0[i0]+tmp_ep1pp*_rmp[i0]*_p1[i0];//eta0post*_rmp[imp]*_p0[ip0];
    tmpq       = c1lam0*tmp; 
    _v0p0[i0] += tmpq;
    _v0m1[i0] -= tmpq;
    //    _v0[i0]   += tmpq;//eta0post*_rmp[imp]*c1lam0*_p0[ip0];
    //    _v0[i0-1] -= tmpq;//eta0post*_rmp[imp]*c1lam0*_p0[ip0];
    tmpq       = c2lam0*tmp;
    _v0p1[i0] += tmpq;
    _v0m2[i0] -= tmpq;
    //    _v0[i0+1] += tmpq;//eta0post*_rmp[imp]*c2lam0*_p0[ip0];
    //    _v0[i0-2] -= tmpq;//eta0post*_rmp[imp]*c2lam0*_p0[ip0];
    tmpq       = c1lam1*tmp;
    _v1p0[i0] += tmpq;//eta0post*_rmp[imp]*c1lam1*_p0[ip0];
    _v1m1[i0] -= tmpq;//eta0post*_rmp[imp]*c1lam1*_p0[ip0];
    tmpq       = c2lam1*tmp;
    _v1p1[i0] += tmpq;//eta0post*_rmp[imp]*c2lam1*_p0[ip0];
    _v1m2[i0] -= tmpq;//eta0post*_rmp[imp]*c2lam1*_p0[ip0];
  }

  for (i0=0;i0<n0;i0++) {

    sdiv = 
      (c1lam0*(_rv0[i0]-_rv0[i0-1])+c2lam0*(_rv0[i0+1]-_rv0[i0-2]) +
       c1lam1*(_rv1p0[i0]-_rv1m1[i0])+c2lam1*(_rv1p1[i0]-_rv1m2[i0]));

    _mp[i0]   += (ep0pp[i0]*_p0[i0]+tmp_ep1pp*_p1[i0])*sdiv;
    _p0[i0]    = ep0pp[i0]*ep0p[i0]*_p0[i0];//(eta0post * eta0pre - REAL_ONE) * _p0[ip0];
    _p1[i0]    = tmp_ep1pp*tmp_ep1p*_p1[i0];
    
  }
}

void av0loop(int n0,
	     ireal c1lam0,
	     ireal c2lam0,
	     ireal * restrict _p0m1,
	     ireal * restrict _p0p0,
	     ireal * restrict _p0p1,
	     ireal * restrict _p0p2,
	     ireal * restrict _mv0,
	     ireal * restrict _rp0,
	     ireal * restrict _rmv0,
	     ireal * restrict _ev0,
	     ireal * restrict _v0,
	     ireal * restrict ev0p,
	     ireal * restrict ev0pp) {

  int i0;
  ireal tmp;
  ireal tmpq;

  for (i0=0; i0<n0; i0++) {

    tmp        = ev0p[i0]*_rmv0[i0]*_v0[i0];
    tmpq       = c1lam0*tmp;
    _p0p1[i0] += tmpq;
    _p0p0[i0] -= tmpq;
    //    _p0[i0+1] += tmpq;
    //    _p0[i0  ] -= tmpq;
    tmpq       = c2lam0*tmp;
    _p0p2[i0] += tmpq;
    _p0m1[i0] -= tmpq;
    //    _p0[i0+2] += tmpq;
    //    _p0[i0-1] -= tmpq;
    
    _mv0[i0]  += ev0p[i0]*_v0[i0]*
      (c1lam0*(_rp0[i0+1]-_rp0[i0]) +
       c2lam0*(_rp0[i0+2]-_rp0[i0-1]));
    _v0[i0]    = ev0pp[i0]*_v0[i0]; 
    
  }
}

void av1loop(int n0,
	     ireal c1lam1,
	     ireal c2lam1,
	     ireal tmp_ev1p,
	     ireal tmp_ev1pp,
	     ireal * restrict _p1p0,
	     ireal * restrict _p1p1,
	     ireal * restrict _p1p2,
	     ireal * restrict _p1m1,
	     ireal * restrict _mv1,
	     ireal * restrict _rp1p0,
	     ireal * restrict _rp1p1,
	     ireal * restrict _rp1p2,
	     ireal * restrict _rp1m1,
	     ireal * restrict _rmv1,
	     ireal * restrict _v1) {

  int i0;
  ireal tmp;
  ireal tmpq;

  for (i0=0; i0<n0; i0++) {
    tmp        = tmp_ev1p*_rmv1[i0]*_v1[i0];
    tmpq       = c1lam1*tmp;
    _p1p1[i0] += tmpq;
    _p1p0[i0] -= tmpq;
    tmpq       = c2lam1*tmp;
    _p1p2[i0] += tmpq;
    _p1m1[i0] -= tmpq;
    
    _mv1[i0]   += tmp_ev1p*_v1[i0]*
      (c1lam1*(_rp1p1[i0]-_rp1p0[i0])+
       c2lam1*(_rp1p2[i0]-_rp1m1[i0]));
    _v1[i0]     = tmp_ev1pp*_v1[i0]; 
  }
}
