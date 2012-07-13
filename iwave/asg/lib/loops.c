#include "sgn.h"

/* ####################### P 2D ############################ */

void ploop1_2d(int n0,
	       ireal c1lam0, 
	       ireal c1lam1,
	       ireal tmp_ep1p,
	       ireal tmp_ep1pp,
	       ireal * restrict _p0,
	       ireal * restrict _p1,
	       ireal * restrict _mp,
	       ireal * restrict _v0,
	       ireal * restrict _v1p0,
	       ireal * restrict _v1m1,
	       ireal * restrict ep0p,
	       ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0]-_v0[i0-1]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
  }
}

void ploop2_2d(int n0,
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

void ploop4_2d(int n0,
	       ireal c1lam0, 
	       ireal c1lam1,
	       ireal c2lam0,
	       ireal c2lam1,
	       ireal c3lam0,
	       ireal c3lam1,
	       ireal c4lam0,
	       ireal c4lam1,
	       ireal tmp_ep1p,
	       ireal tmp_ep1pp,
	       ireal * restrict _p0,
	       ireal * restrict _p1,
	       ireal * restrict _mp,
	       ireal * restrict _v0,
	       ireal * restrict _v1p0,
	       ireal * restrict _v1p1,
	       ireal * restrict _v1p2,
	       ireal * restrict _v1p3,
	       ireal * restrict _v1m1,
	       ireal * restrict _v1m2,
	       ireal * restrict _v1m3,
	       ireal * restrict _v1m4,
	       ireal * restrict ep0p,
	       ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0  ]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
       c3lam0*(_v0[i0+2]-_v0[i0-3])+c4lam0*(_v0[i0+3]-_v0[i0-4]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]) +
       c3lam1*(_v1p2[i0]-_v1m3[i0])+c4lam1*(_v1p3[i0]-_v1m4[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
  }
}

/* ####################### P 3D ############################ */

void ploop1_3d(int n0,
	       ireal c1lam0, 
	       ireal c1lam1,
	       ireal c1lam2,
	       ireal tmp_ep1p,
	       ireal tmp_ep1pp,
	       ireal tmp_ep2p,
	       ireal tmp_ep2pp,
	       ireal * restrict _p0,
	       ireal * restrict _p1,
	       ireal * restrict _p2,
	       ireal * restrict _mp,
	       ireal * restrict _v0,
	       ireal * restrict _v1p0,
	       ireal * restrict _v1m1,
	       ireal * restrict _v2p0,
	       ireal * restrict _v2m1,
	       ireal * restrict ep0p,
	       ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0]-_v0[i0-1]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0]) +
       c1lam2*(_v2p0[i0]-_v2m1[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
    _p2[i0] = (_p2[i0]*tmp_ep2p + sdiv)*tmp_ep2pp;
  }
}

void ploop2_3d(int n0,
	       ireal c1lam0, 
	       ireal c1lam1,
	       ireal c1lam2,
	       ireal c2lam0,
	       ireal c2lam1,
	       ireal c2lam2,
	       ireal tmp_ep1p,
	       ireal tmp_ep1pp,
	       ireal tmp_ep2p,
	       ireal tmp_ep2pp,
	       ireal * restrict _p0,
	       ireal * restrict _p1,
	       ireal * restrict _p2,
	       ireal * restrict _mp,
	       ireal * restrict _v0,
	       ireal * restrict _v1p0,
	       ireal * restrict _v1p1,
	       ireal * restrict _v1m1,
	       ireal * restrict _v1m2,
	       ireal * restrict _v2p0,
	       ireal * restrict _v2p1,
	       ireal * restrict _v2m1,
	       ireal * restrict _v2m2,
	       ireal * restrict ep0p,
	       ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]) +
       c1lam2*(_v2p0[i0]-_v2m1[i0])+c2lam2*(_v2p1[i0]-_v2m2[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
    _p2[i0] = (_p2[i0]*tmp_ep2p + sdiv)*tmp_ep2pp;
  }
}

void ploop4_3d(int n0,
	       ireal c1lam0, 
	       ireal c1lam1,
	       ireal c1lam2,
	       ireal c2lam0,
	       ireal c2lam1,
	       ireal c2lam2,
	       ireal c3lam0,
	       ireal c3lam1,
	       ireal c3lam2,
	       ireal c4lam0,
	       ireal c4lam1,
	       ireal c4lam2,
	       ireal tmp_ep1p,
	       ireal tmp_ep1pp,
	       ireal tmp_ep2p,
	       ireal tmp_ep2pp,
	       ireal * restrict _p0,
	       ireal * restrict _p1,
	       ireal * restrict _p2,
	       ireal * restrict _mp,
	       ireal * restrict _v0,
	       ireal * restrict _v1p0,
	       ireal * restrict _v1p1,
	       ireal * restrict _v1p2,
	       ireal * restrict _v1p3,
	       ireal * restrict _v1m1,
	       ireal * restrict _v1m2,
	       ireal * restrict _v1m3,
	       ireal * restrict _v1m4,
	       ireal * restrict _v2p0,
	       ireal * restrict _v2p1,
	       ireal * restrict _v2p2,
	       ireal * restrict _v2p3,
	       ireal * restrict _v2m1,
	       ireal * restrict _v2m2,
	       ireal * restrict _v2m3,
	       ireal * restrict _v2m4,
	       ireal * restrict ep0p,
	       ireal * restrict ep0pp) {
  
  int i0;
  register ireal sdiv;

  for (i0=0;i0<n0;i0++) {
    sdiv = _mp[i0]*
      (c1lam0*(_v0[i0  ]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
       c3lam0*(_v0[i0+2]-_v0[i0-3])+c4lam0*(_v0[i0+3]-_v0[i0-4]) +
       c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]) +
       c3lam1*(_v1p2[i0]-_v1m3[i0])+c4lam1*(_v1p3[i0]-_v1m4[i0]) +
       c1lam2*(_v2p0[i0]-_v2m1[i0])+c2lam2*(_v2p1[i0]-_v2m2[i0]) +
       c3lam2*(_v2p2[i0]-_v2m3[i0])+c4lam2*(_v2p3[i0]-_v2m4[i0]));
    _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv)*ep0pp[i0];
    _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv)*tmp_ep1pp;
    _p2[i0] = (_p2[i0]*tmp_ep2p + sdiv)*tmp_ep2pp;
  }
}

/* ####################### V ############################ */

void v0loop1(int n0,
	     ireal c1lam0,
	     ireal * restrict _p0,
	     ireal * restrict _mv0,
	     ireal * restrict _v0,
	     ireal * restrict ev0p,
	     ireal * restrict ev0pp) {
  
  int i0;

  for (i0=0; i0<n0; i0++) {   
    _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
      _mv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0]));
  }
}

void v0loop2(int n0,
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

void v0loop4(int n0,
	     ireal c1lam0,
	     ireal c2lam0,
	     ireal c3lam0,
	     ireal c4lam0,
	     ireal * restrict _p0,
	     ireal * restrict _mv0,
	     ireal * restrict _v0,
	     ireal * restrict ev0p,
	     ireal * restrict ev0pp) {
  
  int i0;

  for (i0=0; i0<n0; i0++) {   
    _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
      _mv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0])+ 
		c2lam0*(_p0[i0+2]-_p0[i0-1])+ 
		c3lam0*(_p0[i0+3]-_p0[i0-2])+ 
		c4lam0*(_p0[i0+4]-_p0[i0-3]));
  }
}

// generic v loop - order 2
void vloop1(int n0,
	    ireal c1lam,
	    ireal tmp_evp,
	    ireal tmp_evpp,
	    ireal * restrict _pp0,
	    ireal * restrict _pp1,
	    ireal * restrict _mv,
	    ireal * restrict _v) {
  int i0;

  for (i0=0; i0<n0; i0++) {
    _v[i0] = tmp_evpp*_v[i0] + tmp_evp*      
      _mv[i0]*(c1lam*(_pp1[i0]-_pp0[i0]));
  }
}

// generic v loop - order 4
void vloop2(int n0,
	    ireal c1lam,
	    ireal c2lam,
	    ireal tmp_evp,
	    ireal tmp_evpp,
	    ireal * restrict _pp0,
	    ireal * restrict _pp1,
	    ireal * restrict _pp2,
	    ireal * restrict _pm1,
	    ireal * restrict _mv,
	    ireal * restrict _v) {
  int i0;

  for (i0=0; i0<n0; i0++) {
    _v[i0] = tmp_evpp*_v[i0] + tmp_evp*      
      _mv[i0]*(c1lam*(_pp1[i0]-_pp0[i0])+
	       c2lam*(_pp2[i0]-_pm1[i0]));
  }
}

// generic v loop - order 8
void vloop4(int n0,
	    ireal c1lam,
	    ireal c2lam,
	    ireal c3lam,
	    ireal c4lam,
	    ireal tmp_evp,
	    ireal tmp_evpp,
	    ireal * restrict _pp0,
	    ireal * restrict _pp1,
	    ireal * restrict _pp2,
	    ireal * restrict _pp3,
	    ireal * restrict _pp4,
	    ireal * restrict _pm1,
	    ireal * restrict _pm2,
	    ireal * restrict _pm3,
	    ireal * restrict _mv,
	    ireal * restrict _v) {
  int i0;

  for (i0=0; i0<n0; i0++) {
    _v[i0] = tmp_evpp*_v[i0] + tmp_evp*      
      _mv[i0]*(c1lam*(_pp1[i0]-_pp0[i0])+
	       c2lam*(_pp2[i0]-_pm1[i0])+
	       c3lam*(_pp3[i0]-_pm2[i0])+
	       c4lam*(_pp4[i0]-_pm3[i0]));
  }
}

