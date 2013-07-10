#include "sgn.h"

#define C1 ( -27.0e0/24.0e0 )
#define C2 ( 1.0e0/24.0e0 )

//#define LOOPFUN
#undef LOOPFUN

extern void dploop(int n0,
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
	    ireal * restrict ep0pp);

extern void dv0loop(int n0,
		    ireal c1lam0,
		    ireal c2lam0,
		    ireal * restrict _p0,
		    ireal * restrict _mv0,
		    ireal * restrict _rp0,
		    ireal * restrict _rmv0,
		    ireal * restrict _v0,
		    ireal * restrict ev0p,
		    ireal * restrict ev0pp);

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
	     ireal * restrict _v1);

extern void aploop(int n0,
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
		   ireal * restrict ep0pp);

extern void av0loop(int n0,
		    ireal c1lam0,
		    ireal c2lam0,
		    ireal * restrict _p0m1,
		    ireal * restrict _p0p0,
		    ireal * restrict _p0p1,
		    ireal * restrict _p0p2,
		    ireal * restrict _mv0,
		    ireal * restrict _rp0,
		    ireal * restrict _rmv0,
		    ireal * restrict _v0,
		    ireal * restrict ev0p,
		    ireal * restrict ev0pp);

extern void av1loop(int n0,
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
		    ireal * restrict _v1);

int duhasgfm24_2d_p(RDOM * dom, RDOM * rdom, void *pars) {

  // scaled Courant numbers
  ireal c1lam0, c1lam1, c2lam0, c2lam1;
  // counters
  int i0, i1;
  // inner loop length
  int n0;
  // scaled velocity divergence
  register ireal * restrict sdiv;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffmp, ioffv0, ioffv1;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1p0;
  register ireal * restrict _v1p1;
  register ireal * restrict _v1m1;
  register ireal * restrict _v1m2;
  register ireal * restrict _rmp;
  register ireal * restrict _rv0;
  register ireal * restrict _rv1p0;
  register ireal * restrict _rv1p1;
  register ireal * restrict _rv1m1;
  register ireal * restrict _rv1m2;

  register ireal * restrict ep0p;
  register ireal * restrict ep0pp;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ep1p __attribute__ ((__aligned__(16)));
  ireal tmp_ep1pp __attribute__ ((__aligned__(16)));

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // size of computational domain for P0 - same as for P1
  rd_gse(dom, D_P0, gsc, gec);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_EP0, offep0, ipntbuf);    
  rd_a_gse(dom, D_EP1, offep1, ipntbuf);    
  rd_a_gse(dom, D_P0, offp0, ipntbuf);    
  rd_a_gse(dom, D_P1, offp1, ipntbuf);    
  rd_a_gse(dom, D_MP0, offmp, ipntbuf);    
  rd_a_gse(dom, D_V0, offv0, ipntbuf);    
  rd_a_gse(dom, D_V1, offv1, ipntbuf);    

  // strides
  rd_a_size(dom, D_EP0, nep0);
  rd_a_size(dom, D_EP1, nep1);
  rd_a_size(dom, D_P0, np0);
  rd_a_size(dom, D_P1, np1);
  rd_a_size(dom, D_MP0, nmp);
  rd_a_size(dom, D_V0, nv0);
  rd_a_size(dom, D_V1, nv1);
  
  // field update loop
  n0   = gec[0]-gsc[0]+1;

  sdiv = (ireal *)usermalloc_(n0*sizeof(ireal));

  ep0p  = &((sgnpars->ep0_p)[gsc[0]-offep0[0]]);
  ep0pp = &((sgnpars->ep0_pp)[gsc[0]-offep0[0]]);

  for (i1=gsc[1];i1<gec[1]+1;i1++) {

    ioffp0 = -offp0[0] + (i1-offp0[1])*np0[0];
    ioffp1 = -offp1[0] + (i1-offp1[1])*np1[0];
    ioffmp = -offmp[0] + (i1-offmp[1])*nmp[0];
    ioffv0 = -offv0[0] + (i1-offv0[1])*nv0[0];
    ioffv1 = -offv1[0] + (i1-offv1[1])*nv1[0];

    tmp_ep1p  = (sgnpars->ep1_p)[i1-offep1[0]];
    tmp_ep1pp = (sgnpars->ep1_pp)[i1-offep1[0]];

    _p0   = &(((dom->_s)[D_P0 ]._s0)[gsc[0]+ioffp0]);
    _p1   = &(((dom->_s)[D_P1 ]._s0)[gsc[0]+ioffp1]);
    _mp   = &(((dom->_s)[D_MP0]._s0)[gsc[0]+ioffmp]);
    _v0   = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
    _v1p0 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
    _v1p1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
    _v1m1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
    _v1m2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);
    _rmp  = &(((rdom->_s)[D_MP0]._s0)[gsc[0]+ioffmp]);
    _rv0  = &(((rdom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
    _rv1p0= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
    _rv1p1= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
    _rv1m1= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
    _rv1m2= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);

#ifdef LOOPFUN
    dploop(n0,
	    c1lam0,
	    c1lam1,
	    c2lam0,
	    c2lam1,
	    tmp_ep1p,
	    tmp_ep1pp,
	    _p0,
	    _p1,
	    _mp,
	    _v0,
	    _v1p0,
	    _v1p1,
	    _v1m1,
	    _v1m2,
	    _rmp,
	    _rv0,
	    _rv1p0,
	    _rv1p1,
	    _rv1m1,
	    _rv1m2,
	    ep0p,
	    ep0pp);
#else 
    for (i0=0;i0<n0;i0++) {
      sdiv[i0] = _rmp[i0]*
	(c1lam0*(_v0[i0]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
	 c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]));			      
    }
    for (i0=0;i0<n0;i0++) {
      sdiv[i0] += _mp[i0]*
	(c1lam0*(_rv0[i0]-_rv0[i0-1])+c2lam0*(_rv0[i0+1]-_rv0[i0-2]) +
	 c1lam1*(_rv1p0[i0]-_rv1m1[i0])+c2lam1*(_rv1p1[i0]-_rv1m2[i0]));			      
    }
    for (i0=0;i0<n0;i0++) {
      _p0[i0] = (_p0[i0]*ep0p[i0] + sdiv[i0])*ep0pp[i0];
      _p1[i0] = (_p1[i0]*tmp_ep1p + sdiv[i0])*tmp_ep1pp;
    }
#endif
  }
  
  userfree_(sdiv);

  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;

  if (sgnpars->lbc[0]) {
    for (i1=0; i1<np0[1]; i1++) {
      _p0[i1*np0[0]]=-_p0[2+i1*np0[0]];
      _p0[1+i1*np0[0]]= REAL_ZERO;
    }
  }
  if (sgnpars->rbc[0]) {
    for (i1=0; i1<np0[1]; i1++) {
      _p0[np0[0]-1+i1*np0[0]]=-_p0[np0[0]-3+i1*np0[0]];
      _p0[np0[0]-2+i1*np0[0]]= REAL_ZERO;
    }
  }
  if (sgnpars->lbc[1]) {
    for (i0=0; i0<np1[0]; i0++) {
      _p1[i0]=-_p1[i0+2*np1[0]]; 
      _p1[i0+np1[0]]= REAL_ZERO;
    }
  }
  if (sgnpars->rbc[1]) {
    for (i0=0; i0<np1[0]; i0++) {
      _p1[i0+(np1[1]-1)*np1[0]]=-_p1[i0+(np1[1]-3)*np1[0]];
      _p1[i0+(np1[1]-2)*np1[0]]= REAL_ZERO;
    }
  }

  return 0;
}

int duhasgfm24_2d_v(RDOM * dom, RDOM * rdom, void * pars) {

  // scaled Courant numbers
  ireal c1lam0, c2lam0;
  ireal c1lam1, c2lam1;
  // counters
  int i0, i1;
  // inner loop limit
  int n0;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffmv0, ioffmv1, ioffv0, ioffv1;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  // workspace
  IPNT ipntbuf;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ev1p __attribute__ ((__aligned__(16)));
  ireal tmp_ev1pp __attribute__ ((__aligned__(16)));

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _mv0;
  register ireal * restrict _rp0;
  register ireal * restrict _rmv0;
  register ireal * restrict _v0;
  register ireal * restrict _p1p0;
  register ireal * restrict _p1p1;
  register ireal * restrict _p1p2;
  register ireal * restrict _p1m1;
  register ireal * restrict _mv1;
  register ireal * restrict _rp1p0;
  register ireal * restrict _rp1p1;
  register ireal * restrict _rp1p2;
  register ireal * restrict _rp1m1;
  register ireal * restrict _rmv1;
  register ireal * restrict _v1;

  register ireal * restrict ev0p;
  register ireal * restrict ev0pp;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // size of computational domain for V1 
  rd_gse(dom, D_V0, gsc0, gec0);    
  rd_gse(dom, D_V1, gsc1, gec1);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_P0,  offp0,  ipntbuf);    
  rd_a_gse(dom, D_MV0, offmv0, ipntbuf);    
  rd_a_gse(dom, D_EV0, offev0, ipntbuf);    
  rd_a_gse(dom, D_V0,  offv0,  ipntbuf);    
  rd_a_gse(dom, D_P1,  offp1,  ipntbuf);    
  rd_a_gse(dom, D_MV1, offmv1, ipntbuf);    
  rd_a_gse(dom, D_EV1, offev1, ipntbuf);    
  rd_a_gse(dom, D_V1,  offv1,  ipntbuf);    

  // strides
  rd_a_size(dom, D_P0,  np0);
  rd_a_size(dom, D_MV0, nmv0);
  rd_a_size(dom, D_EV0, nev0);
  rd_a_size(dom, D_V0,  nv0);
  rd_a_size(dom, D_P1,  np1);
  rd_a_size(dom, D_MV1, nmv1);
  rd_a_size(dom, D_EV1, nev1);
  rd_a_size(dom, D_V1,  nv1);

  ev0p  = &((sgnpars->ev0_p)[gsc0[0]-offev0[0]]);
  ev0pp = &((sgnpars->ev0_pp)[gsc0[0]-offev0[0]]);
  
  n0   = gec0[0]-gsc0[0]+1;

  for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
    ioffp0  = -offp0[0] + (i1 - offp0[1])*np0[0];
    ioffmv0 = -offmv0[0] + (i1 - offmv0[1])*nmv0[0];
    ioffv0  = -offv0[0] + (i1 - offv0[1])*nv0[0];

    _p0   = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0]);
    _mv0  = &(((dom->_s)[D_MV0]._s0)[gsc0[0]+ioffmv0]);
    _rp0  = &(((rdom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0]);
    _rmv0 = &(((rdom->_s)[D_MV0]._s0)[gsc0[0]+ioffmv0]);
    _v0   = &(((dom->_s)[D_V0 ]._s0)[gsc0[0]+ioffv0]);

#ifdef LOOPFUN
    dv0loop(n0,
	    c1lam0,
	    c2lam0,
	    _p0,
	    _mv0,
	    _rp0,
	    _rmv0,
	    _v0,
	    ev0p,
	    ev0pp);
#else
    for (i0=0; i0<n0; i0++) {
      _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
	(_rmv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0]) + 
		    c2lam0*(_p0[i0+2]-_p0[i0-1])) +
	 _mv0[ i0]*(c1lam0*(_rp0[i0+1]-_rp0[i0]) + 
		    c2lam0*(_rp0[i0+2]-_rp0[i0-1])));
    }
#endif
  }

  n0 = gec1[0]-gsc1[0]+1;

  for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {

    tmp_ev1p  = (sgnpars->ev1_p)[i1-offev1[0]];
    tmp_ev1pp = (sgnpars->ev1_pp)[i1-offev1[0]];

    ioffp1    = -offp1[0] + (i1 - offp1[1])*np1[0];
    ioffmv1   = -offmv1[0] + (i1 - offmv1[1])*nmv1[0];
    ioffv1    = -offv1[0] + (i1 - offv1[1])*nv1[0];

    _p1p0 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
    _p1p1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
    _p1p2 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
    _p1m1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-np1[0]]);
    _rp1p0= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
    _rp1p1= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
    _rp1p2= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
    _rp1m1= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-np1[0]]);

    _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
    _rmv1 = &(((rdom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);

    _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);

#ifdef LOOPFUN
    dv1loop(n0,
	    c1lam1,
	    c2lam1,
	    tmp_ev1p,
	    tmp_ev1pp,
	    _p1p0,
	    _p1p1,
	    _p1p2,
	    _p1m1,
	    _mv1,
	    _rp1p0,
	    _rp1p1,
	    _rp1p2,
	    _rp1m1,
	    _rmv1,
	    _v1);
#else
    for (i0=0; i0<n0; i0++) {
      _v1[i0] = tmp_ev1pp*_v1[i0] + tmp_ev1p*      
	(_rmv1[i0]*(c1lam1*(_p1p1[i0]-_p1p0[i0]) +
		    c2lam1*(_p1p2[i0]-_p1m1[i0])) +
	 _mv1[ i0]*(c1lam1*(_rp1p1[i0]-_rp1p0[i0])+
		    c2lam1*(_rp1p2[i0]-_rp1m1[i0])));
    }
#endif
  }

  _v0   = (dom->_s)[D_V0 ]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;

  if (sgnpars->lbc[0]) {
    for (i1=0; i1<nv0[1]; i1++) {    
      _v0[i1*nv0[0]]=_v0[1+i1*nv0[0]];
    }
  }
  if (sgnpars->rbc[0]) {
    for (i1=0; i1<nv0[1]; i1++) {    
      _v0[nv0[0]-1+i1*nv0[0]]=_v0[nv0[0]-2+i1*nv0[0]];
    }
  }

  if (sgnpars->lbc[1]) {
    for (i0=0; i0<nv1[0]; i0++) {    
      _v1[i0]=_v1[i0+nv1[0]];
    }
  }
  if (sgnpars->rbc[1]) {
    for (i0=0; i0<nv1[0]; i0++) {    
      _v1[i0+(nv1[1]-1)*nv1[0]]=_v1[i0+(nv1[1]-2)*nv1[0]];
    }
  }
  
  return 0;
}

int duhasgam24_2d_p(RDOM * dom, RDOM * rdom, void *pars) {

  // scaled Courant numbers
  ireal c1lam0, c1lam1, c2lam0, c2lam1;
  // counters
  int i0, i1;
  // inner loop length
  int n0;
  // workspace for inner loop hoisting
  ireal tmp, tmpq, sdiv;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffmp, ioffv0, ioffv1;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  //  register ireal * restrict _v0;
  register ireal * restrict _v0p0;
  register ireal * restrict _v0p1;
  register ireal * restrict _v0m1;
  register ireal * restrict _v0m2;
  register ireal * restrict _v1p0;
  register ireal * restrict _v1p1;
  register ireal * restrict _v1m1;
  register ireal * restrict _v1m2;
  register ireal * restrict _rmp;
  register ireal * restrict _rv0;
  register ireal * restrict _rv1p0;
  register ireal * restrict _rv1p1;
  register ireal * restrict _rv1m1;
  register ireal * restrict _rv1m2;

  register ireal * restrict ep0p;
  register ireal * restrict ep0pp;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ep1p __attribute__ ((__aligned__(16)));
  ireal tmp_ep1pp __attribute__ ((__aligned__(16)));

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // size of computational domain for P0 - same as for P1
  rd_gse(dom, D_P0, gsc, gec);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_EP0, offep0, ipntbuf);    
  rd_a_gse(dom, D_EP1, offep1, ipntbuf);    
  rd_a_gse(dom, D_P0, offp0, ipntbuf);    
  rd_a_gse(dom, D_P1, offp1, ipntbuf);    
  rd_a_gse(dom, D_MP0, offmp, ipntbuf);    
  rd_a_gse(dom, D_V0, offv0, ipntbuf);    
  rd_a_gse(dom, D_V1, offv1, ipntbuf);    

  // strides
  rd_a_size(dom, D_EP0, nep0);
  rd_a_size(dom, D_EP1, nep1);
  rd_a_size(dom, D_P0, np0);
  rd_a_size(dom, D_P1, np1);
  rd_a_size(dom, D_MP0, nmp);
  rd_a_size(dom, D_V0, nv0);
  rd_a_size(dom, D_V1, nv1);
  
  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;

  if (sgnpars->lbc[0]) {
    for (i1=0; i1<np0[1]; i1++) {
      //   _p0[i1*np0[0]]=-_p0[2+i1*np0[0]];
      //   _p0[1+i1*np0[0]]= REAL_ZERO;
      _p0[2+i1*np0[0]] -= _p0[i1*np0[0]];
      _p0[1+i1*np0[0]]  = REAL_ZERO;
      _p0[0+i1*np0[0]]  = REAL_ZERO;
    }
  }
  if (sgnpars->rbc[0]) {
    for (i1=0; i1<np0[1]; i1++) {
      //      _p0[np0[0]-1+i1*np0[0]]=-_p0[np0[0]-3+i1*np0[0]];
      //      _p0[np0[0]-2+i1*np0[0]]= REAL_ZERO;
      _p0[np0[0]-3+i1*np0[0]] -= _p0[np0[0]-1+i1*np0[0]];
      _p0[np0[0]-2+i1*np0[0]]  = REAL_ZERO;
      _p0[np0[0]-1+i1*np0[0]]  = REAL_ZERO;
    }
  }
  if (sgnpars->lbc[1]) {
    for (i0=0; i0<np1[0]; i0++) {
      //      _p1[i0]=-_p1[i0+2*np1[0]]; 
      //      _p1[i0+np1[0]]= REAL_ZERO;
      _p1[i0+2*np1[0]] -= _p1[i0];
      _p1[i0+1*np1[0]]  = REAL_ZERO;
      _p1[i0+0*np1[0]]  = REAL_ZERO;
    }
  }
  if (sgnpars->rbc[1]) {
    for (i0=0; i0<np1[0]; i0++) {
      //      _p1[i0+(np1[1]-1)*np1[0]]=-_p1[i0+(np1[1]-3)*np1[0]];
      //      _p1[i0+(np1[1]-2)*np1[0]]= REAL_ZERO;
      _p1[i0+(np1[1]-3)*np1[0]] -= _p1[i0+(np1[1]-1)*np1[0]];
      _p1[i0+(np1[1]-2)*np1[0]]  = REAL_ZERO;
      _p1[i0+(np1[1]-1)*np1[0]]  = REAL_ZERO;
    }
  }
  
  n0   = gec[0]-gsc[0]+1;

  ep0p  = &((sgnpars->ep0_p)[gsc[0]-offep0[0]]);
  ep0pp = &((sgnpars->ep0_pp)[gsc[0]-offep0[0]]);

  for (i1=gsc[1];i1<gec[1]+1;i1++) {

    ioffp0 = -offp0[0] + (i1-offp0[1])*np0[0];
    ioffp1 = -offp1[0] + (i1-offp1[1])*np1[0];
    ioffmp = -offmp[0] + (i1-offmp[1])*nmp[0];
    ioffv0 = -offv0[0] + (i1-offv0[1])*nv0[0];
    ioffv1 = -offv1[0] + (i1-offv1[1])*nv1[0];

    tmp_ep1p  = (sgnpars->ep1_p)[i1-offep1[0]];
    tmp_ep1pp = (sgnpars->ep1_pp)[i1-offep1[0]];

    _p0   = &(((dom->_s)[D_P0 ]._s0)[gsc[0]+ioffp0]);
    _p1   = &(((dom->_s)[D_P1 ]._s0)[gsc[0]+ioffp1]);
    _mp   = &(((dom->_s)[D_MP0]._s0)[gsc[0]+ioffmp]);
    //    _v0   = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
    _v0p0 = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
    _v0p1 = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0+1]);
    _v0m1 = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0-1]);
    _v0m2 = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0-2]);
    _v1p0 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
    _v1p1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
    _v1m1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
    _v1m2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);
    _rmp  = &(((rdom->_s)[D_MP0]._s0)[gsc[0]+ioffmp]);
    _rv0  = &(((rdom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
    _rv1p0= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
    _rv1p1= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
    _rv1m1= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
    _rv1m2= &(((rdom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);

#ifdef LOOPFUN
    aploop(n0,
	   c1lam0,
	   c1lam1,
	   c2lam0,
	   c2lam1,
	   tmp_ep1p,
	   tmp_ep1pp,
	   _p0,
	   _p1,
	   _mp,
	   _v0p0,
	   _v0p1,
	   _v0m1,
	   _v0m2,
	   _v1p0,
	   _v1p1,
	   _v1m1,
	   _v1m2,
	   _rmp,
	   _rv0,
	   _rv1p0,
	   _rv1p1,
	   _rv1m1,
	   _rv1m2,
	   ep0p,
	   ep0pp);
#else
    for (i0=0;i0<n0;i0++) {
      
      sdiv = 
	(c1lam0*(_rv0[i0]-_rv0[i0-1])+c2lam0*(_rv0[i0+1]-_rv0[i0-2]) +
	 c1lam1*(_rv1p0[i0]-_rv1m1[i0])+c2lam1*(_rv1p1[i0]-_rv1m2[i0]));

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

      _mp[i0]   += (ep0pp[i0]*_p0[i0]+tmp_ep1pp*_p1[i0])*sdiv;
      _p0[i0]    = ep0pp[i0]*ep0p[i0]*_p0[i0];//(eta0post * eta0pre - REAL_ONE) * _p0[ip0];
      _p1[i0]    = tmp_ep1pp*tmp_ep1p*_p1[i0];

    }
#endif
  }
  
  return 0;
}

int duhasgam24_2d_v(RDOM * dom, RDOM * rdom, void * pars) {

  // scaled Courant numbers
  ireal c1lam0, c2lam0;
  ireal c1lam1, c2lam1;
  // counters
  int i0, i1;
  // length of inner loop
  int n0;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffmv0, ioffmv1, ioffv0, ioffv1;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  // workspace
  IPNT ipntbuf;
  // workspace for inner loop hoists
  ireal tmp, tmpq;

  // field pointers - allocated arrays
  //  register ireal * restrict _p0;
  register ireal * restrict _p0p0;  
  register ireal * restrict _p0p1;  
  register ireal * restrict _p0p2;  
  register ireal * restrict _p0m1;
  register ireal * restrict _rp0;
  register ireal * restrict _mv0;
  register ireal * restrict _rmv0;
  register ireal * restrict _v0;
  register ireal * restrict _p1p0;
  register ireal * restrict _p1p1;
  register ireal * restrict _p1p2;
  register ireal * restrict _p1m1;
  register ireal * restrict _rp1p0;
  register ireal * restrict _rp1p1;
  register ireal * restrict _rp1p2;
  register ireal * restrict _rp1m1;
  register ireal * restrict _mv1;
  register ireal * restrict _rmv1;
  register ireal * restrict _v1;

  register ireal * restrict ev0p;
  register ireal * restrict ev0pp;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ev1p __attribute__ ((__aligned__(16)));
  ireal tmp_ev1pp __attribute__ ((__aligned__(16)));

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // size of computational domain for V1 
  rd_gse(dom, D_V0, gsc0, gec0);    
  rd_gse(dom, D_V1, gsc1, gec1);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_P0,  offp0,  ipntbuf);    
  rd_a_gse(dom, D_MV0, offmv0, ipntbuf);    
  rd_a_gse(dom, D_EV0, offev0, ipntbuf);    
  rd_a_gse(dom, D_V0,  offv0,  ipntbuf);    
  rd_a_gse(dom, D_P1,  offp1,  ipntbuf);    
  rd_a_gse(dom, D_MV1, offmv1, ipntbuf);    
  rd_a_gse(dom, D_EV1, offev1, ipntbuf);    
  rd_a_gse(dom, D_V1,  offv1,  ipntbuf);    

  // strides
  rd_a_size(dom, D_P0,  np0);
  rd_a_size(dom, D_MV0, nmv0);
  rd_a_size(dom, D_EV0, nev0);
  rd_a_size(dom, D_V0,  nv0);
  rd_a_size(dom, D_P1,  np1);
  rd_a_size(dom, D_MV1, nmv1);
  rd_a_size(dom, D_EV1, nev1);
  rd_a_size(dom, D_V1,  nv1);

  // bcs - dir 0
  _v0   = (dom->_s)[D_V0 ]._s0;

  if (sgnpars->lbc[0]) {
    for (i1=0; i1<nv0[1]; i1++) {    
      _v0[1+i1*nv0[0]] += _v0[i1*nv0[0]];
      _v0[i1*nv0[0]]    = REAL_ZERO;
    }
  }
  if (sgnpars->rbc[0]) {
    for (i1=0; i1<nv0[1]; i1++) {    
      _v0[nv0[0]-2+i1*nv0[0]] += _v0[nv0[0]-1+i1*nv0[0]];
      _v0[nv0[0]-1+i1*nv0[0]]  = REAL_ZERO;
    }
  }

  // body loop - dir 0
  ev0p  = &((sgnpars->ev0_p)[gsc0[0]-offev0[0]]);
  ev0pp = &((sgnpars->ev0_pp)[gsc0[0]-offev0[0]]);
  
  n0   = gec0[0]-gsc0[0]+1;

  for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
    ioffp0  = -offp0[0] + (i1 - offp0[1])*np0[0];
    ioffmv0 = -offmv0[0] + (i1 - offmv0[1])*nmv0[0];
    ioffv0  = -offv0[0] + (i1 - offv0[1])*nv0[0];

    //    _p0   = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0]);
    _p0m1 = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0-1]);
    _p0p0 = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0  ]);
    _p0p1 = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0+1]);
    _p0p2 = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0+2]);
    _rp0  = &(((rdom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0]);
    _mv0  = &(((dom->_s)[D_MV0]._s0)[gsc0[0]+ioffmv0]);
    _rmv0 = &(((rdom->_s)[D_MV0]._s0)[gsc0[0]+ioffmv0]);
    _v0   = &(((dom->_s)[D_V0 ]._s0)[gsc0[0]+ioffv0]);

#ifdef LOOPFUN 
    av0loop(n0,
	    c1lam0,
	    c2lam0,
	    _p0m1,
	    _p0p0,
	    _p0p1,
	    _p0p2,
	    _mv0,
	    _rp0,
	    _rmv0,
	    _v0,
	    ev0p,
	    ev0pp);
#else
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
#endif

  }

  // bcs - dir 1
  _v1   = (dom->_s)[D_V1 ]._s0;

  if (sgnpars->lbc[1]) {
    for (i0=0; i0<nv1[0]; i0++) {    
      _v1[i0+nv1[0]] += _v1[i0];
      _v1[i0]         = REAL_ZERO;
    }
  }
  if (sgnpars->rbc[1]) {
    for (i0=0; i0<nv1[0]; i0++) {
      _v1[i0+(nv1[1]-2)*nv1[0]] += _v1[i0+(nv1[1]-1)*nv1[0]];     
      _v1[i0+(nv1[1]-1)*nv1[0]]  = REAL_ZERO;
    }
  }

  // body loop - dir 1
  n0 = gec1[0]-gsc1[0]+1;

  for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {

    tmp_ev1p  = (sgnpars->ev1_p)[i1-offev1[0]];
    tmp_ev1pp = (sgnpars->ev1_pp)[i1-offev1[0]];

    ioffp1    = -offp1[0] + (i1 - offp1[1])*np1[0];
    ioffmv1   = -offmv1[0] + (i1 - offmv1[1])*nmv1[0];
    ioffv1    = -offv1[0] + (i1 - offv1[1])*nv1[0];

    _p1p0 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
    _p1p1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
    _p1p2 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
    _p1m1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-np1[0]]);
    _rp1p0= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
    _rp1p1= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
    _rp1p2= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
    _rp1m1= &(((rdom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-np1[0]]);
    _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
    _rmv1 = &(((rdom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
    _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);

#ifdef LOOPFUN
    av1loop(n0,
	    c1lam1,
	    c2lam1,
	    tmp_ev1p,
	    tmp_ev1pp,
	    _p1p0,
	    _p1p1,
	    _p1p2,
	    _p1m1,
	    _mv1,
	    _rp1p0,
	    _rp1p1,
	    _rp1p2,
	    _rp1m1,
	    _rmv1,
	    _v1);
#else
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
#endif
  }

  return 0;
}
    
int duhasgfm24_2d(RDOM *dom, RDOM * rdom, int iv, void *pars) {

  if ( iv == 0 ) return duhasgfm24_2d_p(dom, rdom, pars);
  if ( iv == 1 ) return duhasgfm24_2d_v(dom, rdom, pars);

  return 0;
}

int duhasgam24_2d(RDOM *dom, RDOM * rdom, int iv, void *pars) {

  if ( iv == 0 ) return duhasgam24_2d_p(dom, rdom, pars);
  if ( iv == 1 ) return duhasgam24_2d_v(dom, rdom, pars);

  return 0;
}



