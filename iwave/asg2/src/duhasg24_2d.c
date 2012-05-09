#include "sgn.h"

#define C1 ( -27.0e0/24.0e0 )
#define C2 ( 1.0e0/24.0e0 )

#define LOOPFUN

#define PMLVEC
//#undef PMLVEC

extern void ploop(int n0,
		  ireal c1lam0, 
		  ireal c1lam1,
		  ireal c2lam0,
		  ireal c2lam1,
		  ireal tmp_ep1p,
		  ireal tmp_ep1pp,
		  ireal * restrict p0,
		  ireal * restrict p1,
		  ireal * restrict mp,
		  ireal * restrict v0,
		  ireal * restrict v1p0,
		  ireal * restrict v1p1,
		  ireal * restrict v1m1,
		  ireal * restrict v1m2,
		  ireal * restrict ep0p,
		  ireal * restrict ep0pp);

extern void v0loop(int,
		   ireal,
		   ireal,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict);

extern void v1loop(int,
		   ireal,
		   ireal,
		   ireal,
		   ireal,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict,
		   ireal * restrict);

int duhasg24_2d_p(RDOM * dom, void *pars) {

  //hi max
  // scaled Courant numbers
  ireal c1lam0 __attribute__ ((__aligned__(16)));
  ireal c1lam1 __attribute__ ((__aligned__(16)));
  ireal c2lam0 __attribute__ ((__aligned__(16)));
  ireal c2lam1 __attribute__ ((__aligned__(16)));
  // scaled velocity divergence
  register ireal * restrict sdiv;
  // counters
  int i0, i1;
  // inner loop length
  int n0;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffmp, ioffv0, ioffv1;
  // strides = allocated rarray sizes, loop lengths
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // half time step
  ireal dt2;

  // workspace
  IPNT ipntbuf;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ep1p __attribute__ ((__aligned__(16)));
  ireal tmp_ep1pp __attribute__ ((__aligned__(16)));

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1p0;
  register ireal * restrict _v1p1;
  register ireal * restrict _v1m1;
  register ireal * restrict _v1m2;

  register ireal * restrict ep0p;
  register ireal * restrict ep0pp;

  // fd parameter struct
  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // half timestep
  dt2 = sgnpars->dt / 2.0;

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
  /* version 1  */
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

#ifdef LOOPFUN
    ploop(n0,
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
	  ep0p,
	  ep0pp);
#else
    for (i0=0;i0<n0;i0++) {
      sdiv[i0] = _mp[i0]*
	(c1lam0*(_v0[i0]-_v0[i0-1])+c2lam0*(_v0[i0+1]-_v0[i0-2]) +
	 c1lam1*(_v1p0[i0]-_v1m1[i0])+c2lam1*(_v1p1[i0]-_v1m2[i0]));
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

int duhasg24_2d_v(RDOM * dom, void * pars) {

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

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ev1p __attribute__ ((__aligned__(16)));
  ireal tmp_ev1pp __attribute__ ((__aligned__(16)));

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _mv0;
  register ireal * restrict _v0;
  register ireal * restrict _p1p0;
  register ireal * restrict _p1p1;
  register ireal * restrict _p1p2;
  register ireal * restrict _p1m1;
  register ireal * restrict _mv1;
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
    _v0   = &(((dom->_s)[D_V0 ]._s0)[gsc0[0]+ioffv0]);
  
#ifdef LOOPFUN
    v0loop(n0,
	   c1lam0,
	   c2lam0,
	   _p0,
	   _mv0,
	   _v0,
	   ev0p,
	   ev0pp);
#else 
    for (i0=0; i0<n0; i0++) {

      _v0[i0] = ev0pp[i0]*_v0[i0] + ev0p[i0] * 
	_mv0[i0]*(c1lam0*(_p0[i0+1]-_p0[i0])+ 
		  c2lam0*(_p0[i0+2]-_p0[i0-1]));
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
    _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
    _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);
#ifdef LOOPFUN
    v1loop(n0,
	   c1lam1,
	   c2lam1,
	   tmp_ev1p,
	   tmp_ev1pp,
	   _p1p0,
	   _p1p1,
	   _p1p2,
	   _p1m1,
	   _mv1,
	   _v1);
#else
    for (i0=0; i0<n0; i0++) {
      _v1[i0] = tmp_ev1pp*_v1[i0] + tmp_ev1p*      
	_mv1[i0]*(c1lam1*(_p1p1[i0]-_p1p0[i0])+
		  c2lam1*(_p1p2[i0]-_p1m1[i0]));
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
    
int duhasg24_2d(RDOM *dom, int iv, void *pars) {

  if ( iv == 0 ) return duhasg24_2d_p(dom, pars);
  if ( iv == 1 ) return duhasg24_2d_v(dom, pars);

  return 0;
}




