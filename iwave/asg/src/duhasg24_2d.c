#include "sgn.h"

#define C1 ( -27.0e0/24.0e0 )
#define C2 ( 1.0e0/24.0e0 )

int duhasg24_2d_p(RDOM * dom, void *pars) {

  // scaled Courant numbers
  ireal c1lam0, c1lam1, c2lam0, c2lam1;
  // scaled velocity divergence
  //  ireal sdiv;
  // index variables for field update loop
  int iep0, iep1, ip0, ip1, iv0, iv1, imp, div0, div1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _ep0;
  register ireal * restrict _ep1;
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _ep0  = (dom->_s)[D_EP0]._s0;
  _ep1  = (dom->_s)[D_EP1]._s0;
  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _mp   = (dom->_s)[D_MP0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;

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
  for (i1=gsc[1];i1<gec[1]+1;i1++) {
    iep1= i1-offep1[0];
    eta1pre = (REAL_ONE - _ep1[iep1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ep1[iep1]*dt2);
    for (i0=gsc[0];i0<gec[0]+1;i0++) {
      iep0= i0-offep0[0];
      eta0pre = (REAL_ONE - _ep0[iep0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ep0[iep0]*dt2);
      ip0 = i0-offp0[0] + (i1-offp0[1])*np0[0];
      ip1 = i0-offp1[0] + (i1-offp1[1])*np1[0];
      imp = i0-offmp[0] + (i1-offmp[1])*nmp[0];
      iv0 = i0-offv0[0] + (i1-offv0[1])*nv0[0];
      iv1 = i0-offv1[0] + (i1-offv1[1])*nv1[0];
      div0 = 1;
      div1 = nv1[0];

      /*      
      sdiv = _mp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]));
      _p0[ip0] *= eta0pre;
      _p1[ip1] *= eta1pre;
      _p0[ip0] += sdiv;
      _p1[ip1] += sdiv;
      _p0[ip0] *= eta0post;
      _p1[ip1] *= eta1post;      
      */

      // inefficient but more transparent for tranformation

      _p0[ip0] += (eta0post * eta0pre - REAL_ONE) * _p0[ip0] + eta0post * _mp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]));

      _p1[ip1] += (eta1post * eta1pre - REAL_ONE) * _p1[ip1] + eta1post * _mp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]));
			      
    }
  }
  
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
  // index variables for field update loop
  int ip0, iev0, iv0, imv0, dip0;
  int ip1, iev1, iv1, imv1, dip1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _mv0;
  register ireal * restrict _ev0;
  register ireal * restrict _v0;
  register ireal * restrict _p1;
  register ireal * restrict _mv1;
  register ireal * restrict _ev1;
  register ireal * restrict _v1;


  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _p0   = (dom->_s)[D_P0 ]._s0;
  _mv0  = (dom->_s)[D_MV0]._s0;
  _ev0  = (dom->_s)[D_EV0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _mv1  = (dom->_s)[D_MV1]._s0;
  _ev1  = (dom->_s)[D_EV1]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // half timestep
  dt2 = sgnpars->dt / 2.0;

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

  for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
    for (i0=gsc0[0]; i0<gec0[0]+1; i0++) {
      iev0 = i0-offev0[0];
      eta0pre = (REAL_ONE - _ev0[iev0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ev0[iev0]*dt2);      
      ip0  = i0-offp0[0] + (i1 - offp0[1])*np0[0];
      dip0 = 1;
      imv0 = i0-offmv0[0] + (i1 - offmv0[1])*nmv0[0];
      iv0  = i0-offv0[0] + (i1 - offv0[1])*nv0[0];

      _v0[iv0] += (eta0pre*eta0post - REAL_ONE)*_v0[iv0] + eta0post * 
	_mv0[imv0]*(c1lam0*(_p0[ip0+dip0]-_p0[ip0])+
		    c2lam0*(_p0[ip0+2*dip0]-_p0[ip0-dip0]));
    }
  }

  for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {
    iev1 = i1-offev1[0];
    eta1pre = (REAL_ONE - _ev1[iev1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ev1[iev1]*dt2);      
    for (i0=gsc1[0]; i0<gec1[0]+1; i0++) {
      
      ip1  = i0-offp1[0] + (i1 - offp1[1])*np1[0];
      dip1 = np1[0];
      imv1 = i0-offmv1[0] + (i1 - offmv1[1])*nmv1[0];
      iv1  = i0-offv1[0] + (i1 - offv1[1])*nv1[0];

      _v1[iv1] += (eta1pre*eta1post-REAL_ONE)*_v1[iv1] + eta1post*
	_mv1[imv1]*(c1lam1*(_p1[ip1+dip1]-_p1[ip1])+
		    c2lam1*(_p1[ip1+2*dip1]-_p1[ip1-dip1]));
    }
  }

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
    
int duhasgfm24_2d_p(RDOM * dom, RDOM * rdom, void *pars) {

  // scaled Courant numbers
  ireal c1lam0, c1lam1, c2lam0, c2lam1;
  // scaled velocity divergence
  //  ireal sdiv;
  // index variables for field update loop
  int iep0, iep1, ip0, ip1, iv0, iv1, imp, div0, div1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _ep0;
  register ireal * restrict _ep1;
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1;
  register ireal * restrict _rmp;
  register ireal * restrict _rv0;
  register ireal * restrict _rv1;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _ep0  = (rdom->_s)[D_EP0]._s0;
  _ep1  = (rdom->_s)[D_EP1]._s0;
  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _mp   = (dom->_s)[D_MP0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;
  _rmp  = (rdom->_s)[D_MP0]._s0;
  _rv0  = (rdom->_s)[D_V0 ]._s0;
  _rv1  = (rdom->_s)[D_V1 ]._s0;

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
  for (i1=gsc[1];i1<gec[1]+1;i1++) {
    iep1= i1-offep1[0];
    eta1pre = (REAL_ONE - _ep1[iep1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ep1[iep1]*dt2);
    for (i0=gsc[0];i0<gec[0]+1;i0++) {
      iep0= i0-offep0[0];
      eta0pre = (REAL_ONE - _ep0[iep0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ep0[iep0]*dt2);
      ip0 = i0-offp0[0] + (i1-offp0[1])*np0[0];
      ip1 = i0-offp1[0] + (i1-offp1[1])*np1[0];
      imp = i0-offmp[0] + (i1-offmp[1])*nmp[0];
      iv0 = i0-offv0[0] + (i1-offv0[1])*nv0[0];
      iv1 = i0-offv1[0] + (i1-offv1[1])*nv1[0];
      div0 = 1;
      div1 = nv1[0];

      /*      
      sdiv = _mp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]));
      _p0[ip0] *= eta0pre;
      _p1[ip1] *= eta1pre;
      _p0[ip0] += sdiv;
      _p1[ip1] += sdiv;
      _p0[ip0] *= eta0post;
      _p1[ip1] *= eta1post;      
      */

      // inefficient but more transparent for tranformation

      _p0[ip0] += (eta0post * eta0pre - REAL_ONE) * _p0[ip0]  +
	eta0post * _rmp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1])) +
	eta0post * _mp[imp]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));

      _p1[ip1] += (eta1post * eta1pre - REAL_ONE) * _p1[ip1] + 
	eta1post * _rmp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1])) +
	eta1post * _mp[imp]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));
			      
    }
  }
  
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
  // index variables for field update loop
  int ip0, iev0, iv0, imv0, dip0;
  int ip1, iev1, iv1, imv1, dip1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _mv0;
  register ireal * restrict _rp0;
  register ireal * restrict _rmv0;
  register ireal * restrict _ev0;
  register ireal * restrict _v0;
  register ireal * restrict _p1;
  register ireal * restrict _mv1;
  register ireal * restrict _rp1;
  register ireal * restrict _rmv1;
  register ireal * restrict _ev1;
  register ireal * restrict _v1;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _p0   = (dom->_s)[D_P0 ]._s0;
  _rp0  = (rdom->_s)[D_P0 ]._s0;
  _mv0  = (dom->_s)[D_MV0]._s0;
  _rmv0 = (rdom->_s)[D_MV0]._s0;
  _ev0  = (rdom->_s)[D_EV0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _rp1  = (rdom->_s)[D_P1 ]._s0;
  _mv1  = (dom->_s)[D_MV1]._s0;
  _rmv1 = (rdom->_s)[D_MV1]._s0;
  _ev1  = (rdom->_s)[D_EV1]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // half timestep
  dt2 = sgnpars->dt / 2.0;

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

  for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
    for (i0=gsc0[0]; i0<gec0[0]+1; i0++) {
      iev0 = i0-offev0[0];
      eta0pre = (REAL_ONE - _ev0[iev0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ev0[iev0]*dt2);      
      ip0  = i0-offp0[0] + (i1 - offp0[1])*np0[0];
      dip0 = 1;
      imv0 = i0-offmv0[0] + (i1 - offmv0[1])*nmv0[0];
      iv0  = i0-offv0[0] + (i1 - offv0[1])*nv0[0];

      _v0[iv0] += (eta0pre*eta0post - REAL_ONE)*_v0[iv0] + eta0post * 
	_rmv0[imv0]*(c1lam0*(_p0[ip0+dip0]-_p0[ip0])+
		    c2lam0*(_p0[ip0+2*dip0]-_p0[ip0-dip0])) + eta0post * 
	_mv0[imv0]*(c1lam0*(_rp0[ip0+dip0]-_rp0[ip0])+
		    c2lam0*(_rp0[ip0+2*dip0]-_rp0[ip0-dip0]));
    }
  }

  for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {
    iev1 = i1-offev1[0];
    eta1pre = (REAL_ONE - _ev1[iev1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ev1[iev1]*dt2);      
    for (i0=gsc1[0]; i0<gec1[0]+1; i0++) {
      
      ip1  = i0-offp1[0] + (i1 - offp1[1])*np1[0];
      dip1 = np1[0];
      imv1 = i0-offmv1[0] + (i1 - offmv1[1])*nmv1[0];
      iv1  = i0-offv1[0] + (i1 - offv1[1])*nv1[0];

      _v1[iv1] += (eta1pre*eta1post-REAL_ONE)*_v1[iv1] + eta1post*
	_rmv1[imv1]*(c1lam1*(_p1[ip1+dip1]-_p1[ip1])+
		    c2lam1*(_p1[ip1+2*dip1]-_p1[ip1-dip1])) + eta1post*
	_mv1[imv1]*(c1lam1*(_rp1[ip1+dip1]-_rp1[ip1])+
		    c2lam1*(_rp1[ip1+2*dip1]-_rp1[ip1-dip1]));
    }
  }

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
  // scaled velocity divergence
  //  ireal sdiv;
  // index variables for field update loop
  int iep0, iep1, ip0, ip1, iv0, iv1, imp, div0, div1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offp0, offp1, offmp, offv0, offv1;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, np0, np1, nmp, nv0, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _ep0;
  register ireal * restrict _ep1;
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1;
  register ireal * restrict _rmp;
  register ireal * restrict _rv0;
  register ireal * restrict _rv1;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _ep0  = (rdom->_s)[D_EP0]._s0;
  _ep1  = (rdom->_s)[D_EP1]._s0;
  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _mp   = (dom->_s)[D_MP0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;
  _rmp  = (rdom->_s)[D_MP0]._s0;
  _rv0  = (rdom->_s)[D_V0 ]._s0;
  _rv1  = (rdom->_s)[D_V1 ]._s0;

  // half timestep
  dt2 = sgnpars->dt / 2.0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // size of computational domain for P0 - same as for P1
  rd_gse(rdom, D_P0, gsc, gec);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(rdom, D_EP0, offep0, ipntbuf);    
  rd_a_gse(rdom, D_EP1, offep1, ipntbuf);    
  rd_a_gse(rdom, D_P0, offp0, ipntbuf);    
  rd_a_gse(rdom, D_P1, offp1, ipntbuf);    
  rd_a_gse(rdom, D_MP0, offmp, ipntbuf);    
  rd_a_gse(rdom, D_V0, offv0, ipntbuf);    
  rd_a_gse(rdom, D_V1, offv1, ipntbuf);    

  // strides
  rd_a_size(rdom, D_EP0, nep0);
  rd_a_size(rdom, D_EP1, nep1);
  rd_a_size(rdom, D_P0, np0);
  rd_a_size(rdom, D_P1, np1);
  rd_a_size(rdom, D_MP0, nmp);
  rd_a_size(rdom, D_V0, nv0);
  rd_a_size(rdom, D_V1, nv1);
  
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
  
  // field update loop
  for (i1=gsc[1];i1<gec[1]+1;i1++) {
    iep1= i1-offep1[0];
    eta1pre = (REAL_ONE - _ep1[iep1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ep1[iep1]*dt2);
    for (i0=gsc[0];i0<gec[0]+1;i0++) {
      iep0= i0-offep0[0];
      eta0pre = (REAL_ONE - _ep0[iep0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ep0[iep0]*dt2);
      ip0 = i0-offp0[0] + (i1-offp0[1])*np0[0];
      ip1 = i0-offp1[0] + (i1-offp1[1])*np1[0];
      imp = i0-offmp[0] + (i1-offmp[1])*nmp[0];
      iv0 = i0-offv0[0] + (i1-offv0[1])*nv0[0];
      iv1 = i0-offv1[0] + (i1-offv1[1])*nv1[0];
      div0 = 1;
      div1 = nv1[0];

      /*      
      sdiv = _mp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]));
      _p0[ip0] *= eta0pre;
      _p1[ip1] *= eta1pre;
      _p0[ip0] += sdiv;
      _p1[ip1] += sdiv;
      _p0[ip0] *= eta0post;
      _p1[ip1] *= eta1post;      
      */

      // inefficient but more transparent for tranformation

      /*
      _p0[ip0] += (eta0post * eta0pre - REAL_ONE) * _p0[ip0]  +
	eta0post * _rmp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1])) +
	eta0post * _mp[imp]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));
      */

      _v0[iv0]          += eta0post*_rmp[imp]*c1lam0*_p0[ip0];
      _v0[iv0-div0]     -= eta0post*_rmp[imp]*c1lam0*_p0[ip0];
      _v0[iv0+div0]     += eta0post*_rmp[imp]*c2lam0*_p0[ip0];
      _v0[iv0-2*div0]   -= eta0post*_rmp[imp]*c2lam0*_p0[ip0];
      _v1[iv1]          += eta0post*_rmp[imp]*c1lam1*_p0[ip0];
      _v1[iv1-div1]     -= eta0post*_rmp[imp]*c1lam1*_p0[ip0];
      _v1[iv1+div1]     += eta0post*_rmp[imp]*c2lam1*_p0[ip0];
      _v1[iv1-2*div1]   -= eta0post*_rmp[imp]*c2lam1*_p0[ip0];

      _mp[imp]          += eta0post*_p0[ip0]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));
      _p0[ip0]          += (eta0post * eta0pre - REAL_ONE) * _p0[ip0];

      /*
      _p1[ip1] += (eta1post * eta1pre - REAL_ONE) * _p1[ip1] + 
	eta1post * _rmp[imp]*
	(c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	 c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1])) +
	eta1post * _mp[imp]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));
      */		      

      _v0[iv0]          += eta1post*_rmp[imp]*c1lam0*_p1[ip1];
      _v0[iv0-div0]     -= eta1post*_rmp[imp]*c1lam0*_p1[ip1];
      _v0[iv0+div0]     += eta1post*_rmp[imp]*c2lam0*_p1[ip1];
      _v0[iv0-2*div0]   -= eta1post*_rmp[imp]*c2lam0*_p1[ip1];
      _v1[iv1]          += eta1post*_rmp[imp]*c1lam1*_p1[ip1];
      _v1[iv1-div1]     -= eta1post*_rmp[imp]*c1lam1*_p1[ip1];
      _v1[iv1+div1]     += eta1post*_rmp[imp]*c2lam1*_p1[ip1];
      _v1[iv1-2*div1]   -= eta1post*_rmp[imp]*c2lam1*_p1[ip1];

      _mp[imp]          += eta1post*_p1[ip1]*
	(c1lam0*(_rv0[iv0]-_rv0[iv0-div0])+c2lam0*(_rv0[iv0+div0]-_rv0[iv0-2*div0]) +
	 c1lam1*(_rv1[iv1]-_rv1[iv1-div1])+c2lam1*(_rv1[iv1+div1]-_rv1[iv1-2*div1]));
      _p1[ip1]          += (eta1post * eta1pre - REAL_ONE) * _p1[ip1];


    }
  }
  
  return 0;
}

int duhasgam24_2d_v(RDOM * dom, RDOM * rdom, void * pars) {

  // scaled Courant numbers
  ireal c1lam0, c2lam0;
  ireal c1lam1, c2lam1;
  // index variables for field update loop
  int ip0, iev0, iv0, imv0, dip0;
  int ip1, iev1, iv1, imv1, dip1;
  // counters
  int i0, i1;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _mv0;
  register ireal * restrict _rp0;
  register ireal * restrict _rmv0;
  register ireal * restrict _ev0;
  register ireal * restrict _v0;
  register ireal * restrict _p1;
  register ireal * restrict _mv1;
  register ireal * restrict _rp1;
  register ireal * restrict _rmv1;
  register ireal * restrict _ev1;
  register ireal * restrict _v1;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _p0   = (dom->_s)[D_P0 ]._s0;
  _rp0  = (rdom->_s)[D_P0 ]._s0;
  _mv0  = (dom->_s)[D_MV0]._s0;
  _rmv0 = (rdom->_s)[D_MV0]._s0;
  _ev0  = (rdom->_s)[D_EV0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _rp1  = (rdom->_s)[D_P1 ]._s0;
  _mv1  = (dom->_s)[D_MV1]._s0;
  _rmv1 = (rdom->_s)[D_MV1]._s0;
  _ev1  = (rdom->_s)[D_EV1]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);

  // half timestep
  dt2 = sgnpars->dt / 2.0;

  // size of computational domain for V1 
  rd_gse(rdom, D_V0, gsc0, gec0);    
  rd_gse(rdom, D_V1, gsc1, gec1);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(rdom, D_P0,  offp0,  ipntbuf);    
  rd_a_gse(rdom, D_MV0, offmv0, ipntbuf);    
  rd_a_gse(rdom, D_EV0, offev0, ipntbuf);    
  rd_a_gse(rdom, D_V0,  offv0,  ipntbuf);    
  rd_a_gse(rdom, D_P1,  offp1,  ipntbuf);    
  rd_a_gse(rdom, D_MV1, offmv1, ipntbuf);    
  rd_a_gse(rdom, D_EV1, offev1, ipntbuf);    
  rd_a_gse(rdom, D_V1,  offv1,  ipntbuf);    

  // strides
  rd_a_size(rdom, D_P0,  np0);
  rd_a_size(rdom, D_MV0, nmv0);
  rd_a_size(rdom, D_EV0, nev0);
  rd_a_size(rdom, D_V0,  nv0);
  rd_a_size(rdom, D_P1,  np1);
  rd_a_size(rdom, D_MV1, nmv1);
  rd_a_size(rdom, D_EV1, nev1);
  rd_a_size(rdom, D_V1,  nv1);

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

  for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
    for (i0=gsc0[0]; i0<gec0[0]+1; i0++) {
      iev0 = i0-offev0[0];
      eta0pre = (REAL_ONE - _ev0[iev0]*dt2);
      eta0post = REAL_ONE/(REAL_ONE + _ev0[iev0]*dt2);      
      ip0  = i0-offp0[0] + (i1 - offp0[1])*np0[0];
      dip0 = 1;
      imv0 = i0-offmv0[0] + (i1 - offmv0[1])*nmv0[0];
      iv0  = i0-offv0[0] + (i1 - offv0[1])*nv0[0];

      /*
      _v0[iv0] += (eta0pre*eta0post - REAL_ONE)*_v0[iv0] + eta0post * 
      	_rmv0[imv0]*(c1lam0*(_p0[ip0+dip0]-_p0[ip0])+
		    c2lam0*(_p0[ip0+2*dip0]-_p0[ip0-dip0])) + eta0post * 
	_mv0[imv0]*(c1lam0*(_rp0[ip0+dip0]-_rp0[ip0])+
		    c2lam0*(_rp0[ip0+2*dip0]-_rp0[ip0-dip0]));
      */
      _p0[ip0+dip0]       += eta0post*_rmv0[imv0]*c1lam0*_v0[iv0];
      _p0[ip0]            -= eta0post*_rmv0[imv0]*c1lam0*_v0[iv0];
      _p0[ip0+2*dip0]     += eta0post*_rmv0[imv0]*c2lam0*_v0[iv0];
      _p0[ip0-dip0]       -= eta0post*_rmv0[imv0]*c2lam0*_v0[iv0];

      _mv0[imv0]          += eta0post*_v0[iv0]*
	(c1lam0*(_rp0[ip0+dip0]-_rp0[ip0])+
	 c2lam0*(_rp0[ip0+2*dip0]-_rp0[ip0-dip0]));
      _v0[iv0] += (eta0pre*eta0post - REAL_ONE)*_v0[iv0]; 
    }
  }

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

  for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {
    iev1 = i1-offev1[0];
    eta1pre = (REAL_ONE - _ev1[iev1]*dt2);
    eta1post = REAL_ONE/(REAL_ONE + _ev1[iev1]*dt2);      
  
    for (i0=gsc1[0]; i0<gec1[0]+1; i0++) {
      
      ip1  = i0-offp1[0] + (i1 - offp1[1])*np1[0];
      dip1 = np1[0];
      imv1 = i0-offmv1[0] + (i1 - offmv1[1])*nmv1[0];
      iv1  = i0-offv1[0] + (i1 - offv1[1])*nv1[0];

      /*
      _v1[iv1] += (eta1pre*eta1post-REAL_ONE)*_v1[iv1] + eta1post*
	_rmv1[imv1]*(c1lam1*(_p1[ip1+dip1]-_p1[ip1])+
		    c2lam1*(_p1[ip1+2*dip1]-_p1[ip1-dip1])) + eta1post*
	_mv1[imv1]*(c1lam1*(_rp1[ip1+dip1]-_rp1[ip1])+
		    c2lam1*(_rp1[ip1+2*dip1]-_rp1[ip1-dip1]));
      */

      _p1[ip1+dip1]       += eta1post*_rmv1[imv1]*c1lam1*_v1[iv1];
      _p1[ip1]            -= eta1post*_rmv1[imv1]*c1lam1*_v1[iv1];
      _p1[ip1+2*dip1]     += eta1post*_rmv1[imv1]*c2lam1*_v1[iv1];
      _p1[ip1-dip1]       -= eta1post*_rmv1[imv1]*c2lam1*_v1[iv1];

      _mv1[imv1]          += eta1post*_v1[iv1]*
	(c1lam1*(_rp1[ip1+dip1]-_rp1[ip1])+
	 c2lam1*(_rp1[ip1+2*dip1]-_rp1[ip1-dip1]));
      _v1[iv1] += (eta1pre*eta1post - REAL_ONE)*_v1[iv1]; 
    }
  }


  return 0;
}
    
int duhasg24_2d(RDOM *dom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return duhasg24_2d_p(dom, pars);
  if ( iarr == D_V0 ) return duhasg24_2d_v(dom, pars);

  return 0;
}


int duhasgfm24_2d(RDOM *dom, RDOM * rdom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return duhasgfm24_2d_p(dom, rdom, pars);
  if ( iarr == D_V0 ) return duhasgfm24_2d_v(dom, rdom, pars);

  return 0;
}

int duhasgam24_2d(RDOM *dom, RDOM * rdom, int iarr, void *pars) {

  if ( iarr == D_P0 ) return duhasgam24_2d_p(dom, rdom, pars);
  //  if ( iarr == D_P0 ) return duhasgam24_2d_v(dom, rdom, pars);
  //  if ( iarr == D_P0 ) return asg_atsm2d_24p01(dom, rdom, pars);
  if ( iarr == D_V0 ) return duhasgam24_2d_v(dom, rdom, pars);
  //  if ( iarr == D_V0 ) return duhasgam24_2d_p(dom, rdom, pars);
  //  if ( iarr == D_V0 ) return (asg_atsm2d_24v0(dom, rdom, pars) || asg_atsm2d_24v1(dom, rdom, pars));

  return 0;
}


