#include "sgn.h"

/*
 * order 2
 */
extern void ploop1_2d(int n0,
		      ireal c1lam0, 
		      ireal c1lam1,
		      ireal tmp_ep1p,
		      ireal tmp_ep1pp,
		      ireal * restrict p0,
		      ireal * restrict p1,
		      ireal * restrict mp,
		      ireal * restrict v0,
		      ireal * restrict v1p0,
		      ireal * restrict v1m1,
		      ireal * restrict ep0p,
		      ireal * restrict ep0pp);

extern void ploop1_3d(int n0,
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
		      ireal * restrict ep0pp);

extern void v0loop1(int n0,
		    ireal c1lam0,
		    ireal * restrict _p0,
		    ireal * restrict _mv0,
		    ireal * restrict _v0,
		    ireal * restrict ev0p,
		    ireal * restrict ev0pp);

extern void vloop1(int n0,
		   ireal c1lam,
		   ireal tmp_evp,
		   ireal tmp_evpp,
		   ireal * restrict _pp0,
		   ireal * restrict _pp1,
		   ireal * restrict _mv,
		   ireal * restrict _v);

/*
 * order 4
 */
extern void ploop2_2d(int n0,
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

extern void ploop2_3d(int n0,
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
		      ireal * restrict ep0pp);

extern void v0loop2(int n0,
		    ireal c1lam0,
		    ireal c2lam0,
		    ireal * restrict _p0,
		    ireal * restrict _mv0,
		    ireal * restrict _v0,
		    ireal * restrict ev0p,
		    ireal * restrict ev0pp);

extern void vloop2(int n0,
		   ireal c1lam,
		   ireal c2lam,
		   ireal tmp_evp,
		   ireal tmp_evpp,
		   ireal * restrict _pp0,
		   ireal * restrict _pp1,
		   ireal * restrict _pp2,
		   ireal * restrict _pm1,
		   ireal * restrict _mv,
		   ireal * restrict _v);

/*
 * order 8
 */
extern void ploop4_2d(int n0,
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
		      ireal * restrict p0,
		      ireal * restrict p1,
		      ireal * restrict mp,
		      ireal * restrict v0,
		      ireal * restrict v1p0,
		      ireal * restrict v1p1,
		      ireal * restrict v1p2,
		      ireal * restrict v1p3,
		      ireal * restrict v1m1,
		      ireal * restrict v1m2,
		      ireal * restrict v1m3,
		      ireal * restrict v1m4,
		      ireal * restrict ep0p,
		      ireal * restrict ep0pp);

extern void ploop4_3d(int n0,
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
		      ireal * restrict ep0pp);

extern void v0loop4(int n0,
		    ireal c1lam0,
		    ireal c2lam0,
		    ireal c3lam0,
		    ireal c4lam0,
		    ireal * restrict _p0,
		    ireal * restrict _mv0,
		    ireal * restrict _v0,
		    ireal * restrict ev0p,
		    ireal * restrict ev0pp);

extern void vloop4(int n0,
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
		   ireal * restrict _v);


int asg_step_p(RDOM * dom, void *pars) {

  // counters
  int i0, i1, i2;
  // inner loop length
  int n0;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offep2, offp0, offp1, offp2, offmp, offv0, offv1, offv2;
  // inner loop offsets
  int ioffp0, ioffp1, ioffp2, ioffmp, ioffv0, ioffv1, ioffv2;
  // middle loop offsets 3D
  int ioffp0_3d, ioffp1_3d, ioffp2_3d, ioffmp_3d, ioffv0_3d, ioffv1_3d, ioffv2_3d;
  // strides = allocated rarray sizes, loop lengths
  IPNT nep0, nep1, nep2, np0, np1, np2, nmp, nv0, nv1, nv2;
  // half time step
  ireal dt2;

  // workspace
  IPNT ipntbuf;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ep1p __attribute__ ((__aligned__(16)));
  ireal tmp_ep1pp __attribute__ ((__aligned__(16)));
  ireal tmp_ep2p __attribute__ ((__aligned__(16)));
  ireal tmp_ep2pp __attribute__ ((__aligned__(16)));

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _p2;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1p0;
  register ireal * restrict _v1p1;
  register ireal * restrict _v1p2;
  register ireal * restrict _v1p3;
  register ireal * restrict _v1m1;
  register ireal * restrict _v1m2;
  register ireal * restrict _v1m3;
  register ireal * restrict _v1m4;
  register ireal * restrict _v2p0;
  register ireal * restrict _v2p1;
  register ireal * restrict _v2p2;
  register ireal * restrict _v2p3;
  register ireal * restrict _v2m1;
  register ireal * restrict _v2m2;
  register ireal * restrict _v2m3;
  register ireal * restrict _v2m4;

  register ireal * restrict ep0p;
  register ireal * restrict ep0pp;

  // fd parameter struct
  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // sanity-check order for current p loops
  if ((sgnpars->k !=1 ) && (sgnpars->k != 2) && (sgnpars->k != 4)) {
    // this is fatal so go ahead and blast a message
    if (retrieveRank()==0) {
      fprintf(stderr,"ERROR: asg_step_p\n");	
      fprintf(stderr,"half-order param = %d not one of permitted values = 1, 2, 4\n",sgnpars->k);
    }
    return E_BADINPUT;
  }

  // half timestep
  dt2 = sgnpars->dt / 2.0;

  // size of computational domain for P0 - same as for P1, P2
  rd_gse(dom, D_P0, gsc, gec);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_EP0, offep0, ipntbuf);    
  rd_a_gse(dom, D_P0, offp0, ipntbuf);    
  rd_a_gse(dom, D_V0, offv0, ipntbuf);    

  rd_a_size(dom, D_EP0, nep0);
  rd_a_size(dom, D_P0, np0);
  rd_a_size(dom, D_V0, nv0);

  rd_a_gse(dom, D_EP1, offep1, ipntbuf);    
  rd_a_gse(dom, D_P1, offp1, ipntbuf);    
  rd_a_gse(dom, D_V1, offv1, ipntbuf);    

  rd_a_size(dom, D_EP1, nep1);
  rd_a_size(dom, D_P1, np1);
  rd_a_size(dom, D_V1, nv1);

  rd_a_gse(dom, D_EP2, offep2, ipntbuf);    
  rd_a_gse(dom, D_P2, offp2, ipntbuf);    
  rd_a_gse(dom, D_V2, offv2, ipntbuf);    

  rd_a_size(dom, D_EP2, nep2);
  rd_a_size(dom, D_P2, np2);
  rd_a_size(dom, D_V2, nv2);

  rd_a_gse(dom, D_MP0, offmp, ipntbuf);    
  rd_a_size(dom, D_MP0, nmp);
  
  // field update loop
  /* version 1  */
  n0   = gec[0]-gsc[0]+1;

  ep0p  = &((sgnpars->ep0_p)[gsc[0]-offep0[0]]);
  ep0pp = &((sgnpars->ep0_pp)[gsc[0]-offep0[0]]);

  if (sgnpars->ndim==2) {

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

      if (sgnpars->k>0) {
	_v1p0 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
	_v1m1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
      }
      if (sgnpars->k>1) {
	_v1p1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
	_v1m2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);
      }
      if (sgnpars->k>3) {
	_v1p2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+2*nv1[0]]);
	_v1p3 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+3*nv1[0]]);
	_v1m3 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-3*nv1[0]]);
	_v1m4 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-4*nv1[0]]);
      }
    
      if (sgnpars->k==1) {
	ploop1_2d(n0,
		  (sgnpars->c11)[0],
		  (sgnpars->c11)[1],
		  tmp_ep1p,
		  tmp_ep1pp,
		  _p0,
		  _p1,
		  _mp,
		  _v0,
		  _v1p0,
		  _v1m1,
		  ep0p,
		  ep0pp);
      }
      if (sgnpars->k==2) {
	ploop2_2d(n0,
		  (sgnpars->c12)[0], //c1lam0, 
		  (sgnpars->c12)[1], //c1lam1,
		  (sgnpars->c22)[0], //c2lam0,
		  (sgnpars->c22)[1], //c2lam1,
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
      }
      if (sgnpars->k==4) {
	ploop4_2d(n0,
		  (sgnpars->c14)[0],
		  (sgnpars->c14)[1],
		  (sgnpars->c24)[0],
		  (sgnpars->c24)[1],
		  (sgnpars->c34)[0],
		  (sgnpars->c34)[1],
		  (sgnpars->c44)[0],
		  (sgnpars->c44)[1],
		  tmp_ep1p,
		  tmp_ep1pp,
		  _p0,
		  _p1,
		  _mp,
		  _v0,
		  _v1p0,
		  _v1p1,
		  _v1p2,
		  _v1p3,
		  _v1m1,
		  _v1m2,
		  _v1m3,
		  _v1m4,
		  ep0p,
		  ep0pp);
      }
    }
  }
  else if (sgnpars->ndim==3) {
    // fprintf(stderr,"gsc[0]=%d gec[0]=%d\n",gsc[0],gec[0]);
    // fprintf(stderr,"gsc[1]=%d gec[1]=%d\n",gsc[1],gec[1]);
    // fprintf(stderr,"gsc[2]=%d gec[2]=%d\n",gsc[2],gec[2]);
    for (i2=gsc[2];i2<gec[2]+1;i2++) {

      ioffp0_3d = -offp0[0] + (i2-offp0[2])*np0[0]*np0[1];
      ioffp1_3d = -offp1[0] + (i2-offp1[2])*np1[0]*np1[1];
      ioffp2_3d = -offp2[0] + (i2-offp2[2])*np2[0]*np2[1];
      ioffmp_3d = -offmp[0] + (i2-offmp[2])*nmp[0]*nmp[1];
      ioffv0_3d = -offv0[0] + (i2-offv0[2])*nv0[0]*nv0[1];
      ioffv1_3d = -offv1[0] + (i2-offv1[2])*nv1[0]*nv1[1];
      ioffv2_3d = -offv2[0] + (i2-offv2[2])*nv2[0]*nv2[1];

      tmp_ep2p  = (sgnpars->ep2_p)[i2-offep2[0]];
      tmp_ep2pp = (sgnpars->ep2_pp)[i2-offep2[0]];

      for (i1=gsc[1];i1<gec[1]+1;i1++) {
      
	ioffp0 = ioffp0_3d + (i1-offp0[1])*np0[0];
	ioffp1 = ioffp1_3d + (i1-offp1[1])*np1[0];
	ioffp2 = ioffp2_3d + (i1-offp2[1])*np2[0];
	ioffmp = ioffmp_3d + (i1-offmp[1])*nmp[0];
	ioffv0 = ioffv0_3d + (i1-offv0[1])*nv0[0];
	ioffv1 = ioffv1_3d + (i1-offv1[1])*nv1[0];
	ioffv2 = ioffv2_3d + (i1-offv2[1])*nv2[0];
	
	tmp_ep1p  = (sgnpars->ep1_p)[i1-offep1[0]];
	tmp_ep1pp = (sgnpars->ep1_pp)[i1-offep1[0]];
	
	_p0   = &(((dom->_s)[D_P0 ]._s0)[gsc[0]+ioffp0]);
	_p1   = &(((dom->_s)[D_P1 ]._s0)[gsc[0]+ioffp1]);
	_p2   = &(((dom->_s)[D_P2 ]._s0)[gsc[0]+ioffp2]);
	_mp   = &(((dom->_s)[D_MP0]._s0)[gsc[0]+ioffmp]);
	_v0   = &(((dom->_s)[D_V0 ]._s0)[gsc[0]+ioffv0]);
	
	if (sgnpars->k>0) {
	  _v1p0 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1]);
	  _v1m1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-nv1[0]]);
	  _v2p0 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2]);
	  _v2m1 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2-nv2[0]*nv2[1]]);
	}
	if (sgnpars->k>1) {
	  _v1p1 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+nv1[0]]);
	  _v1m2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-2*nv1[0]]);
	  _v2p1 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2+nv2[0]*nv2[1]]);
	  _v2m2 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2-2*nv2[0]*nv2[1]]);
	}
	if (sgnpars->k>2) {
	  _v1p2 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+2*nv1[0]]);
	  _v1m3 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-3*nv1[0]]);
	  _v2p2 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2+2*nv2[0]*nv2[1]]);
	  _v2m3 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2-3*nv2[0]*nv2[1]]);
	}
	if (sgnpars->k>3) {
	  _v1p3 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1+3*nv1[0]]);
	  _v1m4 = &(((dom->_s)[D_V1 ]._s0)[gsc[0]+ioffv1-4*nv1[0]]);
	  _v2p3 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2+3*nv2[0]*nv2[1]]);
	  _v2m4 = &(((dom->_s)[D_V2 ]._s0)[gsc[0]+ioffv2-4*nv2[0]*nv2[1]]);
	}
    
	// fprintf(stderr,"in asgstep: i1=%d i2=%d v0 offset=%d\n",i1,i2,gsc[0]+ioffv0);
	//	fprintf(stderr,"last entry in _p0 = %e\n",_p0[n0-1]);

	if (sgnpars->k==1) {
	  ploop1_3d(n0,
		    (sgnpars->c11)[0],
		    (sgnpars->c11)[1],
		    (sgnpars->c11)[2],
		    tmp_ep1p,
		    tmp_ep1pp,
		    tmp_ep2p,
		    tmp_ep2pp,
		    _p0,
		    _p1,
		    _p2,
		    _mp,
		    _v0,
		    _v1p0,
		    _v1m1,
		    _v2p0,
		    _v2m1,
		    ep0p,
		    ep0pp);
	}
	if (sgnpars->k==2) {
	  ploop2_3d(n0,
		    (sgnpars->c12)[0], //c1lam0, 
		    (sgnpars->c12)[1], //c1lam1,
		    (sgnpars->c12)[2], //c1lam1,
		    (sgnpars->c22)[0], //c2lam0,
		    (sgnpars->c22)[1], //c2lam1,
		    (sgnpars->c22)[2], //c2lam1,
		    tmp_ep1p,
		    tmp_ep1pp,
		    tmp_ep2p,
		    tmp_ep2pp,
		    _p0,
		    _p1,
		    _p2,
		    _mp,
		    _v0,
		    _v1p0,
		    _v1p1,
		    _v1m1,
		    _v1m2,
		    _v2p0,
		    _v2p1,
		    _v2m1,
		    _v2m2,
		    ep0p,
		    ep0pp);
	}
	if (sgnpars->k==4) {
	  ploop4_3d(n0,
		    (sgnpars->c14)[0],
		    (sgnpars->c14)[1],
		    (sgnpars->c14)[2],
		    (sgnpars->c24)[0],
		    (sgnpars->c24)[1],
		    (sgnpars->c24)[2],
		    (sgnpars->c34)[0],
		    (sgnpars->c34)[1],
		    (sgnpars->c34)[2],
		    (sgnpars->c44)[0],
		    (sgnpars->c44)[1],
		    (sgnpars->c44)[2],
		    tmp_ep1p,
		    tmp_ep1pp,
		    tmp_ep2p,
		    tmp_ep2pp,
		    _p0,
		    _p1,
		    _p2,
		    _mp,
		    _v0,
		    _v1p0,
		    _v1p1,
		    _v1p2,
		    _v1p3,
		    _v1m1,
		    _v1m2,
		    _v1m3,
		    _v1m4,
		    _v2p0,
		    _v2p1,
		    _v2p2,
		    _v2p3,
		    _v2m1,
		    _v2m2,
		    _v2m3,
		    _v2m4,
		    ep0p,
		    ep0pp);
	}
      }
    }
  }
  
  _p0   = (dom->_s)[D_P0 ]._s0;
  if (sgnpars->ndim > 1) _p1   = (dom->_s)[D_P1 ]._s0;
  if (sgnpars->ndim > 2) _p2   = (dom->_s)[D_P2 ]._s0;

  if (sgnpars->ndim==2) {
    if (sgnpars->k==2) {
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
    }
    if (sgnpars->k==4) {
      if (sgnpars->lbc[0]) {
	for (i1=0; i1<np0[1]; i1++) {
	  _p0[ 0+i1*np0[0]]=-_p0[6+i1*np0[0]];
	  _p0[ 1+i1*np0[0]]=-_p0[5+i1*np0[0]];
	  _p0[ 2+i1*np0[0]]=-_p0[4+i1*np0[0]];
	  _p0[ 3+i1*np0[0]]= REAL_ZERO;
	}
      }
      if (sgnpars->rbc[0]) {
	for (i1=0; i1<np0[1]; i1++) {
	  _p0[np0[0]-1+i1*np0[0]]=-_p0[np0[0]-7+i1*np0[0]];
	  _p0[np0[0]-2+i1*np0[0]]=-_p0[np0[0]-6+i1*np0[0]];
	  _p0[np0[0]-3+i1*np0[0]]=-_p0[np0[0]-5+i1*np0[0]];
	  _p0[np0[0]-4+i1*np0[0]]= REAL_ZERO;
	}
      }
      if (sgnpars->lbc[1]) {
	for (i0=0; i0<np1[0]; i0++) {
	  _p1[i0+0*np1[0]]=-_p1[i0+6*np1[0]]; 
	  _p1[i0+1*np1[0]]=-_p1[i0+5*np1[0]]; 
	  _p1[i0+2*np1[0]]=-_p1[i0+4*np1[0]]; 
	  _p1[i0+3*np1[0]]= REAL_ZERO;
	}
      }
      if (sgnpars->rbc[1]) {
	for (i0=0; i0<np1[0]; i0++) {
	  _p1[i0+(np1[1]-1)*np1[0]]=-_p1[i0+(np1[1]-7)*np1[0]];
	  _p1[i0+(np1[1]-2)*np1[0]]=-_p1[i0+(np1[1]-6)*np1[0]];
	  _p1[i0+(np1[1]-3)*np1[0]]=-_p1[i0+(np1[1]-5)*np1[0]];
	  _p1[i0+(np1[1]-4)*np1[0]]= REAL_ZERO;
	}
      }
    }
  }
  if (sgnpars->ndim==3) {
    if (sgnpars->k==2) {
      if (sgnpars->lbc[0]) {
	for (i2=0; i2<np0[2]; i2++) {
	  for (i1=0; i1<np0[1]; i1++) {
	    _p0[  (i1+i2*np0[1])*np0[0]]=-_p0[2+(i1+i2*np0[1])*np0[0]];
	    _p0[1+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->rbc[0]) {
	for (i2=0; i2<np0[2]; i2++) {
	  for (i1=0; i1<np0[1]; i1++) {
	    _p0[np0[0]-1+(i1+i2*np0[1])*np0[0]]=-_p0[np0[0]-3+(i1+np0[1])*np0[0]];
	    _p0[np0[0]-2+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->lbc[1]) {
	for (i2=0; i2<np1[2]; i2++) {
	  for (i0=0; i0<np1[0]; i0++) {
	    _p1[i0+i2*np1[0]*np1[1]]=-_p1[i0+i2*np1[0]*np1[1]+2*np1[0]]; 
	    _p1[i0+i2*np1[0]*np1[1]+np1[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->rbc[1]) {
	for (i2=0; i2<np1[2]; i2++) {
	  for (i0=0; i0<np1[0]; i0++) {
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-1)*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+(np1[1]-3)*np1[0]];
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-2)*np1[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->lbc[2]) {
	for (i1=0; i1<np2[1]; i1++) {
	  for (i0=0; i0<np2[0]; i0++) {
	    _p2[i0+i1*np2[0]+0*np2[0]*np2[1]]=-_p2[i0+i1*np2[0]+2*np2[0]*np2[1]];
	      // p2(i0,i1,1)=0
	    _p2[i0+i1*np2[0]+1*np2[0]*np2[1]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->rbc[2]) {
	for (i1=0; i1<np2[1]; i1++) {
	  for (i0=0; i0<np2[0]; i0++) {
	    _p2[i0+i1*np2[0]+(np2[2]-3)*np2[1]*np2[0]]=-_p2[i0+i1*np2[0]+(np2[2]-1)*np2[1]*np2[0]];
	    // p2(i0,i1,np2[2]-2)=0
	    _p2[i0+i1*np2[0]+(np2[2]-2)*np2[1]*np2[0]]= REAL_ZERO;
	  }
	}
      }
	
    }
    if (sgnpars->k==4) {
      if (sgnpars->lbc[0]) {
	for (i2=0; i2<np0[2]; i2++) {
	  for (i1=0; i1<np0[1]; i1++) {
	    _p0[ 0+(i1+i2*np0[1])*np0[0]]=-_p0[6+(i1+i2*np0[1])*np0[0]];
	    _p0[ 1+(i1+i2*np0[1])*np0[0]]=-_p0[5+(i1+i2*np0[1])*np0[0]];
	    _p0[ 2+(i1+i2*np0[1])*np0[0]]=-_p0[4+(i1+i2*np0[1])*np0[0]];
	    _p0[ 3+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->rbc[0]) {
	for (i2=0; i2<np0[2]; i2++) {
	  for (i1=0; i1<np0[1]; i1++) {
	    _p0[np0[0]-1+(i1+i2*np0[1])*np0[0]]=-_p0[np0[0]-7+(i1+i2*np0[1])*np0[0]];
	    _p0[np0[0]-2+(i1+i2*np0[1])*np0[0]]=-_p0[np0[0]-6+(i1+i2*np0[1])*np0[0]];
	    _p0[np0[0]-3+(i1+i2*np0[1])*np0[0]]=-_p0[np0[0]-5+(i1+i2*np0[1])*np0[0]];
	    _p0[np0[0]-4+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->lbc[1]) {
	for (i2=0; i2<np0[2]; i2++) {	
	  for (i0=0; i0<np1[0]; i0++) {
	    _p1[i0+i2*np1[0]*np1[1]+0*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+6*np1[0]]; 
	    _p1[i0+i2*np1[0]*np1[1]+1*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+5*np1[0]]; 
	    _p1[i0+i2*np1[0]*np1[1]+2*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+4*np1[0]]; 
	    _p1[i0+i2*np1[0]*np1[1]+3*np1[0]]= REAL_ZERO;
	    // p1(i0,3,i2) = 0
	  }
	}
      }
      if (sgnpars->rbc[1]) {
	for (i2=0; i2<np0[2]; i2++) {	
	  for (i0=0; i0<np1[0]; i0++) {
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-1)*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+(np1[1]-7)*np1[0]];
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-2)*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+(np1[1]-6)*np1[0]];
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-3)*np1[0]]=-_p1[i0+i2*np1[0]*np1[1]+(np1[1]-5)*np1[0]];
	    // p1(i0,np1[1]-4,i2)=0;
	    _p1[i0+i2*np1[0]*np1[1]+(np1[1]-4)*np1[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->lbc[2]) {
	for (i1=1; i1<np2[1]; i1++) {
	  for (i0=0; i0<np2[0]; i0++) {
	    _p2[i0+i1*np2[0]+0*np2[0]*np1[0]]=-_p2[i0+i1*np2[0]+6*np2[0]*np1[0]]; 
	    _p2[i0+i1*np2[0]+1*np2[0]*np1[0]]=-_p2[i0+i1*np2[0]+5*np2[0]*np1[0]]; 
	    _p2[i0+i1*np2[0]+2*np2[0]*np1[0]]=-_p2[i0+i1*np2[0]+4*np2[0]*np1[0]]; 
	    // p2(i0,i1,3)=0
	    _p2[i0+i1*np2[0]+3*np2[0]*np1[0]]= REAL_ZERO;
	  }
	}
      }
      if (sgnpars->rbc[2]) {
	for (i1=0; i1<np2[0]; i1++) {
	  for (i0=0; i0<np2[0]; i0++) {
	    _p2[i0+i1*np2[0]+(np2[2]-1)*np2[1]*np2[0]]=-_p2[i0+i1*np2[0]+(np2[2]-7)*np2[1]*np2[0]]; 
	    _p2[i0+i1*np2[0]+(np2[2]-2)*np2[1]*np2[0]]=-_p2[i0+i1*np2[0]+(np2[2]-6)*np2[1]*np2[0]]; 
	    _p2[i0+i1*np2[0]+(np2[2]-3)*np2[1]*np2[0]]=-_p2[i0+i1*np2[0]+(np2[2]-5)*np2[1]*np2[0]]; 
	    // p2(i0,i1,np2[2]-4) =0
	    _p2[i0+i1*np2[0]+(np2[2]-4)*np2[1]*np2[0]]= REAL_ZERO;           
	  }
	}
      }
    }
  }

  return 0;
}

int asg_step_v(RDOM * dom, void * pars) {

  // counters
  int i0, i1, i2;
  // length of inner loop
  int n0;
  // loop limits for computational array for velocity components
  IPNT gsc0, gec0; // v0
  IPNT gsc1, gec1; // v1
  IPNT gsc2, gec2; // v2
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  IPNT offp2, offev2, offv2, offmv2;
  // inner loop offsets for arrays
  int ioffp0, ioffp1, ioffp2;
  int ioffmv0, ioffmv1, ioffmv2;
  int ioffv0, ioffv1, ioffv2;
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  IPNT np2, nmv2, nev2, nv2;
  // workspace
  IPNT ipntbuf;

  // tmp scalar storage for outer loop PML coeffs
  ireal tmp_ev1p __attribute__ ((__aligned__(16)));
  ireal tmp_ev1pp __attribute__ ((__aligned__(16)));
  ireal tmp_ev2p __attribute__ ((__aligned__(16)));
  ireal tmp_ev2pp __attribute__ ((__aligned__(16)));

  // field pointers - allocated arrays
  register ireal * restrict _p0;
  register ireal * restrict _p1p0;
  register ireal * restrict _p1p1;
  register ireal * restrict _p1p2;
  register ireal * restrict _p1p3;
  register ireal * restrict _p1p4;
  register ireal * restrict _p1m1;
  register ireal * restrict _p1m2;
  register ireal * restrict _p1m3;
  register ireal * restrict _p2p0;
  register ireal * restrict _p2p1;
  register ireal * restrict _p2p2;
  register ireal * restrict _p2p3;
  register ireal * restrict _p2p4;
  register ireal * restrict _p2m1;
  register ireal * restrict _p2m2;
  register ireal * restrict _p2m3;
  register ireal * restrict _mv0;
  register ireal * restrict _mv1;
  register ireal * restrict _mv2;
  register ireal * restrict _v0;
  register ireal * restrict _v1;
  register ireal * restrict _v2;

  register ireal * restrict ev0p;
  register ireal * restrict ev0pp;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // sanity-check order for current v loops
  if ((sgnpars->k != 1) && (sgnpars->k != 2) && (sgnpars->k != 4)) {
    if (retrieveRank()==0) {
      fprintf(stderr,"ERROR: asg_step_v\n");	
      fprintf(stderr,"half-order param = %d not one of permitted values = 1, 2, 4\n",sgnpars->k);
    }
    return E_BADINPUT;
  }

  // zeroing out these arrays allows to write 1D or 2D loops
  // in 3D form, for at most slight penalty in index arith
  IASN(gsc0,IPNT_0);
  IASN(gec0,IPNT_0);
  IASN(gsc1,IPNT_0);
  IASN(gec1,IPNT_0);
  IASN(gsc2,IPNT_0);
  IASN(gec2,IPNT_0);
  IASN(np0,IPNT_1);
  IASN(nmv0,IPNT_1);
  IASN(nev0,IPNT_1);
  IASN(nv0,IPNT_1);
  IASN(np1,IPNT_1);
  IASN(nmv1,IPNT_1);
  IASN(nev1,IPNT_1);
  IASN(nv1,IPNT_1);
  IASN(np2,IPNT_1);
  IASN(nmv2,IPNT_1);
  IASN(nev2,IPNT_1);
  IASN(nv2,IPNT_1);
  IASN(offp0,IPNT_0);
  IASN(offmv0,IPNT_0);
  IASN(offev0,IPNT_0);
  IASN(offv0,IPNT_0);
  IASN(offp1,IPNT_0);
  IASN(offmv1,IPNT_0);
  IASN(offev1,IPNT_0);
  IASN(offv1,IPNT_0);
  IASN(offp2,IPNT_0);
  IASN(offmv2,IPNT_0);
  IASN(offev2,IPNT_0);
  IASN(offv2,IPNT_0);

  // size of computational domain for V components 
  if (sgnpars->ndim > 0) rd_gse(dom, D_V0, gsc0, gec0);    
  if (sgnpars->ndim > 1) rd_gse(dom, D_V1, gsc1, gec1);    
  if (sgnpars->ndim > 2) rd_gse(dom, D_V2, gsc2, gec2);    
  
  // start indices of allocated domains (index offsets)
  if (sgnpars->ndim > 0) {
    rd_a_gse(dom, D_P0,  offp0,  ipntbuf);    
    rd_a_gse(dom, D_MV0, offmv0, ipntbuf);    
    rd_a_gse(dom, D_EV0, offev0, ipntbuf);    
    rd_a_gse(dom, D_V0,  offv0,  ipntbuf);    
    rd_a_size(dom, D_P0,  np0);
    rd_a_size(dom, D_MV0, nmv0);
    rd_a_size(dom, D_EV0, nev0);
    rd_a_size(dom, D_V0,  nv0);
  }
  if (sgnpars->ndim > 1) {
    rd_a_gse(dom, D_P1,  offp1,  ipntbuf);    
    rd_a_gse(dom, D_MV1, offmv1, ipntbuf);    
    rd_a_gse(dom, D_EV1, offev1, ipntbuf);    
    rd_a_gse(dom, D_V1,  offv1,  ipntbuf);    
    rd_a_size(dom, D_P1,  np1);
    rd_a_size(dom, D_MV1, nmv1);
    rd_a_size(dom, D_EV1, nev1);
    rd_a_size(dom, D_V1,  nv1);
  }
  if (sgnpars->ndim > 2) {
    rd_a_gse(dom, D_P2,  offp2,  ipntbuf);    
    rd_a_gse(dom, D_MV2, offmv2, ipntbuf);    
    rd_a_gse(dom, D_EV2, offev2, ipntbuf);    
    rd_a_gse(dom, D_V2,  offv2,  ipntbuf);    
    rd_a_size(dom, D_P2,  np2);
    rd_a_size(dom, D_MV2, nmv2);
    rd_a_size(dom, D_EV2, nev2);
    rd_a_size(dom, D_V2,  nv2);
  }

  // strides

  ev0p  = &((sgnpars->ev0_p)[gsc0[0]-offev0[0]]);
  ev0pp = &((sgnpars->ev0_pp)[gsc0[0]-offev0[0]]);
  
  if (sgnpars->ndim > 0) {

    n0   = gec0[0]-gsc0[0]+1;

    for (i2=gsc0[2]; i2<gec0[2]+1; i2++) {
      for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
	
	ioffp0  = -offp0[0] + (i1 - offp0[1])*np0[0] + (i2-offp0[2])*np0[0]*np0[1];
	ioffmv0 = -offmv0[0] + (i1 - offmv0[1])*nmv0[0] + (i2-offmv0[2])*nmv0[0]*nmv0[1];
	ioffv0  = -offv0[0] + (i1 - offv0[1])*nv0[0] + (i2-offv0[2])*nv0[0]*nv0[1];
	
	_p0   = &(((dom->_s)[D_P0 ]._s0)[gsc0[0]+ioffp0]);
	_mv0  = &(((dom->_s)[D_MV0]._s0)[gsc0[0]+ioffmv0]);
	_v0   = &(((dom->_s)[D_V0 ]._s0)[gsc0[0]+ioffv0]);
	
	if (sgnpars->k == 1) {
	  v0loop1(n0,
		  (sgnpars->c11)[0],
		  _p0,
		  _mv0,
		  _v0,
		  ev0p,
		  ev0pp);
	}
	if (sgnpars->k == 2) {
	  v0loop2(n0,
		  (sgnpars->c12)[0],
		  (sgnpars->c22)[0],
		  _p0,
		  _mv0,
		  _v0,
		  ev0p,
		  ev0pp);
	}
	if (sgnpars->k == 4) {
	  v0loop4(n0,
		  (sgnpars->c14)[0],
		  (sgnpars->c24)[0],
		  (sgnpars->c34)[0],
		  (sgnpars->c44)[0],
		  _p0,
		  _mv0,
		  _v0,
		  ev0p,
		  ev0pp);
	}
      }
    }
  }

  if (sgnpars->ndim > 1) {

    n0 = gec1[0]-gsc1[0]+1;
    
    for (i2=gsc1[2]; i2<gec1[2]+1; i2++) {
      for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {

	tmp_ev1p  = (sgnpars->ev1_p)[i1-offev1[0]];
	tmp_ev1pp = (sgnpars->ev1_pp)[i1-offev1[0]];

	ioffp1    = -offp1[0] + (i1 - offp1[1])*np1[0] + (i2 - offp1[2])*np1[0]*np1[1];
	ioffmv1   = -offmv1[0] + (i1 - offmv1[1])*nmv1[0] + (i2 - offmv1[2])*nmv1[0]*nmv1[1];
	ioffv1    = -offv1[0] + (i1 - offv1[1])*nv1[0] + (i2 - offv1[2])*nv1[0]*nv1[1];

	if (sgnpars->k == 1) {

	  _p1p0 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
	  _p1p1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
	  _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
	  _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);
	  
	  vloop1(n0,
		 (sgnpars->c11)[1],
		 tmp_ev1p,
		 tmp_ev1pp,
		 _p1p0,
		 _p1p1,
		 _mv1,
		 _v1);
	}

	if (sgnpars->k == 2) {

	  _p1p0 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
	  _p1p1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
	  _p1p2 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
	  _p1m1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-np1[0]]);
	  _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
	  _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);
	  
	  vloop2(n0,
		 (sgnpars->c12)[1],
		 (sgnpars->c22)[1],
		 tmp_ev1p,
		 tmp_ev1pp,
		 _p1p0,
		 _p1p1,
		 _p1p2,
		 _p1m1,
		 _mv1,
		 _v1);
	}
	if (sgnpars->k == 4) {

	  _p1p0 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1]);
	  _p1p1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+np1[0]]);
	  _p1p2 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+2*np1[0]]);
	  _p1p3 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+3*np1[0]]);
	  _p1p4 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1+4*np1[0]]);
	  _p1m1 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-1*np1[0]]);
	  _p1m2 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-2*np1[0]]);
	  _p1m3 = &(((dom->_s)[D_P1 ]._s0)[gsc1[0]+ioffp1-3*np1[0]]);
	  _mv1  = &(((dom->_s)[D_MV1]._s0)[gsc1[0]+ioffmv1]);
	  _v1   = &(((dom->_s)[D_V1 ]._s0)[gsc1[0]+ioffv1]);
	  
	  vloop4(n0,
		 (sgnpars->c14)[1],
		 (sgnpars->c24)[1],
		 (sgnpars->c34)[1],
		 (sgnpars->c44)[1],
		 tmp_ev1p,
		 tmp_ev1pp,
		 _p1p0,
		 _p1p1,
		 _p1p2,
		 _p1p3,
		 _p1p4,
		 _p1m1,
		 _p1m2,
		 _p1m3,
		 _mv1,
		 _v1);
	}
      }
    }
  }

  if (sgnpars->ndim > 2) {

    n0 = gec2[0]-gsc2[0]+1;
    
    for (i2=gsc2[2]; i2<gec2[2]+1; i2++) {

      tmp_ev2p  = (sgnpars->ev2_p)[i2-offev2[0]];
      tmp_ev2pp = (sgnpars->ev2_pp)[i2-offev2[0]];

      for (i1=gsc2[1]; i1<gec2[1]+1; i1++) {

	ioffp2    = -offp2[0] + (i1 - offp2[1])*np2[0] + (i2 - offp2[2])*np2[0]*np2[1];
	ioffmv2   = -offmv2[0] + (i1 - offmv2[1])*nmv2[0] + (i2 - offmv2[2])*nmv2[0]*nmv2[1];
	ioffv2    = -offv2[0] + (i1 - offv2[1])*nv2[0] + (i2 - offv2[2])*nv2[0]*nv2[1];

	if (sgnpars->k == 1) {

	  _p2p0 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2]);
	  _p2p1 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+np2[0]*np2[1]]);
	  _mv2  = &(((dom->_s)[D_MV2]._s0)[gsc2[0]+ioffmv2]);
	  _v2   = &(((dom->_s)[D_V2 ]._s0)[gsc2[0]+ioffv2]);
	  
	  vloop1(n0,
		 (sgnpars->c11)[2],
		 tmp_ev2p,
		 tmp_ev2pp,
		 _p2p0,
		 _p2p1,
		 _mv2,
		 _v2);
	}

	if (sgnpars->k == 2) {

	  _p2p0 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2]);
	  _p2p1 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+1*np2[0]*np2[1]]);
	  _p2p2 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+2*np2[0]*np2[1]]);
	  _p2m1 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2-1*np2[0]*np2[1]]);
	  _mv2  = &(((dom->_s)[D_MV2]._s0)[gsc2[0]+ioffmv2]);
	  _v2   = &(((dom->_s)[D_V2 ]._s0)[gsc2[0]+ioffv2]);
	  
	  vloop2(n0,
		 (sgnpars->c12)[2],
		 (sgnpars->c22)[2],
		 tmp_ev2p,
		 tmp_ev2pp,
		 _p2p0,
		 _p2p1,
		 _p2p2,
		 _p2m1,
		 _mv2,
		 _v2);
	}
	if (sgnpars->k == 4) {

	  _p2p0 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2]);
	  _p2p1 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+1*np2[0]*np2[1]]);
	  _p2p2 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+2*np2[0]*np2[1]]);
	  _p2p3 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+3*np2[0]*np2[1]]);
	  _p2p4 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2+4*np2[0]*np2[1]]);
	  _p2m1 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2-1*np2[0]*np2[1]]);
	  _p2m2 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2-2*np2[0]*np2[1]]);
	  _p2m3 = &(((dom->_s)[D_P2 ]._s0)[gsc2[0]+ioffp2-3*np2[0]*np2[1]]);
	  _mv2  = &(((dom->_s)[D_MV2]._s0)[gsc2[0]+ioffmv2]);
	  _v2   = &(((dom->_s)[D_V2 ]._s0)[gsc2[0]+ioffv2]);
	  
	  vloop4(n0,
		 (sgnpars->c14)[2],
		 (sgnpars->c24)[2],
		 (sgnpars->c34)[2],
		 (sgnpars->c44)[2],
		 tmp_ev2p,
		 tmp_ev2pp,
		 _p2p0,
		 _p2p1,
		 _p2p2,
		 _p2p3,
		 _p2p4,
		 _p2m1,
		 _p2m2,
		 _p2m3,
		 _mv2,
		 _v2);
	}
      }
    }
  }

  if (sgnpars->ndim > 0) {

    _v0   = (dom->_s)[D_V0 ]._s0;

    if (sgnpars->k == 2) {

      if (sgnpars->lbc[0]) {
	for (i2=0; i2<nv0[2]; i2++) {
	  for (i1=0; i1<nv0[1]; i1++) {    
	    _v0[i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[1+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	  }
	}
      }
      if (sgnpars->rbc[0]) {
	for (i2=0; i2<nv0[2]; i2++) {
	  for (i1=0; i1<nv0[1]; i1++) {    
	    _v0[nv0[0]-1+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[nv0[0]-2+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	  }
	}
      }
    }

    if (sgnpars->k == 4) {
      if (sgnpars->lbc[0]) {
	for (i2=0; i2<nv0[2]; i2++) {
	  for (i1=0; i1<nv0[1]; i1++) {    
	    _v0[0+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[5+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	    _v0[1+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[4+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	    _v0[2+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[3+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	  }
	}
      }
      if (sgnpars->rbc[0]) {
	for (i2=0; i2<nv0[2]; i2++) {
	  for (i1=0; i1<nv0[1]; i1++) {    
	    _v0[nv0[0]-1+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[nv0[0]-6+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	    _v0[nv0[0]-2+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[nv0[0]-5+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	    _v0[nv0[0]-3+i1*nv0[0]+i2*nv0[0]*nv0[1]]=_v0[nv0[0]-4+i1*nv0[0]+i2*nv0[0]*nv0[1]];
	  }
	}
      }
    }      
  }

  if (sgnpars->ndim > 1) {

    _v1   = (dom->_s)[D_V1 ]._s0;

    if (sgnpars->k == 2) {
      if (sgnpars->lbc[1]) {
	for (i2=0; i2<nv1[2]; i2++) {
	  for (i0=0; i0<nv1[0]; i0++) {    
	    // v1[i1,0,i2]=vi[i1,1,i2]
	    _v1[i0+i2*nv1[0]*nv1[1]]=_v1[i0+nv1[0]+i2*nv1[0]*nv1[1]];
	  }
	}
      }
      
      if (sgnpars->rbc[1]) {
	for (i2=0; i2<nv1[2]; i2++) {
	  for (i0=0; i0<nv1[0]; i0++) {    
	    _v1[i0+(nv1[1]-1)*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+(nv1[1]-2)*nv1[0]+i2*nv1[0]*nv1[1]];
	  }
	}
      }
    }
    if (sgnpars->k == 4) {
      if (sgnpars->lbc[1]) {
	for (i2=0; i2<nv1[2]; i2++) {
	  for (i0=0; i0<nv1[0]; i0++) {    
	    _v1[i0+0*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+5*nv1[0]+i2*nv1[0]*nv1[1]];
	    _v1[i0+1*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+4*nv1[0]+i2*nv1[0]*nv1[1]];
	    // v1[i0,2,i2]=vi[i0,3,i2]
	    _v1[i0+2*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+3*nv1[0]+i2*nv1[0]*nv1[1]];
	  }
	}
      }
      
      if (sgnpars->rbc[1]) {
	for (i2=0; i2<nv1[2]; i2++) {
	  for (i0=0; i0<nv1[0]; i0++) {    
	    _v1[i0+(nv1[1]-1)*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+(nv1[1]-6)*nv1[0]+i2*nv1[0]*nv1[1]];
	    _v1[i0+(nv1[1]-2)*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+(nv1[1]-5)*nv1[0]+i2*nv1[0]*nv1[1]];
	    _v1[i0+(nv1[1]-3)*nv1[0]+i2*nv1[0]*nv1[1]]=_v1[i0+(nv1[1]-4)*nv1[0]+i2*nv1[0]*nv1[1]];
	  }
	}
      }
    }
  }

  if (sgnpars->ndim > 2) {

    _v2   = (dom->_s)[D_V2 ]._s0;

    if (sgnpars->k==2) {
      if (sgnpars->lbc[2]) {
	for (i1=0; i1<nv2[1]; i1++) {
	  for (i0=0; i0<nv2[0]; i0++) {    
	    // v2[i0,i1,0]=v2[i0,i1,1]
	    _v2[i0+i1*nv2[0]]=_v2[i0+i1*nv2[0]+nv2[0]*nv2[1]];
	  }
	}
      }
      
      if (sgnpars->rbc[2]) {
	for (i1=0; i1<nv2[1]; i1++) {
	  for (i0=0; i0<nv1[0]; i0++) {
	    // v2[i0,i1,nv2[2]-2) = v2[i0,i1,nv2[2]-1)    
	    _v2[i0+i1*nv2[0]+(nv2[2]-1)*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+(nv2[2]-2)*nv2[0]*nv2[1]];
	  }
	}
      }
    }

    if (sgnpars->k==4) {

      if (sgnpars->lbc[2]) {
	for (i1=0; i1<nv2[1]; i1++) {
	  for (i0=0; i0<nv2[0]; i0++) {    
	    _v2[i0+i1*nv2[0]+0*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+5*nv2[0]*nv2[1]];
	    _v2[i0+i1*nv2[0]+1*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+4*nv2[0]*nv2[1]];
	    _v2[i0+i1*nv2[0]+2*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+3*nv2[0]*nv2[1]];
	  }
	}
      }
      
      if (sgnpars->rbc[2]) {
	for (i1=0; i1<nv2[1]; i1++) {
	  for (i0=0; i0<nv1[0]; i0++) {
	    _v2[i0+i1*nv2[0]+(nv2[2]-1)*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+(nv2[2]-6)*nv2[0]*nv2[1]];
	    _v2[i0+i1*nv2[0]+(nv2[2]-2)*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+(nv2[2]-5)*nv2[0]*nv2[1]];
	    _v2[i0+i1*nv2[0]+(nv2[2]-3)*nv2[0]*nv2[1]]=_v2[i0+i1*nv2[0]+(nv2[2]-4)*nv2[0]*nv2[1]];
	  }
	}
      }
    }
    
  }

  return 0;
}
    
int asg_step(RDOM *dom, int iv, void *pars) {

  if ( iv == 0 ) return asg_step_p(dom, pars);
  if ( iv == 1 ) return asg_step_v(dom, pars);

  return 0;
}




