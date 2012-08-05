#include "sgn.h"

#define C1 ( -27.0e0/24.0e0 )
#define C2 ( 1.0e0/24.0e0 )

int asg24_3d_p(RDOM * dom, void *pars) {

  //hi max
  // scaled Courant numbers
  ireal c1lam0, c1lam1, c1lam2, c2lam0, c2lam1, c2lam2;
  // scaled velocity divergence
  ireal sdiv;
  // index variables for field update loop
  int iep0, iep1, iep2, ip0, ip1, ip2, iv0, iv1, iv2, imp, div0, div1, div2;
  // counters
  int i0, i1, i2;
  // loop limits for computational array for pressure
  IPNT gsc, gec;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offep0, offep1, offep2;
  IPNT offp0, offp1, offp2, offmp;
  IPNT offv0, offv1, offv2;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffp1, ioffp2, ioffmp;
  int ioffv0, ioffv1, ioffv2;
  int ioffp0_2, ioffp1_2, ioffp2_2, ioffmp_2;
  int ioffv0_2, ioffv1_2, ioffv2_2;
  // strides = allocated rarray sizes
  IPNT nep0, nep1, nep2;
  IPNT np0, np1, np2, nmp;
  IPNT nv0, nv1, nv2;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  ireal eta2pre, eta2post;
  // workspace
  IPNT ipntbuf;

  // field pointers - allocated arrays
  register ireal * restrict _ep0;
  register ireal * restrict _ep1;
  register ireal * restrict _ep2;
  register ireal * restrict _p0;
  register ireal * restrict _p1;
  register ireal * restrict _p2;
  register ireal * restrict _mp;
  register ireal * restrict _v0;
  register ireal * restrict _v1;
  register ireal * restrict _v2;

  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  

  // assign field pointers 
  _ep0  = (dom->_s)[D_EP0]._s0;
  _ep1  = (dom->_s)[D_EP1]._s0;
  _ep2  = (dom->_s)[D_EP2]._s0;
  _p0   = (dom->_s)[D_P0 ]._s0;
  _p1   = (dom->_s)[D_P1 ]._s0;
  _p2   = (dom->_s)[D_P2 ]._s0;
  _mp   = (dom->_s)[D_MP0]._s0;
  _v0   = (dom->_s)[D_V0 ]._s0;
  _v1   = (dom->_s)[D_V1 ]._s0;
  _v2   = (dom->_s)[D_V2 ]._s0;

  // half timestep
  dt2 = sgnpars->dt / 2.0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c1lam2 = C1*(sgnpars->lam[2]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c2lam1 = C2*(sgnpars->lam[1]);
  c2lam2 = C2*(sgnpars->lam[2]);

  // size of computational domain for P0 - same as for P1, P2
  rd_gse(dom, D_P0, gsc, gec);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_EP0, offep0, ipntbuf);    
  rd_a_gse(dom, D_EP1, offep1, ipntbuf);    
  rd_a_gse(dom, D_EP2, offep2, ipntbuf);    
  rd_a_gse(dom, D_P0, offp0, ipntbuf);    
  rd_a_gse(dom, D_P1, offp1, ipntbuf);    
  rd_a_gse(dom, D_P2, offp2, ipntbuf);    
  rd_a_gse(dom, D_MP0, offmp, ipntbuf);    
  rd_a_gse(dom, D_V0, offv0, ipntbuf);    
  rd_a_gse(dom, D_V1, offv1, ipntbuf);    
  rd_a_gse(dom, D_V2, offv2, ipntbuf);    

  // strides
  rd_a_size(dom, D_EP0, nep0);
  rd_a_size(dom, D_EP1, nep1);
  rd_a_size(dom, D_EP2, nep2);
  rd_a_size(dom, D_P0, np0);
  rd_a_size(dom, D_P1, np1);
  rd_a_size(dom, D_P2, np2);
  rd_a_size(dom, D_MP0, nmp);
  rd_a_size(dom, D_V0, nv0);
  rd_a_size(dom, D_V1, nv1);
  rd_a_size(dom, D_V2, nv2);
  
  // field update loop
  /* version 1  */
  div0 = 1;
  div1 = nv1[0];
  div2 = nv2[0]*nv2[1];

  for (i2=gsc[2];i2<gec[2]+1;i2++) {
    iep2= i2-offep2[0];
    eta2pre  = (REAL_ONE - _ep2[iep2]*dt2);
    eta2post = REAL_ONE/(REAL_ONE + _ep2[iep2]*dt2);
    ioffp0_2 = -offp0[1] + (i2-offp0[2])*np0[1];
    ioffp1_2 = -offp1[1] + (i2-offp1[2])*np1[1];
    ioffp2_2 = -offp2[1] + (i2-offp2[2])*np2[1];
    ioffmp_2 = -offmp[1] + (i2-offmp[2])*nmp[1];
    ioffv0_2 = -offv0[1] + (i2-offv0[2])*nv0[1];
    ioffv1_2 = -offv1[1] + (i2-offv1[2])*nv1[1];
    ioffv2_2 = -offv2[1] + (i2-offv2[2])*nv2[1];

    for (i1=gsc[1];i1<gec[1]+1;i1++) {
      iep1     = i1-offep1[0];
      eta1pre  = (REAL_ONE - _ep1[iep1]*dt2);
      eta1post = REAL_ONE/(REAL_ONE + _ep1[iep1]*dt2);
      ioffp0   = -offp0[0] + (i1+ioffp0_2)*np0[0];
      ioffp1   = -offp1[0] + (i1+ioffp1_2)*np1[0];
      ioffp2   = -offp2[0] + (i1+ioffp2_2)*np2[0];
      ioffmp   = -offmp[0] + (i1+ioffmp_2)*nmp[0];
      ioffv0   = -offv0[0] + (i1+ioffv0_2)*nv0[0];
      ioffv1   = -offv1[0] + (i1+ioffv1_2)*nv1[0];
      ioffv2   = -offv2[0] + (i1+ioffv2_2)*nv2[0];
      
      for (i0=gsc[0];i0<gec[0]+1;i0++) {
	iep0= i0-offep0[0];
	eta0pre = (REAL_ONE - _ep0[iep0]*dt2);
	eta0post = REAL_ONE/(REAL_ONE + _ep0[iep0]*dt2);
	ip0 = i0+ioffp0;
	ip1 = i0+ioffp1;
	ip2 = i0+ioffp2;
	imp = i0+ioffmp;
	iv0 = i0+ioffv0;
	iv1 = i0+ioffv1;
	iv2 = i0+ioffv2;
	
	sdiv = _mp[imp]*
	  (c1lam0*(_v0[iv0]-_v0[iv0-div0])+c2lam0*(_v0[iv0+div0]-_v0[iv0-2*div0]) +
	   c1lam1*(_v1[iv1]-_v1[iv1-div1])+c2lam1*(_v1[iv1+div1]-_v1[iv1-2*div1]) +
	   c1lam2*(_v2[iv2]-_v2[iv2-div2])+c2lam2*(_v2[iv2+div2]-_v2[iv2-2*div2]));

      _p0[ip0] = (_p0[ip0]*eta0pre + sdiv)*eta0post;
      _p1[ip1] = (_p1[ip1]*eta1pre + sdiv)*eta1post;
      _p2[ip2] = (_p2[ip2]*eta2pre + sdiv)*eta2post;
  
      }
    }
  }

  /* version 1 */  
  if (sgnpars->lbc[0]) {
    for (i2=0;i2<np0[2]; i2++) {
      for (i1=0; i1<np0[1]; i1++) {
	_p0[(i1+i2*np0[1])*np0[0]]=-_p0[2+(i1+i2*np0[1])*np0[0]];
	_p0[1+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
      }
    }
  }
  if (sgnpars->rbc[0]) {
    for (i2=0;i2<np0[2]; i2++) {
      for (i1=0; i1<np0[1]; i1++) {
	_p0[np0[0]-1+(i1+i2*np0[1])*np0[0]]=-_p0[np0[0]-3+(i1+i2*np0[1])*np0[0]];
	_p0[np0[0]-2+(i1+i2*np0[1])*np0[0]]= REAL_ZERO;
      }
    }
  }

  if (sgnpars->lbc[1]) {
    for (i2=0;i2<np1[2];i2++) {
      for (i0=0; i0<np1[0]; i0++) {
	// p1(i0,0,i2) = -p1(i0,2,i2)
	_p1[i0+i2*np1[0]*np1[1]]=-_p1[i0+2*np1[0]+i2*np1[0]*np1[1]];
	// p1(i0,1,i2) = 0
	_p1[i0+np1[0]+i2*np1[0]*np1[1]]= REAL_ZERO;
      }
    }
  }
  if (sgnpars->rbc[1]) {
    for (i2=0; i2<np1[0]; i2++) {
      for (i0=0; i0<np1[0]; i0++) {
	// p(i0,n1-1,i2) = -p(i0,n1-3,i2)
	_p1[i0+(np1[1]-1)*np1[0]+i2*np1[0]*np1[1]]=-_p1[i0+(np1[1]-3)*np1[0]+i2*np1[0]*np1[1]];
	// p(i0,n1-2,i2) = 0
	_p1[i0+(np1[1]-2)*np1[0]+i2*np1[0]*np1[1]]= REAL_ZERO;
      }
    }
  }

  if (sgnpars->lbc[2]) {
    for (i1=0; i1<np2[1]; i1++) {
      for (i0=0; i0<np2[0]; i0++) {
	// p(i0,i1,0) = -p(i0,i1,2)
	_p2[i0+i1*np2[0]]=-_p2[i0+i1*np2[0]+2*np2[0]*np2[1]]; 
	// p(i0,i1,1) = 0
	_p2[i0+i1*np2[0]+np2[0]*np2[1]]= REAL_ZERO;
      }
    }
  }
  if (sgnpars->rbc[2]) {
    for (i1=0; i1<np2[1]; i1++) {
      for (i0=0; i0<np2[0]; i0++) {
	// p(i0,i1,n2-1) = - p(i0,i1,n2-3)
	_p2[i0+i1*np2[0]+(np2[2]-1)*np2[1]*np2[0]]=-_p2[i0+i1*np2[0]+(np2[2]-3)*np2[1]*np2[0]];
	// p(i0,i1,n2-2) = 0
	_p1[i0+i1*np2[0]+(np2[2]-2)*np2[1]*np2[0]]= REAL_ZERO;
      }
    }
  }

  return 0;
}

int asg24_3d_v(RDOM * dom, void * pars) {

  // scaled Courant numbers
  ireal c1lam0, c2lam0;
  ireal c1lam1, c2lam1;
  ireal c1lam2, c2lam2;
  // index variables for field update loop
  int ip0, iev0, iv0, imv0, dip0;
  int ip1, iev1, iv1, imv1, dip1;
  int ip2, iev2, iv2, imv2, dip2;
  // counters
  int i0, i1, i2;
  // loop limits for computational array for pressure
  IPNT gsc0, gec0;
  IPNT gsc1, gec1;
  IPNT gsc2, gec2;
  // offsets = indices of allocated rarray origins (gs0)
  IPNT offp0, offev0, offv0, offmv0;
  IPNT offp1, offev1, offv1, offmv1;
  IPNT offp2, offev2, offv2, offmv2;
  // inner loop offsets for 2D arrays
  int ioffp0, ioffmv0, ioffv0;
  int ioffp0_2, ioffmv0_2, ioffv0_2;
  int ioffp1, ioffmv1, ioffv1;
  int ioffp1_2, ioffmv1_2, ioffv1_2;
  int ioffp2, ioffmv2, ioffv2;  
  int ioffp2_2, ioffmv2_2, ioffv2_2;  
  // strides = allocated rarray sizes
  IPNT np0, nmv0, nev0, nv0;
  IPNT np1, nmv1, nev1, nv1;
  IPNT np2, nmv2, nev2, nv2;
  // half time step
  ireal dt2;
  // PML scale factors
  ireal eta0pre, eta0post;
  ireal eta1pre, eta1post;
  ireal eta2pre, eta2post;
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
  register ireal * restrict _p2;
  register ireal * restrict _mv2;
  register ireal * restrict _ev2;
  register ireal * restrict _v2;

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
  _p2   = (dom->_s)[D_P2 ]._s0;
  _mv2  = (dom->_s)[D_MV2]._s0;
  _ev2  = (dom->_s)[D_EV2]._s0;
  _v2   = (dom->_s)[D_V2 ]._s0;

  // assign scaled Courant numbers 
  c1lam0 = C1*(sgnpars->lam[0]);
  c2lam0 = C2*(sgnpars->lam[0]);
  c1lam1 = C1*(sgnpars->lam[1]);
  c2lam1 = C2*(sgnpars->lam[1]);
  c1lam2 = C1*(sgnpars->lam[2]);
  c2lam2 = C2*(sgnpars->lam[2]);

  // half timestep
  dt2 = sgnpars->dt / 2.0;

  // size of computational domain for V1 
  rd_gse(dom, D_V0, gsc0, gec0);    
  rd_gse(dom, D_V1, gsc1, gec1);    
  rd_gse(dom, D_V2, gsc2, gec2);    
  
  // start indices of allocated domains (index offsets)
  rd_a_gse(dom, D_P0,  offp0,  ipntbuf);    
  rd_a_gse(dom, D_MV0, offmv0, ipntbuf);    
  rd_a_gse(dom, D_EV0, offev0, ipntbuf);    
  rd_a_gse(dom, D_V0,  offv0,  ipntbuf);    
  rd_a_gse(dom, D_P1,  offp1,  ipntbuf);    
  rd_a_gse(dom, D_MV1, offmv1, ipntbuf);    
  rd_a_gse(dom, D_EV1, offev1, ipntbuf);    
  rd_a_gse(dom, D_V1,  offv1,  ipntbuf);    
  rd_a_gse(dom, D_P2,  offp2,  ipntbuf);    
  rd_a_gse(dom, D_MV2, offmv2, ipntbuf);    
  rd_a_gse(dom, D_EV2, offev2, ipntbuf);    
  rd_a_gse(dom, D_V2,  offv2,  ipntbuf);    

  // strides
  rd_a_size(dom, D_P0,  np0);
  rd_a_size(dom, D_MV0, nmv0);
  rd_a_size(dom, D_EV0, nev0);
  rd_a_size(dom, D_V0,  nv0);
  rd_a_size(dom, D_P1,  np1);
  rd_a_size(dom, D_MV1, nmv1);
  rd_a_size(dom, D_EV1, nev1);
  rd_a_size(dom, D_V1,  nv1);
  rd_a_size(dom, D_P2,  np2);
  rd_a_size(dom, D_MV2, nmv2);
  rd_a_size(dom, D_EV2, nev2);
  rd_a_size(dom, D_V2,  nv2);

  dip0 = 1;
  dip1 = np1[0];
  dip2 = np2[0]*np2[1];

  for (i2=gsc0[2]; i1<gec0[2]+1; i2++) {
    ioffp0_2  = -offp0[1] + (i2-offp0[2]) *np0[1];
    ioffmv0_2 = -offmv0[1] + (i2-offmv0[2])*nmv0[1];
    ioffv0_2  = -offv0[1] + (i2-offv0[2]) *nv0[1];
    for (i1=gsc0[1]; i1<gec0[1]+1; i1++) {
      ioffp0  = -offp0[0]  + (i1+ioffp0_2) *np0[0];
      ioffmv0 = -offmv0[0] + (i1+ioffmv0_2)*nmv0[0];
      ioffv0  = -offv0[0]  + (i1+ioffv0_2) *nv0[0];
      for (i0=gsc0[0]; i0<gec0[0]+1; i0++) {
	iev0 = i0-offev0[0];
	eta0pre = (REAL_ONE - _ev0[iev0]*dt2);
	eta0post = REAL_ONE/(REAL_ONE + _ev0[iev0]*dt2);      
	ip0  = i0+ioffp0;
	imv0 = i0+ioffmv0;
	iv0  = i0+ioffv0;
	
	_v0[iv0] += (eta0pre*eta0post - REAL_ONE)*_v0[iv0] + eta0post * 
	  _mv0[imv0]*(c1lam0*(_p0[ip0+dip0]   - _p0[ip0])+
		      c2lam0*(_p0[ip0+2*dip0] - _p0[ip0-dip0]));
      }
    }
  }

  for (i2=gsc1[2]; i2<gec1[2]+1; i2++) {
    ioffp1_2  = -offp1[1]  + (i2 - offp1[2]) *np1[1];
    ioffmv1_2 = -offmv1[1] + (i2 - offmv1[2])*nmv1[1];
    ioffv1_2  = -offv1[1]  + (i2 - offv1[2]) *nv1[1];
    for (i1=gsc1[1]; i1<gec1[1]+1; i1++) {
      iev1     = i1-offev1[0];
      eta1pre  = (REAL_ONE - _ev1[iev1]*dt2);
      eta1post = REAL_ONE/(REAL_ONE + _ev1[iev1]*dt2);      
      ioffp1   = -offp1[0]  + (i1 + ioffp1_2) *np1[0];
      ioffmv1  = -offmv1[0] + (i1 + ioffmv1_2)*nmv1[0];
      ioffv1   = -offv1[0]  + (i1 + ioffv1_2) *nv1[0];
      for (i0=gsc1[0]; i0<gec1[0]+1; i0++) {
	ip1  = i0+ioffp1;
	imv1 = i0+ioffmv1;
	iv1  = i0+ioffv1;
	_v1[iv1] += (eta1pre*eta1post-REAL_ONE)*_v1[iv1] + eta1post*
	  _mv1[imv1]*(c1lam1*(_p1[ip1+dip1]-_p1[ip1])+
		      c2lam1*(_p1[ip1+2*dip1]-_p1[ip1-dip1]));
      }
    }
  }

  for (i2=gsc2[2]; i2<gec2[2]+1; i2++) {
    iev2      = i2-offev2[0];
    eta2pre   = (REAL_ONE - _ev2[iev2]*dt2);
    eta2post  = REAL_ONE/(REAL_ONE + _ev2[iev2]*dt2);      
    ioffp2_2  = -offp2[1]  + (i2 - offp2[2]) *np2[1];
    ioffmv2_2 = -offmv2[1] + (i2 - offmv2[2])*nmv2[1];
    ioffv2_2  = -offv2[1]  + (i2 - offv2[2]) *nv2[1];
    for (i1=gsc2[1]; i1<gec2[1]+1; i1++) {
      ioffp2  = -offp2[0]  + (i1 + ioffp2_2) *np2[0];
      ioffmv2 = -offmv2[0] + (i1 + ioffmv2_2)*nmv2[0];
      ioffv2  = -offv2[0]  + (i1 + ioffv2_2) *nv2[0];
      for (i0=gsc2[0]; i0<gec2[0]+1; i0++) {
	ip2  = i0+ioffp2;
	imv2 = i0+ioffmv2;
	iv2  = i0+ioffv2;
	_v2[iv2] += (eta2pre*eta2post-REAL_ONE)*_v2[iv2] + eta2post*
	  _mv2[imv2]*(c1lam2*(_p2[ip2+dip2]-_p2[ip2])+
		      c2lam2*(_p2[ip2+2*dip2]-_p2[ip2-dip2]));
      }
    }
  }

  if (sgnpars->lbc[0]) {
    for (i2=0; i2<nv0[2]; i2++) {
      for (i1=0; i1<nv0[1]; i1++) {  
	// v(0,i1,i2)=v(1,i1,i2)
	_v0[(i1+i2*nv0[1])*nv0[0]]=_v0[1+(i1+i2*nv0[1])*nv0[0]];
      }
    }
  }

  if (sgnpars->rbc[0]) {
    for (i2=0; i2<nv0[2]; i2++) {
      for (i1=0; i1<nv0[1]; i1++) {    
	// v(n0-1,i1,i2) = v(n0-2,i1,i2)
	_v0[nv0[0]-1+(i1+i2*nv0[1])*nv0[0]]=_v0[nv0[0]-2+(i1+i2*nv0[1])*nv0[0]];
      }
    }
  }

  if (sgnpars->lbc[1]) {
    for (i2=0; i2<nv1[2]; i2++) {
      for (i0=0; i0<nv1[0]; i0++) {    
	// v(i0,0,i2)=v(i0,1,i2);
	_v1[i0+i2*nv1[1]*nv1[0]]=_v1[i0+nv1[0]+i2*nv1[1]*nv1[0]];
      }
    }
  }
  if (sgnpars->rbc[1]) {
    for (i2=0; i2<nv1[2]; i2++) {
      for (i0=0; i0<nv1[0]; i0++) {    
	// v(i0,n1-1,i2) = v(i0,n1-2,i2)
	_v1[i0+(nv1[1]-1)*nv1[0]+i2*nv1[1]*nv1[0]]=_v1[i0+(nv1[1]-2)*nv1[0]+i2*nv1[1]*nv1[0]];
      }
    }
  }

  if (sgnpars->lbc[2]) {
    for (i1=0;i1<nv2[1];i1++) {
      for (i0=0;i0<nv2[0];i0++) {
	// v(i0,i1,0) = v(i0,i1,1)
	_v2[i0+i1*nv2[0]] = _v2[i0+i1*nv2[0]+nv2[1]*nv2[0]];
      }
    }
  }

  if (sgnpars->rbc[2]) {
    for (i1=0;i1<nv2[1];i1++) {
      for (i0=0;i0<nv2[0];i0++) {
	// v(i0,i1,n2-1) = v(i0,i1,n2-2)
	_v2[i0+i1*nv2[0]+(nv2[2]-1)*nv2[1]*nv2[0]] = _v2[i0+i1*nv2[0]+(nv2[2]-2)*nv2[1]*nv2[0]];
      }
    }
  }

  return 0;
}
        
int asg24_3d(RDOM *dom, int iv, void *pars) {

  if ( iv == 0 ) return asg24_3d_p(dom, pars);
  if ( iv == 1 ) return asg24_3d_v(dom, pars);

  return 0;
}

